#include "pdp1.h"

#define BC1    1
#define BC2    2
#define BC3    3
#define BC4    4

SCALAR source(void), bc1(void), bc2(void), bc3(void), bc4(void), (*bcptr)(void);
void  circuit1(void), circuit2(void), circuit3(void), circuit4(void), (*circuitptr)(void);
SCALAR *gam, *tri_a, *tri_b, *tri_c, se, jp;
int   bc_flag;
SCALAR cs;
void field_init(void);

/*************************************************/

void fields(void)
{
  static int init_flag=1;
  register int j, isp;
  SCALAR s, bet;

  if(init_flag) {
    field_init();
    init_flag = 0;
  }
	
	/* save old field on the end for injection */
	eold[0]=e[0];
	eold[1]=e[nc];

  /********************************************************/
  /* calculate total rho(x).  If implicit, also calculate */
  /* the numerical (implicit) susceptibility chi.  If     */
  /* explicit, sp_chi_scale=0 for all species => chi=0.   */

  for (j=0; j<=nc; j++) {
    rho[j] = rhoback;
    chi[j] = 0.0;
  }
  
  for (isp=0; isp<nsp; isp++)
    for (j=0; j<=nc; j++) {
      chi[j] += sp_n[isp][j]*sp_chi_scale[isp];
      rho[j] += sp_n[isp][j]*q[isp]/dx;
    }
  
  /**********************************************************/
  /* recalculate the coeffs. for the tridiagnol field solve */
  
  for (j=1; j<nc; j++) {
    tri_a[j] = 1.0 + 0.5*(chi[j] + chi[j-1]);
    tri_c[j] = 1.0 + 0.5*(chi[j] + chi[j+1]);
    tri_b[j] = -(tri_a[j] + tri_c[j]);
  }
  if(bc_flag == BC1 || bc_flag == BC2) {
    tri_b[0] = -tri_a[1];
    tri_c[0] =  tri_a[1];
  }
  else if(bc_flag == BC4) {
    tri_b[0] = 2.25*extl/dt/dt + 1.5*extr/dt + 1/extc;
    tri_b[0] = -(tri_a[1] + dx/tri_b[0]/area/epsilon);
    tri_c[0] =  tri_a[1];
  }

  /*******************************************/
  /* solve a tridiagnal matrix to get Phi[i] */
  
  jp = backj;

	for (isp=0; isp<nsp; isp++){
		jp += jwall[isp];
		// printf("%d \t %lf \n",isp,jwall[isp]);
	}

//	printf("total j %lf \n",jp);

	cs = source();
	
  if (nstrt) phi[0] = cs;
  bet = tri_b[nstrt];
  phi[nstrt] = (*bcptr)()/bet;
  
  for (j=nstrt+1; j<nc; j++) {
    gam[j] = tri_c[j-1]/bet;
    bet = tri_b[j] - tri_a[j]*gam[j];
    phi[j] = -(se*rho[j] + tri_a[j]*phi[j-1])/bet;
  }
  phi[nc] = 0.;
  for(j=nc-2; j>= nstrt; j--) phi[j] -= gam[j+1]*phi[j+1];
  
  /*******************************/
  /* calculate E(x) from phi(x)  */
  
  s = 0.5/dx;
  for (j=1; j<nc; j++) e[j] = (phi[j-1] - phi[j+1])*s;

	e[0]  = (1.0+0.5*(chi[0]+chi[1]))*(phi[0] -phi[1])/dx  -rho[0]*0.5*dx/epsilon;
	e[nc] = (1.0+0.5*(chi[nc-1]+chi[nc]))*(phi[nc-1] -phi[nc])/dx +rho[nc]*0.5*dx/epsilon;

  /****************************************/
  /* calculate the external circuit stuff */
  
  (*circuitptr)();
}

/***************************************************************/

SCALAR source(void)
{
  if (dcramped) {
    if(t < risetime) {
      if(ramp) return (ramp*t);
      else     return (.5*dcbias*(1- cos(t*w0)));
    }
    else       return (dcbias);
  }
  return (dcbias + ramp*t + acbias*sin(t*w0 + theta0));
}

/***************************************************************/
/* Initalizing the arrays and parameters for the field solve   */

void field_init(void)
{     
  /* Setting up the arrays for the poisson solve */
  gam  = (SCALAR *)malloc(ng*sizeof(SCALAR));
  tri_a= (SCALAR *)malloc(ng*sizeof(SCALAR));
  tri_b= (SCALAR *)malloc(ng*sizeof(SCALAR));
  tri_c= (SCALAR *)malloc(ng*sizeof(SCALAR));
  
  /******************************************/
  /* Deciding which circuit solver to use.. */
  if(extc < 1e-30) {   /* when C -> 0.0, OPEN CIRCUIT */
    if(fabs(dcbias)+fabs(ramp)+fabs(acbias) > 0)
      puts("START: Active external source is ignored since C < 1E-30\n");
    bc_flag = BC1;
    bcptr = bc1;
    circuitptr = circuit1;
    nstrt = 0;
  }
  else if(src== 'I' || src== 'i') {  /* When current source is applied */
    bc_flag = BC2;
    bcptr = bc2;
    circuitptr = circuit2;
    nstrt = 0;
  }
  else if(extl < 1e-30 && extr < 1e-30 && extc >= 1e5*epsilon*area/length) {
    /* When L=R=0.0 and C -> infinity, SHORT CIRCUIT */
    bc_flag = BC3;
    bcptr = bc3;
    circuitptr = circuit3;
    nstrt = 1;
  }
  else {             /* The general case with external voltage source */
    bc_flag = BC4;
    bcptr = bc4;
    circuitptr = circuit4;
    nstrt = 0;
  }
  se = dx*dx/epsilon;   /* setting se = dx*dx/epsilon */
}

/***************************************************************/
/* When C = 0: OPEN CIRCUIT */

SCALAR bc1()
{
  return (-(sigma + dt*jp + 0.5*rho[0]*dx)*dx/epsilon);
}

void circuit1()
{
  oldsigma= sigma;
  sigma  += jp*dt;
}

/***************************************************************/
/* When CURRENT SOURCE is applied */

SCALAR bc2()
{
  return (-(sigma + dt*(cs/area + jp) + 0.5*rho[0]*dx)*dx/epsilon);
}

void circuit2()
{
  exti    = cs;
  oldsigma= sigma;
  sigma  += dt*(jp + exti/area);
}

/***************************************************************/
/* When R=L=0, C -> infinity: SHORT CIRCUIT */

SCALAR bc3()
{
  return (-rho[1]*dx*dx/epsilon - tri_a[1]*cs);
}

void circuit3()
{
  oldsigma = sigma;
  sigma= epsilon*e[0];
	if (t==0)
		oldsigma=sigma;
  exti = area*((sigma - oldsigma)/dt- jp);
}

/***************************************************************/
/* The General case */

SCALAR bc4()
{
  static SCALAR a0, a1, a2, a3, a4;
  SCALAR k;
  
  if(!a0) {
    a0 = 2.25*extl/dt/dt + 1.5*extr/dt + 1/extc;
    a1 = -6*extl/dt/dt - 2*extr/dt;
    a2 = 5.5*extl/dt/dt + .5*extr/dt;
    a3 = -2*extl/dt/dt;
    a4 = .25*extl/dt/dt;
  }
  k = (a1*extq + a2*extq_1 + a3*extq_2 + a4*extq_3)/a0;
  return (-(0.5*rho[0]*dx + sigma + dt*jp
	    + (cs/a0 - k - extq)/area)*dx/epsilon);
}

void circuit4() {
  static SCALAR a0, a1, a2, a3, a4;
  SCALAR k;
  
  if(!a0) {
    a0 = 2.25*extl/dt/dt + 1.5*extr/dt + 1/extc;
    a1 = -6*extl/dt/dt - 2*extr/dt;
    a2 = 5.5*extl/dt/dt + .5*extr/dt;
    a3 = -2*extl/dt/dt;
    a4 = .25*extl/dt/dt;
  }
  k = (a1*extq + a2*extq_1 + a3*extq_2 + a4*extq_3)/a0;
  extq_3 = extq_2;
  extq_2 = extq_1;
  extq_1 = extq;
  extq = (cs - phi[0])/a0 - k;
  
  exti = (extq - extq_1)/dt;
  oldsigma= sigma;
  sigma  += (jp + exti/area)*dt;
}

/***************************************************************/
