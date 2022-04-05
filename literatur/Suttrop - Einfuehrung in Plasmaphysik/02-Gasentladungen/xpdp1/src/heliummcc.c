#include "pdp1.h"
#include "xgrafix.h"

#define   NEMAX    500
#define   NIMAX    500

void hnewvel(SCALAR, SCALAR, SCALAR *, SCALAR *, SCALAR *, int);
SCALAR hsigma1(SCALAR), hsigma2(SCALAR), hsigma3(SCALAR), hsigma4(SCALAR), hsigma5(SCALAR);
SCALAR gden, vgth, extengy0, ionengy0;
int ecolsp, icolsp, ionsp;

/*********************************************************/
/* Monte Carlo Collisions for electron- and ion-neutrals */

void heliummcc(int isp)
{
  static int init_flag=1;
  static SCALAR ecol_extra, icol_extra;
  static SCALAR max_sigmav_e, max_sigmav_i;
  static SCALAR col_prob_e, col_prob_i;

  register int j, k, index;
  int N, nnp, s;
  SCALAR random, dum, vel, sigma_total, phi1, cosphi, sinphi, coschi, sinchi;
  SCALAR engy, rengy, del;
  SCALAR sum_sigma, temp, vneutx, vneuty, vneutz;
  
  if (init_flag) {
    gden = NperTORR*pressure/(gtemp+DBL_MIN);  /* Calculating the gas density */
    
    ecolsp= ecollisional-1;
    icolsp= icollisional-1;
    ionsp = ionspecies-1;     /* Fixing the indices into the array of species */
    
/* mindgame: This is needed for both ecollisional and icollisional cases. */
    if(ionsp>=0) vgth= sqrt(gtemp/Escale[ionsp]);
    /*****************************************/
    /* Calculating the null collision prob.  */
    
    if(ecollisional) {
      
      extengy0 = 19.8;
      ionengy0 = 24.5;
      
      max_sigmav_e = 0.0;
      for (engy=0; engy<NEMAX; engy += 0.1)
				max_sigmav_e= max(max_sigmav_e,
													sqrt(2*1.602e-19*engy/m[ecolsp])
													*(hsigma1(engy)+hsigma2(engy)+hsigma3(engy)));
      
      col_prob_e = 1 -exp(-max_sigmav_e*gden*dt*sp_k[ecolsp]);
      ecol_extra = 0.5;
    }
    if(icollisional) {
      max_sigmav_i = 0.0;
      for (engy=0; engy<NIMAX; engy += 0.1)
	max_sigmav_i= max(max_sigmav_i, sqrt(2*1.602e-19*engy/m[icolsp])*(hsigma4(engy)+hsigma5(engy)));
      
      col_prob_i = 1 -exp(-max_sigmav_i*gden*dt*sp_k[icolsp]);
      icol_extra = 0.5;
    }
    init_flag = 0;
  }
  
  /********************************/
  /* Electron collisions with He  */

  if(ecollisional && isp==ecolsp) {
    /**** First Clear the diagnostics  ******/
    if(theRunWithXFlag) {
      for (j=0; j<ng; j++)
				mccrate[0][j] = mccrate[1][j] = mccrate[2][j] = 0.0;
    }
    /**** Now do the collisions *****/
    ecol_extra += np[ecolsp]*col_prob_e;
    N = ecol_extra;
    ecol_extra -= N;
    
    nnp = np[ecolsp];
    for(j=0; j< N; j++) {
      index= nnp*frand();
      nnp--;
      temp = x[ecolsp][nnp];
      x[ecolsp][nnp] = x[ecolsp][index];
      x[ecolsp][index] = temp;
      
      temp = vx[ecolsp][nnp];
      vx[ecolsp][nnp] = vx[ecolsp][index];
      vx[ecolsp][index] = temp;
      
      temp = vy[ecolsp][nnp];
      vy[ecolsp][nnp] = vy[ecolsp][index];
      vy[ecolsp][index] = temp;
      
      temp = vz[ecolsp][nnp];
      vz[ecolsp][nnp] = vz[ecolsp][index];
      vz[ecolsp][index] = temp;
    }
    
    for(j=nnp; j<nnp+N; j++) {
      dum = (vx[ecolsp][j]*vx[ecolsp][j] +vy[ecolsp][j]*vy[ecolsp][j]
						 +vz[ecolsp][j]*vz[ecolsp][j]);
      engy= Escale[ecolsp]*dum;
      vel = sqrt(dum);
      sigma_total = max_sigmav_e/(vel*vscale);
      random= frand();
      
      /*********************************************************************/
      /* determine the type of collision and calculate new velocities, etc */
      
      /*******************************/
      /* if the collision is elastic */
      
      if (random <= (sum_sigma =hsigma1(engy))/sigma_total) { 
				/* first normalize vel */
				vx[ecolsp][j] /= vel;
				vy[ecolsp][j] /= vel;
				vz[ecolsp][j] /= vel;
	
				/* scatter the electron */
				hnewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 1);
	
				/* collect diagnostics */
				if(theRunWithXFlag) {
					s= x[ecolsp][j];
					del= x[ecolsp][j] - s;
					mccrate[0][s]  += (!s) ? 2*(1-del) : 1-del;
					mccrate[0][s+1]+= (s== nc-1) ? 2*del : del;
				}
      }
      /**********************************/
      /* if the collision is excitation */
      
      else if (engy>= extengy0 && random <= (sum_sigma +=hsigma2(engy))/sigma_total) {
				/* first normalize vel */
				vx[ecolsp][j] /= vel;
				vy[ecolsp][j] /= vel;
				vz[ecolsp][j] /= vel;
	
				engy -= extengy0;
				vel   = sqrt(fabs(engy)/Escale[ecolsp]);
	
				/* scatter the electron */
				hnewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
				/* collect diagnostics */
				if(theRunWithXFlag) {
					s= x[ecolsp][j];
					del= x[ecolsp][j] - s;
					mccrate[1][s]  += (!s) ? 2*(1-del) : 1-del;
					mccrate[1][s+1]+= (s== nc-1) ? 2*del : del;
				}
      } 
      /***************************************************************************/
      /* if the collision is ionization, add the released electron and ion first */
      
      else if(engy >= ionengy0 && random <= (sum_sigma + hsigma3(engy))/sigma_total) {
				/* first normalize vel */
				vx[ecolsp][j] /= vel;
				vy[ecolsp][j] /= vel;
				vz[ecolsp][j] /= vel;
	
				/********************************/
				/* subtract the ion. energy and */
				/* partition the remaining energy */
	
				engy -= ionengy0;
				rengy = 10.0*tan(frand()*atan(engy/20.0));
				engy -= rengy;
	
				/********************************/
				/* scatter the created electron */
	
				vel = sqrt(fabs(rengy)/Escale[ecolsp]);
				k = np[ecolsp];
				vx[ecolsp][k] = vx[ecolsp][j];
				vy[ecolsp][k] = vy[ecolsp][j];
				vz[ecolsp][k] = vz[ecolsp][j];
				hnewvel(rengy, vel, &vx[ecolsp][k], &vy[ecolsp][k], &vz[ecolsp][k], 0);
				x[ecolsp][k] = x[ecolsp][j];
	
				/****************************************/
				/* assign velocities to the created ion */
	
				k = np[ionsp];
				maxwellv(&vx[ionsp][k], &vy[ionsp][k], &vz[ionsp][k], vgth);
				x[ionsp][k] = x[ecolsp][j];
	
				s = x[ionsp][k];
				del = x[ionsp][k] - s;
				sp_n_mcc[ionsp][s]  += 1. - del;
				sp_n_mcc[ionsp][s+1]+= del;
	
				/*****************************************/
				/* finally scatter the incident electron */
	
				vel = sqrt(fabs(engy)/Escale[ecolsp]);
				hnewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
				if(++np[ionsp] >maxnp[ionsp] || ++np[ecolsp] >maxnp[ecolsp]) {
					puts("ADJUST(Ionization): too many particles. MUST EXIT!");
					exit(-1);
				}
	
				/* collect diagnostics */
				if(theRunWithXFlag) {
					mccrate[2][s]  += (!s) ? 2*(1-del) : 1-del;
					mccrate[2][s+1]+= (s== nc-1) ? 2*del : del;
				}
      }
    } 
  }
  
  /**************************************/
  /* He+ + He -> .....                  */
  
  if(icollisional && isp==icolsp) {
    /**** First Clear the diagnostics  ******/
    if(theRunWithXFlag) {
      for (j=0; j<ng; j++) mccrate[3][j] = mccrate[4][j] = 0.0;
    }
    /**** Now do the collisions *****/
    icol_extra += np[icolsp]*col_prob_i;
    N = icol_extra;
    icol_extra -= N;
    
    nnp = np[icolsp];
    for(j=0; j< N; j++) {
      index= nnp*frand();
      nnp--;
      temp = x[icolsp][nnp];
      x[icolsp][nnp] = x[icolsp][index];
      x[icolsp][index] = temp;
      
      temp = vx[icolsp][nnp];
      vx[icolsp][nnp] = vx[icolsp][index];
      vx[icolsp][index] = temp;
      
      temp = vy[icolsp][nnp];
      vy[icolsp][nnp] = vy[icolsp][index];
      vy[icolsp][index] = temp;
      
      temp = vz[icolsp][nnp];
      vz[icolsp][nnp] = vz[icolsp][index];
      vz[icolsp][index] = temp;
    }
    
    for(j=nnp; j<nnp+N; j++) {
      maxwellv(&vneutx, &vneuty, &vneutz, vgth);
      vx[icolsp][j] -= vneutx;
      vy[icolsp][j] -= vneuty;
      vz[icolsp][j] -= vneutz;
      dum = (vx[icolsp][j]*vx[icolsp][j] +vy[icolsp][j]*vy[icolsp][j]
	     +vz[icolsp][j]*vz[icolsp][j]);
      if (dum){
	engy= Escale[icolsp]*dum;
	vel = sqrt(dum);
	sigma_total = max_sigmav_i/(vel*vscale);
	random= frand();
      
      /***************************************/
      /* if the collision is charge exchange */
   
      if (random <= (sum_sigma =hsigma5(engy))/sigma_total) {
	
	newivel(engy, vel, &vx[icolsp][j], &vy[icolsp][j], &vz[icolsp][j]);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[icolsp][j];
	  del= x[icolsp][j] - s;
	  mccrate[3][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[3][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      
      else if (random <= (sum_sigma +=hsigma4(engy))/sigma_total) {
	vx[icolsp][j] = vy[icolsp][j] = vz[icolsp][j] = 0.0;
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[icolsp][j];
	  del= x[icolsp][j] - s;
	  mccrate[4][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[4][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      }
      vx[icolsp][j] += vneutx;
      vy[icolsp][j] += vneuty;
      vz[icolsp][j] += vneutz;
    }
  }
}
/**************************************************************/

void hnewvel(SCALAR energy, SCALAR vel, SCALAR *vx, SCALAR *vy, SCALAR *vz, int e_flag)
{
  SCALAR phi1, cosphi, sinphi, coschi, sinchi, up1, up2, up3;
  SCALAR mag, r11, r12, r13, r21, r22, r23, r31, r32, r33;
  
  /* if(energy < 1e-30)  coschi = 1;
     else  coschi = (energy +2 -2*pow(energy +1,frand()))/energy; */
  
  coschi = 1 -2.0*frand();
  sinchi= sqrt(1. -fabs(coschi*coschi));
    
  phi1  = 2*M_PI*frand();
  cosphi= cos(phi1);
  sinphi= sin(phi1);
  
  if(e_flag)  vel *= sqrt(1 - 2*m[ecolsp]*(1-coschi)/m[ionsp]);
  
  r13 = *vx;
  r23 = *vy;
  r33 = *vz;
  
  if(r33 == 1.0) { up1= 0;  up2= 1;  up3= 0; }
  else           { up1= 0;  up2= 0;  up3= 1; }
  
  r12 = r23*up3 -r33*up2;
  r22 = r33*up1 -r13*up3;
  r32 = r13*up2 -r23*up1;
  mag = sqrt(r12*r12 + r22*r22 + r32*r32);
  
  r12/= mag;
  r22/= mag;
  r32/= mag;
  
  r11 = r22*r33 -r32*r23;
  r21 = r32*r13 -r12*r33;
  r31 = r12*r23 -r22*r13;
  
  *vx= vel*(r11*sinchi*cosphi +r12*sinchi*sinphi +r13*coschi);
  *vy= vel*(r21*sinchi*cosphi +r22*sinchi*sinphi +r23*coschi); 
  *vz= vel*(r31*sinchi*cosphi +r32*sinchi*sinphi +r33*coschi);
}

/***********************************/
/*  e + He -> e + He  Elastic      */

SCALAR hsigma1(SCALAR energy)
{
  return (8.5e-19/(pow(energy+10.0, 1.1)));
}

/************************************/
/*  e + He -> e + He  Excitation    */

SCALAR hsigma2(SCALAR energy)
{
  if(energy < extengy0) return (0.0);
  else if(extengy0 <= energy && energy <27.0)
    return (2.08e-22*(energy -extengy0));
  return (3.4e-19/(energy +200));
}

/************************************/
/*  e + He -> e + e + He+  Ion.     */

SCALAR hsigma3(SCALAR energy)
{
  if(energy < ionengy0) return (0.0);
	return(1e-17*(energy -ionengy0)/((energy +50)*pow(energy+300.0, 1.2)));
}

/************************************/
/*  He + He+ -> He+ + He  Charge X  */

SCALAR hsigma4(SCALAR energy)
{
  SCALAR crossx;
  const SCALAR emin = 0.01;
  const SCALAR cutoff = 377.8;

  if(energy < emin) { energy = emin; }
  if(energy < cutoff)
    {
      crossx = 1.2996e-19  - 7.8872e-23*energy + 1.9873e-19/sqrt(energy);
    }
  else
    {
      crossx = (1.5554e-18*log(energy)/energy) + 8.6407e-20;
    }
  return crossx;
 
}

/************************************/
/*  He + He+ -> He + He+  Elastic  */

SCALAR hsigma5(SCALAR energy)
{
  SCALAR crossx;
  const SCALAR emin = 0.01;
  if(energy < emin) { energy = emin; }

  crossx = 3.6463e-19/sqrt(energy) - 7.9897e-21;
  //  fprintf(stderr, "\nenergy=%lf, crossx=%lf", energy, crossx);
  if(crossx < 0) { return 0.0; }
  return crossx;
}

/**************************************************************/

void hmakethefile(void)
{
  SCALAR e;
  FILE *DMPFile;
  
  DMPFile = fopen("xsections", "w"); 
  
  for(e=0; e<300; e+=.1)
    fprintf(DMPFile, "%lf %lf %lf %lf %lf\n", e+1e-30,
            hsigma1(e)+1e-30, hsigma2(e)+1e-30, hsigma3(e)+1e-30, hsigma4(e)+1e-30);
  
  fclose(DMPFile);
}

/**************************************************************/


