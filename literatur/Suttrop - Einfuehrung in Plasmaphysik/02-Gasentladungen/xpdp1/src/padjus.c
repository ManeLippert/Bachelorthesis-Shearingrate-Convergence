#include "pdp1.h"
#include "xgrafix.h"

/****************************************************************/
/* Routine to adjust (initialize, re-pack and inject) particles */
/*  to the desired boundary conditions                           */

void injection_push(int species, int particle, SCALAR part_time);
void injection_push_oldpdp(int species, int particle, SCALAR part_time);
void injection_push_FTSBP(int species, int particle, SCALAR part_time);
void injection_push_MBP(int species, int particle, SCALAR part_time);
void injection_push_FG(int species, int particle, SCALAR part_time);

void sterm(int number);
void ionization(int number);

void adjust(int isp)
{
  static int ionsp, npold;
  static int secsp, init_flag=1;
  static SCALAR extra[NSMAX][2];
  static SCALAR this_time;
  
  register int i, ii, j, k;
  int nnp, secountl=0, secountr=0;
  int nreflux[NSMAX], s;
  SCALAR dum, del_t=0;
  
  /* INITIALIZE array for computing positions of injected particles */
  if (init_flag) {

    ionsp = ionspecies-1;
    secsp = secondary -1;
    
    /* "enter" is now the no. of particles injected each "dt" */
    for (i=0; i<nsp; i++) {
      if (fabs(enter[i][0]) > 0.) extra[i][0] = 0.5123123;
      else extra[i][0] = 0.;
      if (fabs(enter[i][1]) > 0.) extra[i][1] = 0.5123123;
      else extra[i][1] = 0.;
    }
    init_flag = 0;
  }
  
  if (psource && (isp == ionsp))
    npold = np[ionsp]; /* Save the number of ions before adjust */
	
  if (isp==0)  // secondaries add to jwall[secsp] not jwall[isp]
		 for (i=0; i<nsp; i++)
		   jwall[i]=0;
	
  if(secondary) secountl = secountr = 0;
  nreflux[isp]= 0;    /* nreflux[isp] equals zero */
  
  if (np[isp] > 0) {
    nnp = np[isp] -1;
    i = 0;
    
    /* eliminate "outsiders", allow for secondary electron emission, */
    /* and if it left thru LH plate, ADD charge there to sigma       */
    /* (plate surface density). */

    /* note: the system starts at 0 (includes 0) and end before xnc   */
    /******************** DOES NOT INCLUDE xnc ************************/


    do 
      { 
	if (x[isp][i] >= xnc) 
	  {
	    if(!wraparound)
	      {
		x[isp][i] = x[isp][nnp];
		vx[isp][i] = vx[isp][nnp];
		vy[isp][i] = vy[isp][nnp];
		vz[isp][i] = vz[isp][nnp];
		nnp--;
		if(reflux)  nreflux[isp]++;
		else if(secondary) 
		  {
		    if (frand() < seec[isp]) 
		      {
			secountr++;
		      }
		  }
	      }
	    else
	      {
		x[isp][i] -= xnc;
	      }
	  }
	else if (x[isp][i] < 0) 
	  {
	    if(!wraparound) 
	      {
		/********** LHS wall diagnostics   ****************/
		if(theRunWithXFlag) 
		  {
		    dum = (vx[isp][i]*vx[isp][i] + vy[isp][i]*vy[isp][i]
			   +vz[isp][i]*vz[isp][i] - emin[isp])/de[isp];
		    s = dum;
		    if (s<nbin[isp]-1 && dum>=0) 
		      {
			dum -= s;
			fe[isp][s]  += (!s) ? 2*(1-dum) : 1-dum;
			fe[isp][s+1]+= (s==nbin[isp]-2) ? 2*dum : dum;
		      }
		    dum = -atan(sqrt(vy[isp][i]*vy[isp][i]+vz[isp][i]*vz[isp][i])
				/vx[isp][i])/dtheta[isp];
		    s = dum;
		    dum -= s;
		    ftheta[isp][s]  += (!s) ? 2*(1-dum) : 1-dum;
		    ftheta[isp][s+1]+= (s==nbin[isp]-2) ? 2*dum : dum;
		  }
		x[isp][i]  = x[isp][nnp];
		vx[isp][i] = vx[isp][nnp];
		vy[isp][i] = vy[isp][nnp];
		vz[isp][i] = vz[isp][nnp];
		nnp--;
		jwall[isp] += jnorm[isp];
		if (secondary) 
		  {
		    if (frand() < seec[isp]) 
		      {
			secountl++;
			jwall[secsp] -= jnorm[secsp];
		      }
		  }
	      }
	    else
	      {
		x[isp][i] += xnc;
	      }
	  }
	else {
	/******  MID system diagnostics  ************/
	if(theRunWithXFlag) {
	  if (xs_mid[isp] <= x[isp][i] && x[isp][i] <= xf_mid[isp]) {
	    s = (vx[isp][i]*vx[isp][i] + vy[isp][i]*vy[isp][i]
		 +vz[isp][i]*vz[isp][i] - emin_mid[isp])/de_mid[isp];
						
	    if (0 <= s && s< nbin_mid[isp]) {
	      fe_mid[isp][s] += 1;
	      sp_fe[isp][s] += 1;
	    }
	  }
	}
	i++;
      }
    } while (i <= nnp);
    np[isp] = nnp + 1;
  }
  
  /********************************************************/
  /* INJECT new particles at walls, one species at a time */
  
  if (inject[isp] || nreflux[isp]) {
    for(k=0; k<2; k++) {
      extra[isp][k] += enter[isp][k];
      if (enter[isp][k])
	del_t = 1/enter[isp][k];
      if (reflux && k == 1) {
	extra[isp][k] += nreflux[isp];
	if ((enter[isp][k]+nreflux[isp])!=0.0)
	  del_t = 1/(enter[isp][k] + nreflux[isp]);
      }
      while (extra[isp][k] >= 1.0) {
	extra[isp][k] -= 1.0;
	ii = np[isp];
	np[isp]++;
				
	if (ii >= maxnp[isp]) {    /* Move array boundaries here */
	  printf("ADJUST: too many particles, species %d",isp);
	  exit(1);
	}
				
	/* Choose V's */
	vx[isp][ii] = distribution_flux(k,isp,vxloader[isp][k]);
	vy[isp][ii] = v0y[isp]+vty[isp]*maxwellian(vyloader[isp]);
	vz[isp][ii] = v0z[isp]+vtz[isp]*maxwellian(vzloader[isp]);
	if (k)
	  x[isp][ii] = xnc-FLOAT_MIN;
	else
	  x[isp][ii] = 0;
	/* Adjust Vx,x for effect of E and B field 
	   for a partial timestep push*/
		injection_push(isp,ii,extra[isp][k]*del_t); 
/*	 injection_push_oldpdp(isp,ii,extra[isp][k]*del_t); */
	if  ((x[isp][ii]>=xnc)||(x[isp][ii]<0)){
	  np[isp]--;
	}
				
	if (!k) jwall[isp] -= jnorm[isp];
      }
    } 
  } 
  
  if(secondary) {
    i = np[secsp];
    //    np[secsp] += secountl +secountr;
    
    //    if(np[secsp] >= maxnp[isp])
      //      printf("ADJUST(Secondaries): too many particles. MUST EXIT!");
    
    for(j=i; j< i+secountl; j++) {
      x[secsp][j]= 0;
      vx[secsp][j]= distribution_flux(0,secsp,vxloader[secsp][0]);
      vy[secsp][j] = v0y[secsp]+vty[secsp]*maxwellian(vyloader[secsp]);
      vz[secsp][j] = v0z[secsp]+vtz[secsp]*maxwellian(vzloader[secsp]);
      injection_push(secsp,j,frand());
      if  (x[secsp][j]<0){
	j--;
	secountl--;
      }
    }
    i += secountl;
    for(j=i; j< i+secountr; j++) {
      x[secsp][j]= xnc-FLOAT_MIN;
      vx[secsp][j]= distribution_flux(1,secsp,vxloader[secsp][1]);
      vy[secsp][j] = v0y[secsp]+vty[secsp]*maxwellian(vyloader[secsp]);
      vz[secsp][j] = v0z[secsp]+vtz[secsp]*maxwellian(vzloader[secsp]);
      injection_push(secsp,j,frand());
      if  (x[secsp][j]>=xnc){
	j--;
	secountr--;
      }
    }

    np[secsp] += secountl +secountr;
    
    if(np[secsp] >= maxnp[secsp])
      printf("ADJUST(Secondaries): too many particles. MUST EXIT!");
  }
  //  jwall[isp] /= sp_k[isp]*dt;

  /* volume source of ionization */
  if (isp==ionsp){
    if (psource){
      if (t>tstrt)
	ionization(npold-np[ionsp]); /* Create ion/elec pair for each ion lost ?*/
    }
    if ((vol_source > 0)){
      this_time += vol_source;
      ionization((int)this_time);
      this_time -= (int)this_time;
    }
  }
		
} /* end ADJUST */

/***************************************************************/

void sterm(int count)
{
  static SCALAR ionization_vel[NSMAX], vgth;
  static int init_flag=1,ionsp;
  int i, ix, j, isp;
  SCALAR r;
  
  /************************************************/

  if (init_flag){
    vgth= sqrt(gtemp/Escale[ionsp]);	
    ionsp= ionspecies-1;          /* Fixing the indices into the array of species */
    for (isp=0; isp<nsp; isp++)
      ionization_vel[isp] = sqrt(2*gtemp/Escale[ionsp]);
    init_flag = 0;
  }	

  for (isp=0; isp<nsp; isp++)
    if ((np[isp] - count) >= maxnp[isp])
      {
	printf("np[isp] maxnp[isp] %d %d \n  ", np[isp], maxnp[isp]);
	puts("source term: too many particles, species ");
	putchar(isp+49);
	exit(1);
      }

  if (count < 0) /* add only when there is a loss */
    for (j=0; j<abs(count); j++)
      {
	r=frand();
	/* Loading the ions (ionsp)  randomly in the plasma */
	ix = np[ionsp];
	x[ionsp][ix]= (endpts[1]-endpts[0])*r+endpts[0];

	maxwellv(&vx[ionsp][ix], &vy[ionsp][ix], &vz[ionsp][ix], vgth);
				
	if(dde)
	  for(i=0; i<1; i++)
	    x[ionsp][ix] += dde*sin((2*M_PI*x[ionsp][ix]/xnc));

	np[ionsp]++;

	/* Loading the electron (isp=0) at the same position as the ion */
	isp=0;
	ix = np[isp];
	x[isp][ix]= x[ionsp][np[ionsp]-1];
	/* single energy ionization */

	maxwellv(&vx[ionsp][ix], &vy[ionsp][ix], &vz[ionsp][ix], vgth);

	np[isp]++;
      }   
}  /* end STERM  */

void ionization(int count)
{
  static SCALAR ionization_vel[NSMAX], vgth;
  static int init_flag=1,ionsp;
  int i, ix, j, isp;
  SCALAR theta, phi, r;
  SCALAR sintheta, cosphi, costheta, sinphi;
  
  /************************************************/

  if (init_flag){
    ionsp = ionspecies-1;  /* Fixing the indices into the array of species */ 
    if (ionsp<0){
      puts("photo ionization: must have ion species");
      exit(1);
    }
    vgth= sqrt(gtemp/Escale[ionsp]);	
    for (isp=0; isp<nsp; isp++)
      ionization_vel[isp] = m[0]*sqrt(ionization_energy/Escale[isp])/m[isp];
    init_flag = 0;
  }	

  for (isp=0; isp<nsp; isp++)
    if ((np[isp] + count) >= maxnp[isp])
      {
	printf("np[isp] maxnp[isp] %d %d \n  ", np[isp], maxnp[isp]);
	puts("photo ionization: too many particles, species ");
	putchar(isp+49);
	exit(1);
      }
	

  if (count > 0) 
    for (j=0; j<count; j++)
      {
	r=frand();
	/* Loading the ions (ionsp=ionspecies-1)  randomly in the plasma */
	ix = np[ionsp];
	x[ionsp][ix]= (endpts[1]-endpts[0])*r+endpts[0];

	/* single energy ionization with a background gas temperature*/
				
	maxwellv(&vx[ionsp][ix], &vy[ionsp][ix], &vz[ionsp][ix], vgth);
	phi =2*M_PI*frand();
	cosphi = cos(phi);
	sinphi = sin(phi);
	costheta = 1-2*frand();
	sintheta = sqrt(1-costheta*costheta);
	vx[ionsp][ix] += ionization_vel[ionsp]*sintheta*cosphi;
	vy[ionsp][ix] += ionization_vel[ionsp]*sintheta*sinphi;
	vz[ionsp][ix] += ionization_vel[ionsp]*costheta;
				
	if(dde)
	  for(i=0; i<1; i++)
	    x[ionsp][ix] += dde*sin((2*M_PI*x[ionsp][ix]/xnc));

	np[ionsp]++;

	/* Loading the electrons  (isp=0)  randomly in the plasma */
	isp=0;
	ix = np[isp];
	x[isp][ix]= x[ionsp][np[ionsp]-1];
	/* single energy ionization */

	//				theta += M_PI; //con mom
							   costheta *=-1;
	sintheta *=-1;
	vx[isp][ix] = ionization_vel[isp]*sintheta*cosphi;
	vy[isp][ix] = ionization_vel[isp]*sintheta*sinphi;
	vz[isp][ix] = ionization_vel[isp]*costheta;

	np[isp]++;
      }  

}  /* end ionization */

/***************************************************************/

void injection_push(int isp, int i, SCALAR del_t)
{
  int j;
  SCALAR vxinit, vyinit, vzinit, xinit, vxtemp, vytemp, vztemp, xtemp;
  SCALAR vxhalf, vyhalf, vzhalf;
  SCALAR k, a0;
  SCALAR sin2fW;

  SCALAR s;

  SCALAR deltaA, ep, f, fW, A1, A2, A3, A4;

  k = sp_k[isp];
  vxtemp = vxinit = vxhalf = vx[isp][i];
  vytemp = vyinit = vyhalf = vy[isp][i];
  vztemp = vzinit = vzhalf = vz[isp][i];
  xtemp = xinit = x[isp][i];

  j = xinit;
  s = xinit - j;
	
  /* a is for normalized for half a time step */
  /*
    if (s){
    a0 =a[j] + s*(a[j+1] - a[j]);
    deltaA = aold[j] + s*(aold[j+1] - aold[j])- a[j] + s*(a[j+1] - a[j]);
    ep= rho[j]+s*(rho[j+1] - rho[j]);
    }
    else {
    a0 = a[j];
    */
				    
  /*Assume that the particle is injected from the 
    right or left side of the simulation */

  a0=a[j];
  if (j<1){
    deltaA = a_scale[isp]*eold[0]-a0;
    ep = a_scale[isp]*rho[0]*dx/epsilon;
  }
  else{
    deltaA = a_scale[isp]*eold[1]-a0;
    ep = a_scale[isp]*rho[nc]*dx/epsilon;
  }


  if (b>0){
    /*position push*/
    a0*=2;
    ep*=2;
    deltaA*=2;
    f = del_t;
    fW=f*W[isp];
    A1=f*(3*a0+deltaA*(2*f-3)+ep*vxinit*f)*onesixth;
    A2=3*sin(fW)+fW*cos(fW);
    A3=3*cos(fW)-fW*sin(fW);
    x[isp][i] +=f*(vxinit+A1+twothirds*sin(fW)*sin_psi*(A2*(vzinit*cos_psi-(vxinit+0.5*A1)*sin_psi)+A3*vyinit));

    /*velocity push*/
    f = del_t-.5;
    fW=f*W[isp];
    A1=(W[isp]*sin4W[isp]+12*sqr(sin(2*fW)))*one12;
    A2=(W[isp]*sin22W[isp]-3*sin(4*fW))*onethird;
    A3=0.25*(4+ep*sqr(f));
    A4=0.25*f*(2*a0+deltaA*(f-1));
    sin2fW=sin(2*fW);
    vx[isp][i] += a0*f+one24*(deltaA*(12*f*(f-1)-1)-
			      2*sqr(sin_psi)*(a0*W[isp]*sin4W[isp]+24*A4*sqr(sin2fW))+
			      vxinit*(ep*(1+12*sqr(f))-4*sqr(sin_psi)*(12*A1+3*ep*sqr(f*sin2fW))))+
      sin_psi*(-vyinit*A2+
	       vzinit*2.0*cos_psi*A2);
    vy[isp][i] += sin_psi*(-onesixth*W[isp]*a0*cos22W[isp]-A4*sin(4*fW)+
			   vxinit*(onethird*W[isp]*sin22W[isp]-A3*sin(4*fW)))
      -2.0*vyinit*A1-vzinit*A2*cos_psi;
    vz[isp][i] += one12*cos_psi*(sin_psi*(a0*W[isp]*sin4W[isp]+24*sqr(sin2fW)*A4+
					  2*vxinit*(W[isp]*sin4W[isp]+12*sqr(sin2fW)*A3))+
				 12*vyinit*A2-24*vzinit*cos_psi*A1);
  }
  else{
    /*position push*/
    f = del_t;
    A1=f*(3*a0+deltaA*(2*f-3)+ep*vxinit*f)*onesixth;
    x[isp][i] +=f*(vxinit+A1);

    /*velocity push*/
    f = del_t-.5;
    vx[isp][i] += a0*f+one24*(deltaA*(12*f*(f-1)-1)+
			      vxinit*ep*(1+12*sqr(f)));
  } 
}

void injection_push_oldpdp(int isp, int i, SCALAR del_t)
{
  int j;
  SCALAR temp_t;
  SCALAR vxinit, vyinit, vzinit, xinit, vxtemp, vytemp, vztemp, xtemp;
  SCALAR vxhalf, vyhalf, vzhalf;
  SCALAR k, s, a0;
  SCALAR ax;

  k = sp_k[isp];
  vxtemp = vxinit = vxhalf = vx[isp][i];
  vytemp = vyinit = vyhalf = vy[isp][i];
  vztemp = vzinit = vzhalf = vz[isp][i];
  xtemp = xinit = x[isp][i];

  j = xinit;
  s = xinit - j;
			
  if (s)
    a0 =a[j] + s*(a[j+1] - a[j]);
  else a0 = a[j];
  if (b>0){
    temp_t = del_t-.5; /* half step already normalized in a array, tz, and tx.*/
    /****** mag acc ********/

    ax = 2*a0;

    vx[isp][i] += temp_t*ax;
		
    /*********** Update Position **************/
		
    x[isp][i]  += del_t*k*vxinit;
  }
}

void injection_push_FTSBP(int isp, int i, SCALAR del_t)
{
  int j;
  SCALAR temp_t;
  SCALAR vxinit, vyinit, vzinit, xinit, vxtemp, vytemp, vztemp, xtemp;
  SCALAR vxhalf, vyhalf, vzhalf;
  SCALAR k, s, atemp=0, a0;

  SCALAR tx,tz,sx,sz;
  SCALAR t;

  k = sp_k[isp];
  vxtemp = vxinit = vxhalf = vx[isp][i];
  vytemp = vyinit = vyhalf = vy[isp][i];
  vztemp = vzinit = vzhalf = vz[isp][i];
  xtemp = xinit = x[isp][i];

  j = xinit;
  s = xinit - j;
			
  if (s)
    a0 =a[j] + s*(a[j+1] - a[j]);
  else a0 = a[j];

  if (b>0){
    temp_t = del_t-.5;
    atemp = a[j] + s*(a[j+1] - a[j]);

    vx[isp][i] += temp_t*a0;

    t= tan(.5*b*qm[isp]*dt*temp_t*sp_k[isp]);
    s= 2*t/(1+t*t);
		
    tx= t*cosd(psi);
    tz= t*sind(psi);
    sx= s*cosd(psi);
    sz= s*sind(psi);

    /****** Advance velocity ******/
    /****** Boris rotation ********/
    vxhalf= vx[isp][i] +tz*vyinit;
    vyhalf= vyinit -tz*vx[isp][i] +tx*vzinit;
    vzhalf= vzinit -tx*vyinit;
		
    vx[isp][i] += sz*vyhalf;
    vy[isp][i] += -sz*vxhalf +sx*vzhalf;
    vz[isp][i] += -sx*vyhalf;

    vx[isp][i] += temp_t*a0;

    /****** Advance position ******/

    temp_t = del_t/2;
    atemp = a[j] + s*(a[j+1] - a[j]);

    vxtemp = vxinit;
    vxtemp += temp_t*a0;

    t= tan(.5*b*qm[isp]*dt*temp_t*sp_k[isp]);
    s= 2*t/(1+t*t);
		
    tx= t*cosd(psi);
    tz= t*sind(psi);
    sx= s*cosd(psi);
    sz= s*sind(psi);

    /****** Boris rotation ********/
    vyhalf= vyinit -tz*vxtemp +tx*vzinit;
    vxtemp += sz*vyhalf;

    vxtemp += temp_t*a0;
		
    x[isp][i] += vxtemp*del_t;

  }
}

void injection_push_MBP(int isp, int i, SCALAR del_t)
{
  int j;
  SCALAR vxinit, vyinit, vzinit, xinit, vxtemp, vytemp, vztemp, xtemp;
  SCALAR vxhalf, vyhalf, vzhalf;
  SCALAR k, a0;
  SCALAR sin2fW;

  SCALAR s;

  SCALAR dx;

  SCALAR f, fW, A1, A2, A3, A4;

  k = sp_k[isp];
  vxtemp = vxinit = vxhalf = vx[isp][i];
  vytemp = vyinit = vyhalf = vy[isp][i];
  vztemp = vzinit = vzhalf = vz[isp][i];
  xtemp = xinit = x[isp][i];

  j = xinit;
  s = xinit - j;
			
  /*Assume that the particle is injected from 
    the right or left side of the simulation */

  a0=a[j];

  if (b>0){
    /*position push*/
    a0*=2;
    f = del_t;
    fW= f*W[isp];
    A1=f*3*a0*onesixth;
    A2=3*sin(fW)+fW*cos(fW);
    A3=3*cos(fW)-fW*sin(fW);
    dx = f*(vxinit+A1+twothirds*sin(fW)*sin_psi*(A2*(vzinit*cos_psi-(vxinit+0.5*A1)*sin_psi)+A3*vyinit));
    x[isp][i] +=f*(vxinit+A1+twothirds*sin(fW)*sin_psi*(A2*(vzinit*cos_psi-(vxinit+0.5*A1)*sin_psi)+A3*vyinit));

    /*velocity push*/
    f = del_t-.5;
    fW=f*W[isp];
    A1=(W[isp]*sin4W[isp]+12*sqr(sin(2*fW)))*one12;
    A2=(W[isp]*sin22W[isp]-3*sin(4*fW))*onethird;
    A3=0.25*4;
    A4=0.25*f*(2*a0);
    sin2fW=sin(2*fW);
    vx[isp][i] += a0*f+one24*(-2*sqr(sin_psi)*(a0*W[isp]*sin4W[isp]+24*A4*sqr(sin2fW))+
			      vxinit*(-4*sqr(sin_psi)*12*A1))+
      sin_psi*(-vyinit*A2+
	       vzinit*2.0*cos_psi*A2);
    vy[isp][i] += sin_psi*(-onesixth*W[isp]*a0*cos22W[isp]-A4*sin(4*fW)+
			   vxinit*(onethird*W[isp]*sin22W[isp]-A3*sin(4*fW)))
      -2.0*vyinit*A1-vzinit*A2*cos_psi;
    vz[isp][i] += one12*cos_psi*(sin_psi*(a0*W[isp]*sin4W[isp]+24*sqr(sin2fW)*A4+
					  2*vxinit*(W[isp]*sin4W[isp]+12*sqr(sin2fW)*A3))+
				 12*vyinit*A2-24*vzinit*cos_psi*A1);
  }
 
}


void injection_push_FG(int isp, int i, SCALAR del_t)
{
  int j;
  SCALAR vxinit, vyinit, vzinit, xinit, vxtemp, vytemp, vztemp, xtemp;
  SCALAR vxhalf, vyhalf, vzhalf;
  SCALAR k, a0;
  SCALAR sin2fW;

  SCALAR s;

  SCALAR ep, f, fW, A1, A2, A3, A4;

  k = sp_k[isp];
  vxtemp = vxinit = vxhalf = vx[isp][i];
  vytemp = vyinit = vyhalf = vy[isp][i];
  vztemp = vzinit = vzhalf = vz[isp][i];
  xtemp = xinit = x[isp][i];

  j = xinit;
  s = xinit - j;
			
  /* Assume that the particle is injected from 
     the right or left side of the simulation */

  a0=a[j];
  if (j<1){
    ep = a_scale[isp]*rho[0]/epsilon*dx;
  }
  else{
    ep = a_scale[isp]*rho[nc]/epsilon;
  }


  if (b>0){
    /*position push*/
    /* a is for normalized for half a time step	a0*=2; */

    ep*=2;
    f = del_t;
    fW= f*W[isp];
    A1=f*(3*a0+ep*vxinit*f)*onesixth;
    A2=3*sin(fW)+fW*cos(fW);
    A3=3*cos(fW)-fW*sin(fW);
    x[isp][i] +=f*(vxinit+A1+twothirds*sin(fW)*sin_psi*(A2*(vzinit*cos_psi-(vxinit+0.5*A1)*sin_psi)+A3*vyinit));

    /*velocity push*/
    f = del_t-.5;
    fW=f*W[isp];
    A1=(W[isp]*sin4W[isp]+12*sqr(sin(2*fW)))*one12;
    A2=(W[isp]*sin22W[isp]-3*sin(4*fW))*onethird;
    A3=0.25*(4+ep*sqr(f));
    A4=0.5*f*a0;
    sin2fW=sin(2*fW);
    vx[isp][i] += a0*f+one24*(-2*sqr(sin_psi)*(a0*W[isp]*sin4W[isp]+24*A4*sqr(sin2fW))+
			      vxinit*(ep*(1+12*sqr(f))-4*sqr(sin_psi)*(12*A1+3*ep*sqr(f*sin2fW))))+
      sin_psi*(-vyinit*A2+
	       vzinit*2.0*cos_psi*A2);
    vy[isp][i] += sin_psi*(-onesixth*W[isp]*a0*cos22W[isp]-A4*sin(4*fW)+
			   vxinit*(onethird*W[isp]*sin22W[isp]-A3*sin(4*fW)))
      -2.0*vyinit*A1-vzinit*A2*cos_psi;
    vz[isp][i] += one12*cos_psi*(sin_psi*(a0*W[isp]*sin4W[isp]+24*sqr(sin2fW)*A4+
					  2*vxinit*(W[isp]*sin4W[isp]+12*sqr(sin2fW)*A3))+
				 12*vyinit*A2-24*vzinit*cos_psi*A1);
  }
 
}
