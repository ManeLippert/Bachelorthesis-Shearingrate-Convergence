#include "pdp1.h"
#include "xgrafix.h"

/***************************************************************/

void imp_move(int isp, const int EnergyFlag)
{
  register int i, j, k;
  register SCALAR  s;
  SCALAR atemp, ke_temp; 
  
  k = sp_k[isp];
  for (j=0; j<ng; j++) {
    a[j] = e[j]*a_scale[isp];
    if(EnergyFlag) sp_j_x[isp][j] = sp_j_y[isp][j] = sp_j_z[isp][j] = sp_ke_x[isp][j] =  sp_ke_y[isp][j] = sp_ke_z[isp][j] = 0.0;
  }
  
  for(i= 0; i<np[isp]; i++) { 
    /****** Post-push to time n *****/
    j = x[isp][i];
    s = x[isp][i] - j;
    atemp= a[j] + s*(a[j+1] - a[j]);
    
    x[isp][i] += 0.5*k*atemp;
    vx[isp][i]+= atemp;

    /* NOTE: Algorithmic bug: Particle may go out of bounds in prepush. */
    /****** Get diagnostics at time n ****/
    if (x[isp][i]< xnc) { /* exclude outsiders */
    j = x[isp][i];
    s = x[isp][i] - j;

    if(EnergyFlag) {
      ke_temp= vx[isp][i]*vx[isp][i] +vy[isp][i]*vy[isp][i] +vz[isp][i]*vz[isp][i];
      sp_j_x[isp][j]   += (1-s)*vx[isp][i];
      sp_j_x[isp][j+1] += s*vx[isp][i];
      sp_ke_x[isp][j]  += (1-s)*ke_temp;
      sp_ke_x[isp][j+1]+= s*ke_temp;
    }
    
    /****** Pre-push to time n+1  *******/
    atemp= a[j] + s*(a[j+1] - a[j]);
    x[isp][i] += k*(vx[isp][i] +0.5*atemp);
    vx[isp][i]+= atemp;
   }
  }
  
  /****** Fix diagnostics and end points and rescale ****/
  if(EnergyFlag) {
    sp_j_x[isp][0]  *= 2.;
    sp_j_x[isp][nc] *= 2.;
    sp_ke_x[isp][0] *= 2.;
    sp_ke_x[isp][nc]*= 2.;
    for (j=0; j<ng; j++) jdote[isp][j] += sp_j_x[isp][j]*e[j];
  }
}

/***************************************************************/

void exp_move(int isp, const int EnergyFlag)
{
  register int i, j;
  SCALAR k, s, vxtemp, vytemp, vztemp, atemp, ke_temp, vxold, vyold, vzold;
  SCALAR TEp;

  k = sp_k[isp];
  for (j=0; j<ng; j++) {
    a[j] = e[j]*a_scale[isp];
        if(EnergyFlag) sp_j_x[isp][j] = sp_j_y[isp][j] = sp_j_z[isp][j] = 
		     sp_ke_x[isp][j] = sp_ke_y[isp][j] = sp_ke_z[isp][j] =0.0;
		     
  }
    if(EnergyFlag) N_trapped[isp]=N_untrapped[isp]=E_untrapped[isp]=E_trapped[isp]=E_particles[isp]=0; 
  /*********** If external B field **********/
  if(b>0.0) {
    for (i=np[isp]-1; i>=0; i--) {
      j = x[isp][i];
      s = x[isp][i] - j;

      if(EnergyFlag) {
				vxtemp = vx[isp][i];
				vytemp = vy[isp][i];
				vztemp = vz[isp][i];
	}
			
      atemp = a[j] + s*(a[j+1] - a[j]);
      vx[isp][i] += atemp;

      /****** Boris rotation ********/
      vxold= vx[isp][i] +tz[isp]*vy[isp][i];
      vyold= vy[isp][i] -tz[isp]*vx[isp][i] +tx[isp]*vz[isp][i];
      vzold= vz[isp][i] -tx[isp]*vy[isp][i];
      
      vx[isp][i] += sz[isp]*vyold;
      vy[isp][i] += -sz[isp]*vxold +sx[isp]*vzold;
      vz[isp][i] += -sx[isp]*vyold;
      
      /****** Advance velocity and position ******/
      vx[isp][i] += atemp;
      x[isp][i]  += k*vx[isp][i];

      /****** Get diagnostics at time n ****/
      if(EnergyFlag) {
				vxtemp += vx[isp][i];
	vytemp += vy[isp][i];
	vztemp += vz[isp][i];
	sp_j_x[isp][j]   += (1-s)*vxtemp;
	sp_j_x[isp][j+1] += s*vxtemp;
	sp_j_y[isp][j]   += (1-s)*vytemp;
	sp_j_y[isp][j+1] += s*vytemp;
	sp_j_z[isp][j]   += (1-s)*vztemp;
	sp_j_z[isp][j+1] += s*vztemp;
	vxtemp *= vxtemp;
	vytemp *= vytemp;
	vztemp *= vztemp;
	sp_ke_x[isp][j] += (1-s)*vxtemp;
	sp_ke_x[isp][j+1] += s*vxtemp;
	sp_ke_y[isp][j] += (1-s)*vytemp;
	sp_ke_y[isp][j+1] += s*vytemp;
	sp_ke_z[isp][j] += (1-s)*vztemp;
	sp_ke_z[isp][j+1] += s*vztemp;
	
	/** trapping diagnostics **/
	/*
	  TEp=Escale[isp]*ke_temp+phi[0]-((1-s)*phi[j]+s*phi[j+1]);
	  E_particles[isp]+=TEp;
	  if (TEp>0){
	  N_untrapped[isp]++;
	  E_untrapped[isp]+=TEp;
	  }
	  else{
	  N_trapped[isp]++;
	  E_trapped[isp]+=TEp;
	  }
	  */
		}
	}
}
  /*********** If NO external B field **********/
  else {
    for (i=0;i<np[isp]; i++)
      {
	j = x[isp][i];
	s = x[isp][i] - j;
	
	if(EnergyFlag)
	  vxtemp = vx[isp][i];
	  
	/****** Advance velocity and position ******/
	vx[isp][i]+= a[j] + s*(a[j+1] - a[j]);
	x[isp][i] += k*vx[isp][i]; 

	/****** Get diagnostics at time n ****/     
		if(EnergyFlag) {
	  vytemp = vy[isp][i];
	  vztemp = vz[isp][i];
	  vxtemp += vx[isp][i];
	  vytemp += vytemp;
	  vztemp += vztemp;
	  sp_j_x[isp][j]   += (1-s)*vxtemp;
	  sp_j_x[isp][j+1] += s*vxtemp;
	  sp_j_y[isp][j]   += (1-s)*vytemp;
	  sp_j_y[isp][j+1] += s*vytemp;
	  sp_j_z[isp][j]   += (1-s)*vztemp;
	  sp_j_z[isp][j+1] += s*vztemp;
	  vxtemp *= vxtemp;
	  vytemp *= vytemp;
	  vztemp *= vztemp;
	  sp_ke_x[isp][j] += (1-s)*vxtemp;
	  sp_ke_x[isp][j+1] += s*vxtemp;
	  sp_ke_y[isp][j] += (1-s)*vytemp;
	  sp_ke_y[isp][j+1] += s*vytemp;
	  sp_ke_z[isp][j] += (1-s)*vztemp;
	  sp_ke_z[isp][j+1] += s*vztemp;
	}
	
      }
  }
  
  /****** Fix diagnostics and end points and rescale ****/
  if(EnergyFlag) {
    sp_j_x[isp][0]  *= 2.;
    sp_j_x[isp][nc] *= 2.;
    sp_j_y[isp][0]  *= 2.;
    sp_j_y[isp][nc] *= 2.;
    sp_j_z[isp][0]  *= 2.;
    sp_j_z[isp][nc] *= 2.;
    sp_ke_x[isp][0] *= 2.;
    sp_ke_x[isp][nc]*= 2.;
    sp_ke_y[isp][0] *= 2.;
    sp_ke_y[isp][nc]*= 2.;
    sp_ke_z[isp][0] *= 2.;
    sp_ke_z[isp][nc]*= 2.;
    for (j=0; j<ng; j++) jdote[isp][j] += sp_j_x[isp][j]*e[j];
  }
}

/***************************************************************/

