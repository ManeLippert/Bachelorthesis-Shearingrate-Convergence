#include "pdp1.h"
#include "xgrafix.h"

#define   E_ELEC1_MAX     100
#define   E_ELEC2_MAX     100
#define   E_ELEC3_MAX     100
#define   E_ONEG1_MAX      10
#define   E_ONEG2_MAX      10
#define   E_OTWO_MAX      500

void onewvel(SCALAR energy, SCALAR vel, SCALAR *vx, SCALAR *vy, SCALAR *vz, int e_flag);
SCALAR osigma1(SCALAR), osigma2(SCALAR), osigma3(SCALAR), osigma4(SCALAR), osigma5(SCALAR);
SCALAR osigma6(SCALAR), osigma7(SCALAR), osigma8(SCALAR), osigma9(SCALAR), osigma10(SCALAR);
SCALAR osigma11(SCALAR), osigma12(SCALAR), osigma13(SCALAR), osigma14(SCALAR), osigma15(SCALAR);
SCALAR osigma16(SCALAR), osigma17(SCALAR), osigma18(SCALAR), osigma19(SCALAR), osigma20(SCALAR);
SCALAR osigma21(SCALAR), osigma22(SCALAR);
SCALAR gden, vgth;
int ecolsp, onegsp, otwosp;

/****************************************************************/

void oxygenmcc(int isp)
{
  static int init_flag=1;
  static SCALAR ecol_extra,  max_sigmav_e, col_prob_e;
  static SCALAR ecol_extra2, max_sigmav_e2;
  static SCALAR ecol_extra3, max_sigmav_e3;
  static SCALAR oneg_extra,  max_sigmav_oneg, col_prob_oneg;
  static SCALAR oneg_extra2, max_sigmav_oneg2;
  static SCALAR otwo_extra,  max_sigmav_otwo, col_prob_otwo;
  static SCALAR sigma_scale;
  
  register int i, j, k, l, index;
  int s, N, nnp, foundflag;
  SCALAR sum_sigma, chi, sinchi, coschi, phi1, cosphi, sinphi;
  SCALAR temp, engy, eengy, rengy, sigma_total;
  SCALAR vneutx, vneuty, vneutz, vtempx, vtempy, vtempz;
  SCALAR del, dum, random, vel, nu_otwo_max, nu_oneg_max;
  
  SCALAR up1, up2, up3, mag;
  SCALAR r11, r12, r13, r21, r22, r23, r31, r32, r33;
  
  if (init_flag) {
    gden = NperTORR*pressure/(gtemp+DBL_MIN);  /* Calculating the gas density */
    ecolsp = 0;
    otwosp = 1;
    onegsp = 2;
    
    vgth= sqrt(gtemp/Escale[otwosp]);
    sigma_scale = area*dt/nc2p;
    
    /**********************************************/
    /* Calculating the null collision probability */
    
    /*************************/
    max_sigmav_e = 0.0;
    for (engy=0; engy<E_ELEC1_MAX; engy +=0.1)
      max_sigmav_e= max(max_sigmav_e, sqrt(2*1.602e-19*engy/m[ecolsp])
			*(osigma1(engy)+osigma2(engy)+osigma3(engy)+osigma4(engy)
			  +osigma5(engy)+osigma6(engy)+osigma7(engy)+osigma8(engy)
			  +osigma9(engy)+osigma10(engy)+osigma11(engy)+osigma12(engy)
			  +osigma13(engy)+osigma14(engy)+osigma15(engy)));
    
    col_prob_e = 1 -exp(-max_sigmav_e*gden*dt*sp_k[ecolsp]);
    
    /*************************/
    max_sigmav_e2 = 0.0;
    for (engy=0; engy<E_ELEC2_MAX; engy +=0.1)
      max_sigmav_e2= max(max_sigmav_e2, sqrt(2*1.602e-19*engy/m[ecolsp])*osigma21(engy));
    
    /*************************/
    max_sigmav_e3 = 0.0;
    for (engy=0; engy<E_ELEC3_MAX; engy +=0.1)
      max_sigmav_e3= max(max_sigmav_e3, sqrt(2*1.602e-19*engy/m[ecolsp])*osigma22(engy));
    
    /*************************/
    max_sigmav_oneg = 0.0;
    for (engy=0; engy<E_ONEG1_MAX;  engy +=0.1)
      max_sigmav_oneg= max(max_sigmav_oneg, sqrt(2*1.602e-19*engy/m[onegsp])
			   *(osigma16(engy)+osigma17(engy)));
    
    col_prob_oneg = 1 -exp(-max_sigmav_oneg*gden*dt*sp_k[onegsp]);
    
    /*************************/
    max_sigmav_oneg2 = 0.0;
    for (engy=0; engy<E_ONEG2_MAX; engy +=0.1)
      max_sigmav_oneg2= max(max_sigmav_oneg2, sqrt(2*1.602e-19*engy/m[onegsp])*osigma20(engy));
    
    /*************************/
    max_sigmav_otwo = 0.0;
    for (engy=0; engy<E_OTWO_MAX; engy +=0.1)
      max_sigmav_otwo= max(max_sigmav_otwo, sqrt(2*1.602e-19*engy/m[otwosp])*osigma18(engy));
    
    col_prob_otwo = 1 -exp(-max_sigmav_otwo*gden*dt*sp_k[otwosp]);
    
    init_flag = 0;
  }

  /********************************/
  /* electron collisions with O2  */

  if(isp==ecolsp) {
    /**** First Clear the diagnostics  ******/
    if(theRunWithXFlag) {
      for(i=0; i<15; i++)
	for (j=0; j<ng; j++)
	  mccrate[i][j] = 0.0;
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
			if (dum){
				engy= Escale[ecolsp]*dum;
				vel = sqrt(dum);
				sigma_total = max_sigmav_e/(vel*vscale);
				random= frand();
      
      /*********************************************************************/
      /* determine the type of collision and calculate new velocities, etc */
      
      /**************************************************/
      /* momentum transfer                              */
      
      if (random <= (sum_sigma =osigma6(engy))/sigma_total) {
				/* first normalize vel */
				vx[ecolsp][j] /= vel;
				vy[ecolsp][j] /= vel;
				vz[ecolsp][j] /= vel;
				
				/* then scatter the electron */
				onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 1);
				
				/* collect diagnostics */
				if(theRunWithXFlag) {
					s= x[ecolsp][j];
					del= x[ecolsp][j] - s;
					mccrate[0][s]  += (!s) ? 2*(1-del) : 1-del;
					mccrate[0][s+1]+= (s== nc-1) ? 2*del : del;
				}
      }
      /**************************************************/
      /* Rotational excitation       (E loss = 0.02 eV) */
      
      else if (engy >= .02 && random <= (sum_sigma+=osigma5(engy))/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= 0.02;
        vel   = sqrt(engy/Escale[ecolsp]);

	/* then scatter the electron */
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[1][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[1][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /* viberational excitation v=1 (E loss = 0.19 eV) */
      
      else if (engy >= .19 && random <= (sum_sigma+=osigma1(engy))/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= 0.19;
        vel   = sqrt(engy/Escale[ecolsp]);
	
	/* then scatter the electron */
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[2][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[2][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /* viberational excitation v=2 (E loss = 0.38 eV) */
      
      else if (engy >= .38 && random <= (sum_sigma+=osigma2(engy))/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= 0.38;
        vel   = sqrt(engy/Escale[ecolsp]);
	
	/* then scatter the electron */
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[3][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[3][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /* viberational excitation v=3 (E loss = 0.57 eV) */
      
      else if (engy >= .57 && random <= (sum_sigma+=osigma3(engy))/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= 0.57;
        vel   = sqrt(engy/Escale[ecolsp]);

	/* then scatter the electron */
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[4][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[4][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /* viberational excitation v=4 (E loss = 0.75 eV) */
      
      else if (engy >= .75 && random <= (sum_sigma+=osigma4(engy))/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= 0.75;
        vel   = sqrt(engy/Escale[ecolsp]);

	/* then scatter the electron */
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[5][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[5][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /*  O2 SING DELTA             (E LOSS = 0.977 eV) */
      
      else if (engy >= .977 && random <= (sum_sigma+=osigma7(engy))/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= 0.977;
        vel   = sqrt(engy/Escale[ecolsp]);

	/* then scatter the electron */
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[6][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[6][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /*  O2 B SINGLET SIGMA        (E LOSS = 1.627 eV) */
      
      else if (engy >= 1.627 && random <= (sum_sigma+=osigma8(engy))/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= 1.627;
        vel   = sqrt(engy/Escale[ecolsp]);

	/* then scatter the electron */
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[7][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[7][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /*******************************************************/
      /*  e + O2 -> O + O-  Dissociative attachment (Et=4.2) */
      
      else if (engy >= 4.2 && random <= (sum_sigma+=osigma15(engy))/sigma_total) {
	/* create a neg. ion and assign vel */
	k= np[onegsp];
	vx[onegsp][k]= 1-2*frand();
	vy[onegsp][k]= 1-2*frand();
	vz[onegsp][k]= 1-2*frand();
	vel = sqrt((engy -4.2)/(2.0*Escale[onegsp]
				*(vx[onegsp][k]*vx[onegsp][k]
				  +vy[onegsp][k]*vy[onegsp][k]
				  +vz[onegsp][k]*vz[onegsp][k])));
	vx[onegsp][k]*= vel;
	vy[onegsp][k]*= vel;
	vz[onegsp][k]*= vel;
	x[onegsp][k]  = x[ecolsp][j];
	
	s= x[onegsp][k];
	del= x[onegsp][k] - s;
	sp_n_mcc[onegsp][s]  += 1 -del;
	sp_n_mcc[onegsp][s+1]+= del;
	
	/* remove the electron */
	x[ecolsp][j] =  x[ecolsp][nnp+N-1];
	vx[ecolsp][j]= vx[ecolsp][nnp+N-1];
	vy[ecolsp][j]= vy[ecolsp][nnp+N-1];
	vz[ecolsp][j]= vz[ecolsp][nnp+N-1];
	
	/* adjust the indecies and the total numbers */
        /* move particles if ionization occurred */
	if((nnp+N) != np[ecolsp]) {
	  x[ecolsp][nnp+N-1] =  x[ecolsp][np[ecolsp]-1];
	  vx[ecolsp][nnp+N-1]= vx[ecolsp][np[ecolsp]-1];
          vy[ecolsp][nnp+N-1]= vy[ecolsp][np[ecolsp]-1];
          vz[ecolsp][nnp+N-1]= vz[ecolsp][np[ecolsp]-1];
	}
	
	/* adjust the indecies and the total numbers */
	if(++np[onegsp] > maxnp[onegsp]) {
	  puts("OXYGENMCC(Dis. Attachment): too many particles. MUST EXIT!");
	  exit(-1);
	}
	N--;  j--;  np[ecolsp]--;
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  mccrate[8][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[8][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /* 4.5 E loss                   (E LOSS = 4.5 eV) */
      
      else if (engy >= 4.5 && random <= (sum_sigma+=osigma9(engy))/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= 4.5;
        vel   = sqrt(engy/Escale[ecolsp]);

	/* then scatter the electron */
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[9][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[9][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /* 6.0 E loss                   (E LOSS = 6.0 eV) */
      
      else if (engy >= 6.0 && random <= (sum_sigma+=osigma10(engy))/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= 6.0;
        vel   = sqrt(engy/Escale[ecolsp]);

	/* then scatter the electron */
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[10][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[10][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /* 8.4 E loss                   (E LOSS = 8.4 eV) */
      
      else if (engy >= 8.4 && random <= (sum_sigma+=osigma11(engy))/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= 8.4;
        vel   = sqrt(engy/Escale[ecolsp]);
	
	/* then scatter the electron */
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[11][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[11][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /* 9.97 E loss                 (E LOSS = 9.97 eV) */
      
      else if (engy >= 9.97 && random <= (sum_sigma+=osigma12(engy))/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= 9.97;
        vel   = sqrt(engy/Escale[ecolsp]);

	/* then scatter the electron */
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[12][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[12][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /* e + O2 -> e + e + O2+ Ion. (E LOSS = 12.06 eV) */
      
      else if (engy >= 12.06 && random <= (sum_sigma+=osigma13(engy))/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	/* subtract the ion. energy and partition the remaining energy */
	engy -= 12.06;
	rengy = 10.0*tan(frand()*atan(engy/20.0));
	engy -= rengy;
	
	/* scatter the created electron */
	vel = sqrt(fabs(rengy)/Escale[ecolsp]);
	k = np[ecolsp];
	vx[ecolsp][k] = vx[ecolsp][j];
	vy[ecolsp][k] = vy[ecolsp][j];
	vz[ecolsp][k] = vz[ecolsp][j];
	onewvel(rengy, vel, &vx[ecolsp][k], &vy[ecolsp][k], &vz[ecolsp][k], 0);
	x[ecolsp][k]  = x[ecolsp][j];
	
	/* assign velocities to the created ion */
	l = np[otwosp];
	maxwellv(&vx[otwosp][l], &vy[otwosp][l], &vz[otwosp][l], vgth);
	x[otwosp][l] = x[ecolsp][j];
	
	s= x[ecolsp][j];
	del= x[ecolsp][j] - s;
	sp_n_mcc[otwosp][s]  += 1 -del;
	sp_n_mcc[otwosp][s+1]+= del;
	
	/* finally scatter the incident electron */
	vel = sqrt(fabs(engy)/Escale[ecolsp]);
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	if(++np[otwosp] > maxnp[otwosp] || ++np[ecolsp] > maxnp[ecolsp]) {
	  puts("OXYGENMCC(Ionization): too many particles. MUST EXIT!");
	  exit(-1);
	}
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  mccrate[13][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[13][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /* 130 nm line excitation      (E LOSS = 14.7 eV) */
      
      else if (engy >= 14.7 && random <= (sum_sigma+=osigma14(engy))/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= 14.7;
        vel   = sqrt(engy/Escale[ecolsp]);
	
	/* then scatter the electron */
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[14][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[14][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
     }
    }
  }
  /******************************************************/
  /* e + O2+ -> ....                                    */
  
  if(isp==ecolsp) {
    /**** First Clear the diagnostics  ******/
    if(theRunWithXFlag) {
      for (j=0; j<ng; j++)
	mccrate[15][j] = 0.0;
    }
    /**** Now do the collisions *****/
    for(nu_otwo_max= 0.0, i=0; i<ng; i++) nu_otwo_max = max(nu_otwo_max, sp_n[otwosp][i]);
    nu_otwo_max *= max_sigmav_e2*nc2p/area/dx;
    
    ecol_extra2 += np[ecolsp]*(1 -exp(-nu_otwo_max*dt*sp_k[ecolsp]));
    N = ecol_extra2;
    ecol_extra2 -= N;
    
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
      /** choose a positive ion (O2+) **/
      vtempx= vx[otwosp][0];
      vtempy= vy[otwosp][0];
      vtempz= vz[otwosp][0];
      for(i=0, foundflag=1; i <np[otwosp] && foundflag; i++) {
	     if(fabs(x[ecolsp][j]-x[otwosp][i]) <= 1.0) {
	       vtempx= vx[otwosp][i];
	       vtempy= vy[otwosp][i];
	       vtempz= vz[otwosp][i];
	       foundflag= 0; 
	     }
      }
      i--;
      /* if(i== np[otwosp]-1) i=0; */
      /* if O2+ ion is not found near e-, there is no collision.
         So, continue  onto the next colliding e- */
      if (foundflag) continue; 
     
      /** transform into reference frame where v_otwo= 0 **/
      vtempx -= vx[ecolsp][j];
      vtempy -= vy[ecolsp][j];
      vtempz -= vz[ecolsp][j];
      dum = (vtempx*vtempx +vtempy*vtempy +vtempz*vtempz);
      engy= Escale[ecolsp]*dum;
      vel = sqrt(dum);
      k = x[otwosp][i];
      sigma_total = sigma_scale*nu_otwo_max/(vel*sp_n[otwosp][k]);
      random= frand();
      
      /**************************************************/
      /* e + O2+ -> O + O  (Dissociative Recombination) */
      
      if (random <= osigma21(engy)/sigma_total) {
	s = x[otwosp][j];
	del= x[otwosp][j] - s;
	
	/****  remove electron ******/
	x[ecolsp][j] =  x[ecolsp][nnp+N-1];
	vx[ecolsp][j]= vx[ecolsp][nnp+N-1];
	vy[ecolsp][j]= vy[ecolsp][nnp+N-1];
	vz[ecolsp][j]= vz[ecolsp][nnp+N-1];
	
	/****  remove pos. ion ******/
	l = np[otwosp]-1;
	x[otwosp][i] =  x[otwosp][l];
	vx[otwosp][i]= vx[otwosp][l];
	vy[otwosp][i]= vy[otwosp][l];
	vz[otwosp][i]= vz[otwosp][l];
	np[otwosp]--;
	
	sp_n_mcc[otwosp][s]  -= 1 -del;
	sp_n_mcc[otwosp][s+1]-= del;
	
	/** adjust the indecies and total numbers **/
	N--;  j--;  np[ecolsp]--;

	/* collect diagnostics */
	if(theRunWithXFlag) {
	  mccrate[15][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[15][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
    }
  }
  /******************************************************/
  /* e + O- -> ....                                    */
  
  if(isp==ecolsp) {
    /**** First Clear the diagnostics  ******/
    if(theRunWithXFlag) {
      for (j=0; j<ng; j++)
	mccrate[16][j] = 0.0;
    }
    /**** Now do the collisions *****/
    for(nu_oneg_max= 0.0, i=0; i<ng; i++) nu_oneg_max = max(nu_oneg_max, sp_n[onegsp][i]);
    nu_oneg_max *= max_sigmav_e3*nc2p/area/dx;
    
    ecol_extra3 += np[ecolsp]*(1 -exp(-nu_oneg_max*dt*sp_k[ecolsp]));
    N = ecol_extra3;
    ecol_extra3 -= N;
    
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
      /** choose a positive ion (O2+) **/
      vtempx= vx[onegsp][0];
      vtempy= vy[onegsp][0];
      vtempz= vz[onegsp][0];
      for(i=0, foundflag=1; i <np[onegsp] && foundflag; i++) {
	     if(fabs(x[ecolsp][j]-x[onegsp][i]) <= 1.0) {
	       vtempx= vx[onegsp][i];
	       vtempy= vy[onegsp][i];
	       vtempz= vz[onegsp][i];
	       foundflag= 0; 
	     }
      }
      i--;
      /* if(i== np[onegsp]-1) i=0; */
      /* if O- ion is not found near e-, there is no collision.
         So, continue  onto the next colliding e- */
      if (foundflag) continue; 

      /** transform into reference frame where v_oneg= 0 **/
      vx[ecolsp][j] -= vtempx;
      vy[ecolsp][j] -= vtempy;
      vz[ecolsp][j] -= vtempz;
      dum = (vx[ecolsp][j]*vx[ecolsp][j] +vy[ecolsp][j]*vy[ecolsp][j]
	     +vz[ecolsp][j]*vz[ecolsp][j]);
      engy= Escale[ecolsp]*dum;
      vel = sqrt(dum);
      k = x[onegsp][i];
      sigma_total = sigma_scale*nu_oneg_max/(vel*sp_n[onegsp][k]);
      random= frand();
      
      /**************************************************/
      /* e + O- -> O + 2 e (Electron Impact Detachment) */
      
      if (engy >= 1.465 && random <= osigma22(engy)/sigma_total) {
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	/* subtract the threshold energy and partition the remaining energy */
	engy -= 1.465;
	rengy = 10.0*tan(frand()*atan(engy/20.0));
	engy -= rengy;
	
	/*** scatter the ejected electron ***/
	vel = sqrt(fabs(rengy)/Escale[ecolsp]);
	k = np[ecolsp];
	vx[ecolsp][k] = vx[ecolsp][j];
	vy[ecolsp][k] = vy[ecolsp][j];
	vz[ecolsp][k] = vz[ecolsp][j];
	onewvel(rengy, vel, &vx[ecolsp][k], &vy[ecolsp][k], &vz[ecolsp][k], 0);
	/** transform back into lab frame ***/
	vx[ecolsp][k]+= vtempx;
	vy[ecolsp][k]+= vtempy;
	vz[ecolsp][k]+= vtempz;
	x[ecolsp][k]  = x[ecolsp][j];
	
	/*** scatter the incident electron ***/
	vel = sqrt(fabs(engy)/Escale[ecolsp]);
	onewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
	
	/****  remove neg. ion ******/
	s  = x[onegsp][i];
	del= x[onegsp][i] - s;
	sp_n_mcc[onegsp][s]  -= 1 -del;
	sp_n_mcc[onegsp][s+1]-= del;
	
	l = np[onegsp]-1;
	x[onegsp][i] =  x[onegsp][l];
	vx[onegsp][i]= vx[onegsp][l];
	vy[onegsp][i]= vy[onegsp][l];
	vz[onegsp][i]= vz[onegsp][l];
	np[onegsp]--;
	
	/** adjust the indecies and total numbers **/
	if(++np[ecolsp] > maxnp[ecolsp]) {
	  puts("OXYGENMCC(E. Impact Detach.): too many particles. MUST EXIT!");
	  exit(-1);
	}
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  mccrate[16][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[16][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /** transform back into lab frame                **/
      vx[ecolsp][j] += vtempx;
      vy[ecolsp][j] += vtempy;
      vz[ecolsp][j] += vtempz;
    }
  }
  /******************************************************/
  /* 0- + O2 -> ..........                              */
  
  if(isp==onegsp) {
    /**** First Clear the diagnostics  ******/
    if(theRunWithXFlag) {
      for (j=0; j<ng; j++)
	mccrate[17][j] = mccrate[18][j] = 0.0;
    }
    /**** Now do the collisions *****/
    oneg_extra += np[onegsp]*col_prob_oneg;
    N = oneg_extra;
    oneg_extra -= N;
    
    nnp = np[onegsp];
    for(j=0; j< N; j++) {
      index= nnp*frand();
      nnp--;
      temp = x[onegsp][nnp];
      x[onegsp][nnp] = x[onegsp][index];
      x[onegsp][index] = temp;
      
      temp = vx[onegsp][nnp];
      vx[onegsp][nnp] = vx[onegsp][index];
      vx[onegsp][index] = temp;
      
      temp = vy[onegsp][nnp];
      vy[onegsp][nnp] = vy[onegsp][index];
      vy[onegsp][index] = temp;
      
      temp = vz[onegsp][nnp];
      vz[onegsp][nnp] = vz[onegsp][index];
      vz[onegsp][index] = temp;
    }
    
    for(j=nnp; j<nnp+N; j++) {
      /** choose a neutral from bkgd **/
      maxwellv(&vneutx, &vneuty, &vneutz, vgth);

      /** transform into reference frame where v_neutral= 0 **/
      vx[onegsp][j] -= vneutx;
      vy[onegsp][j] -= vneuty;
      vz[onegsp][j] -= vneutz;
      dum = (vx[onegsp][j]*vx[onegsp][j] +vy[onegsp][j]*vy[onegsp][j]
	     +vz[onegsp][j]*vz[onegsp][j]);
      engy= Escale[onegsp]*dum;
      vel = sqrt(dum);
      sigma_total = max_sigmav_oneg/(vel*vscale);
      random= frand();
      
      /**************************************************/
      /*  O- + O2 -> O- + O2      (Momentum Transfer)   */
      
      if (random <= (sum_sigma =osigma17(engy))/sigma_total) {
	/** determine the scattering angles **/
	random= frand();
	chi = atan(sqrt(random*(1-random))/(.75 -random));
	if(chi < 0.0) chi += M_PI;
	coschi= cos(chi);
	sinchi= sin(chi);
	
	phi1  = 2*M_PI*frand();
	cosphi= cos(phi1);
	sinphi= sin(phi1);
	
	/** determine the rotation matrix **/
	r13 = vx[onegsp][j]/vel;
	r23 = vy[onegsp][j]/vel;
	r33 = vz[onegsp][j]/vel;
	
	if(r33 == 1) { up1= 0;  up2= 1;  up3= 0; }
	else         { up1= 0;  up2= 0;  up3= 1; }
	
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
	
	vx[onegsp][j]= r11*sinchi*cosphi +r12*sinchi*sinphi +r13*coschi;
	vy[onegsp][j]= r21*sinchi*cosphi +r22*sinchi*sinphi +r23*coschi; 
	vz[onegsp][j]= r31*sinchi*cosphi +r32*sinchi*sinphi +r33*coschi;      
	
	/** adjusting velocity magnitude **/
	vel *= (coschi +sqrt(coschi*coschi+3))/3.0;
	vx[onegsp][j] *= vel;
	vy[onegsp][j] *= vel;
	vz[onegsp][j] *= vel;
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[onegsp][j];
	  del= x[onegsp][j] - s;
	  mccrate[17][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[17][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /*  O- + O2 -> O2 + O + e    (Detachment)         */ 
      
      else if (engy >= 1.465 && random <= (sum_sigma+=osigma16(engy))/sigma_total) {
	/** determine the energies  **/
	engy -= 1.465;
	if(engy > 2.0) { 
	  eengy= 1.465;        /* oengy= 2*(engy -1.465)*frand()/3.0; */
	}
	else {
	  eengy= engy*frand(); /* oengy= 2*(engy -eengy)*frand()/3.0; */
	}
	
	/**** assign vel to created electron  ***/
	k = np[ecolsp];
	vx[ecolsp][k]= 1-2*frand();
	vy[ecolsp][k]= 1-2*frand();
	vz[ecolsp][k]= 1-2*frand();
	vel = sqrt(eengy/(Escale[ecolsp]
			  *(vx[ecolsp][k]*vx[ecolsp][k] +vy[ecolsp][k]*vy[ecolsp][k]
			    +vz[ecolsp][k]*vz[ecolsp][k])));
	
	vx[ecolsp][k] = vel*vx[ecolsp][k] +vneutx;
	vy[ecolsp][k] = vel*vy[ecolsp][k] +vneuty;
	vz[ecolsp][k] = vel*vz[ecolsp][k] +vneutz;
	x[ecolsp][k]  = x[onegsp][j];
	
	s  = x[ecolsp][k];
	del= x[ecolsp][k] - s;
	sp_n_mcc[ecolsp][s]  += 1 -del;
	sp_n_mcc[ecolsp][s+1]+= del;
	
	/****  remove neg. ion ******/
	x[onegsp][j] =  x[onegsp][nnp+N-1];
	/* Add negative ion in the neutral frame of reference */
	vx[onegsp][j]= vx[onegsp][nnp+N-1] -vneutx;
	vy[onegsp][j]= vy[onegsp][nnp+N-1] -vneuty;
	vz[onegsp][j]= vz[onegsp][nnp+N-1] -vneutz;
	
	/** adjust the indecies and total numbers **/
	if(++np[ecolsp] > maxnp[ecolsp]) {
	  puts("OXYGENMCC(Detachment): too many particles. MUST EXIT!");
	  exit(-1);
	}
	N--;  j--;  np[onegsp]--;
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  mccrate[18][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[18][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /**************************************************/
      /** transform back into lab frame                **/
      vx[onegsp][j] += vneutx;
      vy[onegsp][j] += vneuty;
      vz[onegsp][j] += vneutz;
    }
  }
  /******************************************************/
  /* O- + O2+ -> ....                                   */
  
  if(isp==onegsp) {
    /**** First Clear the diagnostics  ******/
    if(theRunWithXFlag) {
      for (j=0; j<ng; j++)
	mccrate[19][j] = 0.0;
    }
    /**** Now do the collisions *****/
    for(nu_otwo_max= 0.0, i=0; i<ng; i++) nu_otwo_max = max(nu_otwo_max, sp_n[otwosp][i]);
    nu_otwo_max *= max_sigmav_oneg2*nc2p/area/dx; 
    
    oneg_extra2 += np[onegsp]*(1 -exp(-nu_otwo_max*dt*sp_k[onegsp]));
    N = oneg_extra2; 
    oneg_extra2 -= N;
    
    nnp = np[onegsp];
    for(j=0; j< N; j++)
    {
      index= nnp*frand();
      nnp--;
      temp = x[onegsp][nnp];
      x[onegsp][nnp] = x[onegsp][index];
      x[onegsp][index] = temp;
      
      temp = vx[onegsp][nnp];
      vx[onegsp][nnp] = vx[onegsp][index];
      vx[onegsp][index] = temp;
      
      temp = vy[onegsp][nnp];
      vy[onegsp][nnp] = vy[onegsp][index];
      vy[onegsp][index] = temp;
      
      temp = vz[onegsp][nnp];
      vz[onegsp][nnp] = vz[onegsp][index];
      vz[onegsp][index] = temp;
    }
    
    for(j=nnp; j<nnp+N; j++)
    {
      /** choose a positive ion (O2+) **/
      vtempx= vx[otwosp][0];
      vtempy= vy[otwosp][0];
      vtempz= vz[otwosp][0];
      for(i=0, foundflag=1; i <np[otwosp] && foundflag; i++)
      {
	if(fabs(x[onegsp][j]-x[otwosp][i]) <= 1.0)
	{
	  vtempx= vx[otwosp][i];
	  vtempy= vy[otwosp][i];
	  vtempz= vz[otwosp][i];
	  foundflag= 0;
	}
      }
      i--; 
      /* if(i== np[otwosp]-1) i=0; */
      /* if O2+ ion is not found near O-, there is no collision.
         So, continue  onto the next colliding O- */
      if (foundflag) continue; 
      
      /** transform into reference frame where v_otwo= 0 **/
      vtempx -= vx[onegsp][j];
      vtempy -= vy[onegsp][j];
      vtempz -= vz[onegsp][j];
      dum = (vtempx*vtempx +vtempy*vtempy +vtempz*vtempz);
      engy= Escale[onegsp]*dum;
      vel = sqrt(dum);
      k = x[otwosp][i];
      sigma_total = sigma_scale*nu_otwo_max/(vel*sp_n[otwosp][k]);
      random= frand();

      /***********************************************/
      /* O- + O2+ -> O + O2  (Mutual Neutralization) */
      
      if (random <= osigma20(engy)/sigma_total) {
	s = x[otwosp][i];
	del= x[otwosp][i] - s;
	
	/****  remove neg. ion ******/
	x[onegsp][j] =  x[onegsp][nnp+N-1];
	vx[onegsp][j]= vx[onegsp][nnp+N-1];
	vy[onegsp][j]= vy[onegsp][nnp+N-1];
	vz[onegsp][j]= vz[onegsp][nnp+N-1];
	
	/****  remove pos. ion ******/
	l = np[otwosp]-1;
	x[otwosp][i] =  x[otwosp][l];
	vx[otwosp][i]= vx[otwosp][l];
	vy[otwosp][i]= vy[otwosp][l];
	vz[otwosp][i]= vz[otwosp][l];
	sp_n_mcc[otwosp][s]  -= 1 -del;
	sp_n_mcc[otwosp][s+1]-= del;
	np[otwosp]--;
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  mccrate[19][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[19][s+1]+= (s== nc-1) ? 2*del : del;
	}
      	
	/** adjust the indecies and total numbers **/
	N--;  j--;  np[onegsp]--;
      }
    }
  }
  /******************************************************/
  /* 02+ + O2 -> ..........                              */
  
  if(isp==otwosp) {
    /**** First Clear the diagnostics  ******/
    if(theRunWithXFlag) {
      for (j=0; j<ng; j++)
	mccrate[20][j] = 0.0;
    }
    /**** Now do the collisions *****/
    otwo_extra += np[otwosp]*col_prob_otwo;
    N = otwo_extra;
    otwo_extra -= N;
    
    nnp = np[otwosp];
    for(j=0; j< N; j++)
    {
      index= nnp*frand();
      nnp--;
      temp = x[otwosp][nnp];
      x[otwosp][nnp] = x[otwosp][index];
      x[otwosp][index] = temp;
      
      temp = vx[otwosp][nnp];
      vx[otwosp][nnp] = vx[otwosp][index];
      vx[otwosp][index] = temp;
      
      temp = vy[otwosp][nnp];
      vy[otwosp][nnp] = vy[otwosp][index];
      vy[otwosp][index] = temp;
      
      temp = vz[otwosp][nnp];
      vz[otwosp][nnp] = vz[otwosp][index];
      vz[otwosp][index] = temp;
    }
    
    for(j=nnp; j<nnp+N; j++)
    {
      /** choose a neutral from bkgd **/
      maxwellv(&vneutx, &vneuty, &vneutz, vgth);

      /** transform into reference frame where v_neutral= 0 **/
      vx[otwosp][j] -= vneutx;
      vy[otwosp][j] -= vneuty;
      vz[otwosp][j] -= vneutz;
      dum = (vx[otwosp][j]*vx[otwosp][j] +vy[otwosp][j]*vy[otwosp][j]
	     +vz[otwosp][j]*vz[otwosp][j]);
      engy= Escale[otwosp]*dum;
      vel = sqrt(dum);
      sigma_total = max_sigmav_otwo/(vel*vscale);
      random= frand();
      
      /**************************/
      /* O2+ + O2 -> O2 + O2+   */
      
      if (random <= osigma18(engy)/sigma_total) {
	/* bring the ion to rest in the v_neutral frame */
	vx[otwosp][j]= vy[otwosp][j]= vz[otwosp][j]= 0.0;
	
	if(theRunWithXFlag) {
	  s= x[otwosp][j];
	  del= x[otwosp][j] - s;
	  mccrate[20][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[20][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
      /** transform back into lab frame **/
      vx[otwosp][j] += vneutx;
      vy[otwosp][j] += vneuty;
      vz[otwosp][j] += vneutz;
    }
  }
}

/**************************************************************/

void onewvel(SCALAR energy, SCALAR vel, SCALAR *vx, SCALAR *vy, SCALAR *vz, int e_flag)
{
  SCALAR phi1, cosphi, sinphi, coschi, sinchi, up1, up2, up3;
  SCALAR mag, r11, r12, r13, r21, r22, r23, r31, r32, r33;
  
  if(energy < 1e-30)  coschi = 1;
  else  coschi = (energy +2 -2*pow(energy +1,frand()))/energy;
  sinchi= sqrt(fabs(1. -coschi*coschi));
  
  phi1  = 2*M_PI*frand();
  cosphi= cos(phi1);
  sinphi= sin(phi1);
  
  if(e_flag) vel *= sqrt(1 - 2*m[ecolsp]*(1-coschi)/m[otwosp]);
  
  r13 = *vx;
  r23 = *vy;
  r33 = *vz;
  
  if(r33 == 1) { up1= 0;  up2= 1;  up3= 0; }
  else         { up1= 0;  up2= 0;  up3= 1; }
  
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

/**************************************************/
/* viberational excitation v=1 (E loss = 0.19 eV) */

SCALAR osigma1(SCALAR energy)
{
  SCALAR temp1, temp2, temp3;
  SCALAR temp4, temp5, temp6;

  temp6= 0.0;
  if(.19 <= energy && energy <= 1.0) {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;
    temp4 = energy*temp3;
    
    temp6 = -1.508395*temp4 +6.521786*temp3 -9.574636*temp2 +5.092031*temp1 -0.41602*energy -0.066398;
    temp6 *= 1e-20;
    if(temp6 <= 0.0) temp6=0.0;
  }
  else if(1.0 < energy && energy <= 1.67) temp6 = -7.2193e-22*(energy -1.67);
  else if(4.0 <= energy && energy <= 5.5) temp6 = 4.572852e-22*(energy -4.0);
  else if(5.5 < energy && energy <= 16.0)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;
    temp4 = energy*temp3;
    temp5 = energy*temp4;
    
    temp6 = 1.0835e-6*temp5 -9.229e-5*temp4 +0.0030853*temp3 -0.050981*temp2 +0.427934*temp1 -1.6682*energy +2.3919;
    temp6 *= 1e-20;
    if(temp6 <= 0.0) temp6=0.0;
  }
  else if(16.0 < energy && energy <= 25.0) temp6 = -4.098144e-23*(energy -25.0);

  return(temp6);
}

/**************************************************/
/* viberational excitation v=2 (E loss = 0.38 eV) */

SCALAR osigma2(SCALAR energy)
{
  SCALAR temp1, temp2, temp3;
  SCALAR temp4, temp5;

  temp5= 0.0;
  if(.38 <= energy && energy <= 1.67)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;
    
    temp5 = -0.606022*temp3 +3.157773*temp2 -5.933895*temp1 +4.664064*energy -1.233443;
    temp5 *= 1e-20;
    if(temp5 <= 0.0) temp5=0.0;
  }
  else if(4.0 <= energy && energy <= 14.5)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;
    temp4 = energy*temp3;
    
    temp5 = -3.1339358e-6*temp4 +2.0994236e-4*temp3 -.00503577*temp2 +0.0515*temp1 -0.2074798*energy +0.279;
    temp5 *= 1e-20;
    if(temp5 <= 0.0) temp5=0.0;
  }
  else if(14.5 < energy && energy <= 25.0) temp5 = -1.71326e-23*(energy -25.0);

  return(temp5);
}

/**************************************************/
/* viberational excitation v=3 (E loss = 0.57 eV) */

SCALAR osigma3(SCALAR energy)
{
  SCALAR temp1, temp2, temp3;
  SCALAR temp4, temp5;

  temp5= 0.0;
  if(.57 <= energy && energy <= 1.67)
  {
    temp1 = energy*energy;
    
    temp5 = -.055083*temp1 +.12457*energy -.057531;
    temp5 *= 1e-20;
    if(temp5 <= 0.0) temp5=0.0;
  }
  else if(4.0 <= energy && energy <= 15.0)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;
    temp4 = energy*temp3;
    
    temp5 = -7.969385e-6*temp4 +4.78119632e-4*temp3 -.0107124*temp2 +0.1095564*temp1 -0.4962553*energy +0.80444;
    temp5 *= 1e-20;
    if(temp5 <= 0.0) temp5=0.0;
  }
  else if(15.0 < energy && energy <= 20.0) temp5 = -1.76e-23*(energy -20.0);

  return(temp5);
}

/**************************************************/
/* viberational excitation v=4 (E loss = 0.75 eV) */

SCALAR osigma4(SCALAR energy)
{
  SCALAR temp1, temp2, temp3, temp6;

  temp6= 0.0;
  if(.75 <= energy && energy <= 0.85) temp6 = 2.795e-25*(energy - 0.75);
  
  else if(.85 <= energy && energy <= 1.67)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    
    temp6 = -0.049346*temp2 +0.16616*temp1 -0.174061*energy +0.058213;
    temp6 *= 1e-20;
    if(temp6 <= 0.0) temp6=0.0;
  }
  else if(6.0 <= energy && energy <= 15.0)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;
    
    temp6 = -1.3846154e-5*temp3 +8.8449e-4*temp2 -0.020271*temp1 +0.19111*energy -0.589505;
    temp6 *= 1e-20;
    if(temp6 <= 0.0) temp6=0.0;
  }
  return(temp6);
}

/**************************************************/
/*  Rotational excitation      (E loss = 0.02 eV) */

SCALAR osigma5(SCALAR energy)
{
  SCALAR temp1, temp2, temp3;
  SCALAR temp4, temp5, temp6;

  temp6= 0.0;
  if(.07 <= energy && energy <= 1.67)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;
    temp4 = energy*temp3;
    temp5 = energy*temp4;
    
    temp6 = -.0859*temp5 +.4233*temp4 -.7366*temp3 +.5205*temp2 -.1537*temp1 +.0604*energy -.0022;
    temp6 *= 1e-20;
  }
  return(temp6);
}

/**************************************************/
/*  momentum transfer                             */

SCALAR osigma6(SCALAR energy)
{
  SCALAR temp1, temp2, temp3;
  SCALAR temp4, temp5;

  if(0.0 <= energy && energy <= 1.2)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;
    
    temp5 = -12.6132*temp3 +39.2258*temp2 -43.3875*temp1 +23.5722*energy +.4464;
    temp5 *= 1e-20;
  }
  else if(1.2 < energy && energy <= 20.)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;
    temp4 = energy*temp3;
    
    temp5 = -4.0554e-5*temp4 +2.7604e-3*temp3 -.07107*temp2 +.82961*temp1 -3.9163*energy +11.735;
    temp5 *= 1e-20;
  }
  else if(20 < energy && energy <= 100.)
  {
    temp1 = energy*energy;
    temp5 = 1.3874e-4*temp1 -.0417*energy +9.254364;
    temp5 *= 1e-20;
  }
  else temp5= 6.5e-20;

  return(temp5);
}

/**************************************************/
/*  O2 SING DELTA             (E LOSS = 0.977 eV) */

SCALAR osigma7(SCALAR energy)
{
  SCALAR temp1, temp2, temp3;
  SCALAR temp4, temp5;

  temp5= 0.0;
  if(0.977 <= energy && energy <= 10.5)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;
    temp4 = energy*temp3;
    
    temp5 = -3.0913e-6*temp4 +1.436e-4*temp3 -.0022876*temp2 +.0133286*temp1 -.0100266*energy -.0015636;
    temp5 *= 1e-20;
    if(temp5 <= 0.0) temp5=0.0;
  }
  else if(10.5 < energy && energy <= 45.0)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;

    temp5 = -1.0959e-6*temp2 +1.349e-4*temp1 -.005984*energy +.1079;
    temp5 *= 1e-20;
  }
  else if(45 < energy && energy <= 100.0)
  {
    temp5 = -2.18e-4*energy +2.18e-2;
    temp5 *= 1e-20;
  }
  return(temp5);
}

/**************************************************/
/*  O2 B SINGLET SIGMA        (E LOSS = 1.627 eV) */

SCALAR osigma8(SCALAR energy)
{
  SCALAR temp1, temp2, temp3;
  SCALAR temp4, temp5;

  temp5= 0.0;
  if(1.627 <= energy && energy <= 25.0)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;
    temp4 = energy*temp3;
    
    temp5 = 5.34245e-8*temp4 -4.7117e-6*temp3 +1.581e-4*temp2 -2.4783e-3*temp1 +1.70373e-2*energy -2.2343e-2;
    temp5 *= 1e-20;
    if(temp5 <= 0.0) temp5=0.0;
  }
  else if(25.0 < energy && energy <= 45.0)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;

    temp5 = 5.445e-7*temp2 -5.674e-5*temp1 +1.7502e-3*energy -1.0375e-2;
    temp5 *= 1e-20;
  }
  else if(45 < energy && energy <= 100.0)
  {
    temp5 = -5.64e-5*energy +5.64e-3;
    temp5 *= 1e-20;
  }
  return(temp5);
}

/**************************************************/
/* 4.5 E loss                   (E LOSS = 4.5 eV) */

SCALAR osigma9(SCALAR energy)
{
  SCALAR temp1, temp2;
  SCALAR temp5;

  temp5= 0.0;
  if(4.5 <= energy && energy <= 4.8)
  {
    temp5 = 0.01*energy -.045;
    temp5 *= 1e-20;
  }
  else if(4.8 < energy && energy <= 15.0)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;

    temp5 = 6.113e-4*temp2 -2.11e-2*temp1 +.2216*energy -.638556;
    temp5 *= 1e-20;
  }
  return(temp5);
}

/**************************************************/
/* 6.0 E loss                   (E LOSS = 6.0 eV) */

SCALAR osigma10(SCALAR energy)
{
  SCALAR temp1, temp2, temp3;
  SCALAR temp5;

  temp5= 0.0;
  if(6.0<= energy && energy <= 20.)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;

    temp5 = -1.1894e-4*temp3 +6.8655e-3*temp2 -.143425*temp1 +1.26276*energy -3.7338513;
    temp5 *= 1e-20;
    if(temp5 <= 0.0) temp5=0.0;
  }
  else if(20. < energy && energy <= 100.)
  {
    temp1 = energy*energy;

    temp5 = 9.9341e-6*temp1 -1.7857e-3*energy +7.924e-2;
    temp5 *= 1e-20;
    if(temp5 <= 0.0) temp5=0.0;
  }
  return(temp5);
}

/**************************************************/
/* 8.4 E loss                   (E LOSS = 8.4 eV) */

SCALAR osigma11(SCALAR energy)
{
  SCALAR temp1, temp5;

  temp5= 0.0;
  if(8.4<= energy && energy <= 9.4)
  {
    temp5 = energy -8.4;
    temp5 *= 1e-20;
    if(temp5 <= 0.0) temp5=0.0;
  }
  else if(9.4 < energy && energy <= 100.)
  {
    temp1 = energy*energy;

    temp5 = -1.08852e-4*temp1 +1.10145e-2*energy +.92302246;
    temp5 *= 1e-20;
    if(temp5 <= 0.0) temp5=0.0;
  }
  else if(100. < energy) temp5 = 9.4e-21;
  return(temp5);
}

/**************************************************/
/* 9.97 E loss                  (E LOSS = 9.97 eV) */

SCALAR osigma12(SCALAR energy)
{
  SCALAR temp1, temp2, temp3;
  SCALAR temp5;

  temp5= 0.0;
  if(10. <= energy && energy <= 100.)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;

    temp5 = 1.1656e-9*temp3 -3.1555e-7*temp2 +1.8544e-5*temp1 +9.464e-4*energy -.0110422;
    temp5 *= 1e-20;
    if(temp5 <= 0.0) temp5=0.0;
  }
  else if(100. < energy) temp5= 7e-22;
  
  return(temp5);
}

/**************************************************/
/* Ionization                 (E LOSS = 12.06 eV) */

SCALAR osigma13(SCALAR energy)
{
  SCALAR temp1, temp2, temp3;
  SCALAR temp5;

  temp5= 0.0;
  if(12.06 <= energy && energy <= 100.)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;
    temp3 = energy*temp2;

    temp5 = 8.035e-8*temp3 -1.594e-5*temp2 +5.1392e-4*temp1 +.0658*energy -.89892;
    temp5 *= 1e-20;
    if(temp5 <= 0.0) temp5=0.0;
  }
  else if(100. < energy) temp5= 2.9e-20;
  
  return(temp5);
}

/**************************************************/
/* 130 nm line excitation      (E LOSS = 14.7 eV) */

SCALAR osigma14(SCALAR energy)
{
  SCALAR temp1, temp2;
  SCALAR temp5;

  temp5= 0.0;
  if(14.7<= energy && energy <= 100.)
  {
    temp1 = energy*energy;
    temp2 = energy*temp1;

    temp5 = 1.21665e-7*temp2 -3.0483e-5*temp1 +2.51713e-3*energy -.030335;
    temp5 *= 1e-20;
    if(temp5 <= 0.0) temp5=0.0;
  }
  else if(100. < energy) temp5= 3.8e-22;
  
  return(temp5);
}

/*******************************************************/
/*  e + O2 -> O + O-  Dissociative attachment (Et=4.2) */

SCALAR osigma15(SCALAR energy)
{
  register int i;
  SCALAR detach_sigma;
  static int init_flag=1;
  static SCALAR *detach1, *detach2;
  
  /********  Initialization  ********/
  if(init_flag) {
    SCALAR engy;
    detach1 = (SCALAR *) malloc(150*sizeof(SCALAR));
    detach2 = (SCALAR *) malloc(100*sizeof(SCALAR));
	 
    for(i=0; i<150; i++) {
      engy = 0.1*i;
      if(engy <= 7.0)
	detach1[i] = 1.4064e-22*exp(-(engy -6.5)*(engy -6.5)/1.1766);
      else if(7.0 <engy && engy <=15.0)
	detach1[i] = 9.0e-24 -6.0e-25*engy +1.4064e-22*exp(-(engy -6.5)*(engy -6.5)/1.1766);
    }
    for(i=0; i<100; i++) {
      engy = i;
      if(engy <15.0)
	detach2[i] = 0.0;
      else if(15.0<=engy && engy <31.0) 
	detach2[i] = 4.66e-23*exp(-(engy -30.0)*(engy -30.0)/50.0);
      else
	detach2[i] = 5.1745e-23 -1.96e-25*engy; 
    }
    init_flag=0;
  }
  
  /***************************/
  if(energy < 15.0) {
    i= energy*10 +0.5;
    detach_sigma = detach1[i];
  }
  else {
    i = energy;
    if(i >= 100) i = 99;
    detach_sigma = detach2[i];
  }
  return(detach_sigma);
}

/**************************************************/
/*  O- + O2 -> O2 + O + e    (Detachment)         */ 
    
SCALAR osigma16(SCALAR energy)
{
  if(energy < 1.465)                          return (0);
  else if(1.465 <=energy && energy < 28.626)  return (2.761e-21*(energy -1.465));
	return (7.5e-20);
}

/**************************************************/
/*  O- + O2 -> O- + O2      (Momentum Transfer)   */
    
SCALAR osigma17(SCALAR energy)
{
  if (energy <=0.2)  return(5e-19);
	return(1e-19 +1.79e-19/sqrt(energy));
}

/**************************************************/
/*  O2+ + O2 -> O2 + O2+    (Charge echange)      */
    
SCALAR osigma18(SCALAR energy)
{
 if (energy <=10.0)  return(1e-18);
 return(1e-19 +9e-18/energy);
}

/**************************************************/
/*  O + O2 -> O + O2        (Momentum Transfer)   */ 

SCALAR osigma19(SCALAR energy)
{
 if (energy <=0.2)  return(5e-19);
 return(1e-19 +1.79e-19/sqrt(energy));
}

/**************************************************/
/*  O- + O2+ -> O + O2    (Mutual Neutralization) */

SCALAR osigma20(SCALAR energy)
{
  return(2e-16*pow(10.0, -.0182*energy));
}

/**************************************************/
/* e + O2+ -> O + O  (Dissociative Recombination) */

SCALAR osigma21(SCALAR energy)
{
  return(1.8475e-19*pow(10.0, -.4524*energy));
}

/**************************************************/
/*  e + O- -> O + 2e    (Electron impact detach.) */

SCALAR osigma22(SCALAR energy)
{
  int i;
  static int init_flag=1;
  static SCALAR *edetach1, *edetach2;

  /********  Initialization  ********/
  if(init_flag) {
    SCALAR engy;
    edetach1 = (SCALAR *) malloc(300*sizeof(SCALAR));
    edetach2 = (SCALAR *) malloc(432*sizeof(SCALAR));
	 
    for(i=0; i<300; i++) {
      engy = 0.1*i;
      if(engy < 1.465) edetach1[i]= 0.0;
      else             edetach1[i]= -1.0353e-23*engy*engy +1.152e-21*engy +4.113e-20;
    }
    for(i=0; i<432; i++) {
      engy = i;
      if(engy <97.1) edetach2[i]= -1.0353e-23*engy*engy +1.152e-21*engy +4.113e-20;
      else           edetach2[i]= 1.175e-18*log(engy)/engy;
    }
    init_flag =0;
  }

  /***************************/

  if(energy < 30.) {
    i = 10*energy;
    return(edetach1[i]);
  }

	i = energy;
	if(i >= 432) i = 431;
	return(edetach2[i]);
 
}

/****************************************************************/

void omakethefile(void)
{
  SCALAR e;
  FILE *DMPFile;
  
  DMPFile = fopen("xsections", "w"); 
  
  for(e=0; e<3.0; e+= .01)
    fprintf(DMPFile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", e+1e-30,
	    osigma1(e)+1e-30, osigma2(e)+1e-30, osigma3(e)+1e-30, osigma4(e)+1e-30, osigma5(e)+1e-30,
	    osigma6(e)+1e-30, osigma7(e)+1e-30, osigma8(e)+1e-30, osigma9(e)+1e-30, osigma10(e)+1e-30,
	    osigma11(e)+1e-30, osigma12(e)+1e-30, osigma13(e)+1e-30, osigma14(e)+1e-30,
	    osigma15(e)+1e-30, osigma16(e)+1e-30, osigma17(e)+1e-30, osigma18(e)+1e-30,
	    osigma19(e)+1e-30, osigma20(e)+1e-30, osigma21(e)+1e-30, osigma22(e)+1e-30);
  
  for(e=3.0; e<100; e+=.1)
    fprintf(DMPFile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", e+1e-30,
	    osigma1(e)+1e-30, osigma2(e)+1e-30, osigma3(e)+1e-30, osigma4(e)+1e-30, osigma5(e)+1e-30,
	    osigma6(e)+1e-30, osigma7(e)+1e-30, osigma8(e)+1e-30, osigma9(e)+1e-30, osigma10(e)+1e-30,
	    osigma11(e)+1e-30, osigma12(e)+1e-30, osigma13(e)+1e-30, osigma14(e)+1e-30,
	    osigma15(e)+1e-30, osigma16(e)+1e-30, osigma17(e)+1e-30, osigma18(e)+1e-30,
	    osigma19(e)+1e-30, osigma20(e)+1e-30, osigma21(e)+1e-30, osigma22(e)+1e-30);
  
  fclose(DMPFile);
}

/**************************************************************/
