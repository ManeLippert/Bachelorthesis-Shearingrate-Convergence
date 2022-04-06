#include "pdp1.h"
#include "xsect.h"
#include "xgrafix.h"
#include <string.h>    

#define   NEMAX    150
#define   NIMAX    10

SCALAR nsigma0(SCALAR), nsigma1(SCALAR), nsigma2(SCALAR), nsigma3(SCALAR);
SCALAR nsigma4(SCALAR), nsigma5(SCALAR), nsigma6(SCALAR);
SCALAR thrsh_engy1, thrsh_engy2, thrsh_engy3;
SCALAR thrsh_engy4, thrsh_engy5;
SCALAR gden, vgth;
int ecolsp, icolsp, ionsp;
int nTables;           /* number of cross section tables read    */
CrossSect *xSectTable; /* array of cross section tables for Neon */
CrossSect eTotalMom;

void nnewvel(SCALAR energy, SCALAR vel, SCALAR *vx, SCALAR *vy, SCALAR *vz, SCALAR xi, int e_flag);

/*********************************************************/
/* Monte Carlo Collisions for electron- and ion-neutrals */

void neonmcc(int isp)
{
  static int init_flag=1;
  static SCALAR ecol_extra, icol_extra;
  static SCALAR max_sigmav_e, max_sigmav_i;
  static SCALAR col_prob_e, col_prob_i;

  register int i, j, k, index;
  int N, nnp, s;
  SCALAR random, dum, vel, sigma_total;
  SCALAR engy, rengy, del;
  SCALAR sum_sigma, temp, vneutx, vneuty, vneutz;
  char filename[50];
  
  if (init_flag) {
    gden = NperTORR*pressure/(gtemp+DBL_MIN);  /* Calculating the gas density */
    
    printf("gden = %lf \n", gden);
    
    ecolsp= ecollisional-1;
    icolsp= icollisional-1;
    ionsp = ionspecies-1;     /* Fixing the indices into the array of species */
    if(ionsp>=0) vgth= sqrt(gtemp/Escale[ionsp]);
    
    /*****************************************/
    /* Reading the cross section data tables */

    strcpy(filename, PDP1PATH);
    strcat(filename, "/neon.tbl");
    nTables = readXSectionTables(&xSectTable, &eTotalMom, filename);
    
    /*****************************************/
    /* Calculating the null collision prob.  */
    
    if(ecollisional) {
      thrsh_engy1= xSectTable[1].thresholdE;   /*  e + Ne -> e + Ne*   Excitation  */
      thrsh_engy2= xSectTable[2].thresholdE;   /*  e + Ne -> e + Ne**  Excitation  */
      thrsh_engy3= xSectTable[3].thresholdE;   /*  e + Ne -> e + Ne  18.54 E Loss  */
      thrsh_engy4= xSectTable[4].thresholdE;   /*  e + Ne -> e + Ne  18.95 E Loss  */
      thrsh_engy5= xSectTable[5].thresholdE;   /*  e + Ne -> e + e + Ne+  Ioniz'n  */
      
      max_sigmav_e = 0.0;
      for (i = 0; i < eTotalMom.n; i++)
	max_sigmav_e = max(max_sigmav_e, eTotalMom.sigma[i]*sqrt(eTotalMom.E[i]));
      max_sigmav_e *= sqrt(2*1.602E-19/m[isp]);
      col_prob_e = 1 -exp(-max_sigmav_e*gden*dt*sp_k[ecolsp]);
      ecol_extra = 0.5;
      
      printf("e collision prob = %lf, max(sigma-v) = %lf, nsigmavdt= %lf\n",
	     col_prob_e, max_sigmav_e, max_sigmav_e*gden*dt*sp_k[ecolsp]);
    }
    if(icollisional) {
      max_sigmav_i = 0.0;
      for (engy=0; engy<NIMAX; engy += 0.1)
	max_sigmav_i= max(max_sigmav_i,
			  sqrt(2*1.602e-19*engy/m[icolsp])*nsigma6(engy));
      
      col_prob_i = 1 -exp(-max_sigmav_i*gden*dt*sp_k[icolsp]);
      icol_extra = 0.5;
      
      printf("i collision prob = %lf, max(sigma-v) = %lf, nsigmavdt= %lf\n",
	     col_prob_i, max_sigmav_i, max_sigmav_i*gden*dt*sp_k[icolsp]);
    }
    init_flag = 0;
  }
  
  /********************************/
  /* Electron collisions with Ne  */
  
  if(ecollisional && isp==ecolsp) {
    /**** First Clear the diagnostics  ******/
    if(theRunWithXFlag) {
      for (j=0; j<ng; j++) {
	mccrate[0][j] = mccrate[1][j] = mccrate[2][j] = 0.0;
	mccrate[3][j] = mccrate[4][j] = mccrate[5][j] = 0.0;
      }
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
      
      /************************************/
      /* e + Ne -> e + Ne  Mom. Trans.    */
      
      if (random <= (sum_sigma =nsigma0(engy))/sigma_total) { 

	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	/* scatter the electron */
	nnewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], x[ecolsp][j], 1);
	
	/* collect diagnostics*/
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[0][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[0][s+1]+= (s== nc-1) ? 2*del : del;
	}
      } 
      /************************************/
      /*  e + Ne -> e + Ne*   Excitation  */
      
      else if (engy>= thrsh_engy1 && random <= (sum_sigma +=nsigma1(engy))/sigma_total) {

	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= thrsh_engy1;
	vel   = sqrt(fabs(engy)/Escale[ecolsp]);
	
	/* scatter the electron */
	nnewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], x[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[1][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[1][s+1]+= (s== nc-1) ? 2*del : del;
	}
      } 
      /************************************/
      /*  e + Ne -> e + Ne**  Excitation  */
      
      else if (engy>= thrsh_engy2 && random <= (sum_sigma +=nsigma2(engy))/sigma_total) {

	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= thrsh_engy2;
	vel   = sqrt(fabs(engy)/Escale[ecolsp]);
	
	/* scatter the electron */
	nnewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], x[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[2][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[2][s+1]+= (s== nc-1) ? 2*del : del;
	}
      } 
      /************************************/
      /*  e + Ne -> e + Ne 18.54 E Loss   */
      
      else if (engy>= thrsh_engy3 && random <= (sum_sigma +=nsigma3(engy))/sigma_total) {

	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= thrsh_engy3;
	vel   = sqrt(fabs(engy)/Escale[ecolsp]);
	
	/* scatter the electron */
	nnewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], x[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[3][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[3][s+1]+= (s== nc-1) ? 2*del : del;
	}
      } 
      /************************************/
      /*  e + Ne -> e + Ne 18.95 E Loss   */
      
      else if (engy>= thrsh_engy4 && random <= (sum_sigma +=nsigma4(engy))/sigma_total) {

	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	engy -= thrsh_engy4;
	vel   = sqrt(fabs(engy)/Escale[ecolsp]);
	
	/* scatter the electron */
	nnewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], x[ecolsp][j], 0);
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  s= x[ecolsp][j];
	  del= x[ecolsp][j] - s;
	  mccrate[4][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[4][s+1]+= (s== nc-1) ? 2*del : del;
	}
      } 
      /************************************/
      /*  e + Ne -> e + e + Ne+  Ioniz'n  */
    
      else if(engy >= thrsh_engy5 && random <= (sum_sigma +=nsigma5(engy))/sigma_total) {

	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	/********************************/
	/* subtract the ion. energy and */
	/* partition the remaining energy */
	
	engy -= thrsh_engy5;
	rengy = 10.0*tan(frand()*atan(engy/20.0));
	engy -= rengy;
	
	/********************************/
	/* scatter the created electron */

	vel = sqrt(fabs(rengy)/Escale[ecolsp]);
	k = np[ecolsp];
	vx[ecolsp][k] = vx[ecolsp][j];
	vy[ecolsp][k] = vy[ecolsp][j];
	vz[ecolsp][k] = vz[ecolsp][j];
	nnewvel(rengy, vel, &vx[ecolsp][k], &vy[ecolsp][k], &vz[ecolsp][k], x[ecolsp][j], 0);
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
	nnewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], x[ecolsp][j], 0);

	if(++np[ionsp] >maxnp[ionsp] || ++np[ecolsp] >maxnp[ecolsp]) {
	  puts("ADJUST(Ionization): too many particles. MUST EXIT!");
	  exit(-1);
	}
	
	/* collect diagnostics */
	if(theRunWithXFlag) {
	  mccrate[5][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[5][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
    } 
  }
  
  /**************************************/
  /* Ne+ + Ne -> .....                  */
  
  if(icollisional && isp==icolsp) {
    /**** First Clear the diagnostics  ******/
    if(theRunWithXFlag) {
      for (j=0; j<ng; j++) mccrate[6][j] = 0.0;
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
    
    for (j=nnp; j<nnp+N; j++) { 
      maxwellv(&vneutx, &vneuty, &vneutz, vgth);
      vx[icolsp][j] -= vneutx;
      vy[icolsp][j] -= vneuty;
      vz[icolsp][j] -= vneutz;
      dum = (vx[icolsp][j]*vx[icolsp][j] +vy[icolsp][j]*vy[icolsp][j]
				 +vz[icolsp][j]*vz[icolsp][j]);
			if (dum){
				engy= Escale[icolsp]*dum;
				vel = sqrt(dum);
				sigma_total = max_sigmav_i/(vel*vscale+1e-30);
				random= frand();
			}
      else { /* did collide*/
				sigma_total = 1e7;
				random = 0;
			}

      /************************************/
      /*  Ne + Ne+ -> Ne+ + Ne  Charge X  */
      
      if (random <= nsigma6(engy)/sigma_total) {
				vx[icolsp][j] = vy[icolsp][j] = vz[icolsp][j] = 0.0;
				
				/* collect diagnostics */
				if(theRunWithXFlag) {
					s= x[icolsp][j];
					del= x[icolsp][j] - s;
					mccrate[6][s]  += (!s) ? 2*(1-del) : 1-del;
					mccrate[6][s+1]+= (s== nc-1) ? 2*del : del;
				}
      }
      vx[icolsp][j] += vneutx;
      vy[icolsp][j] += vneuty;
      vz[icolsp][j] += vneutz;
    }
  }
}
/**************************************************************/

void nnewvel(SCALAR energy, SCALAR vel, SCALAR *vx, SCALAR *vy, SCALAR *vz, SCALAR xi, int e_flag)
{
  SCALAR phi1, cosphi, sinphi, coschi, sinchi, up1, up2, up3;
  SCALAR mag, r11, r12, r13, r21, r22, r23, r31, r32, r33;
  
  /* if(energy < 1e-30)  coschi = 1;
  else  coschi = (energy +2.0 -2.0*pow(energy+1.0, frand()))/energy; */

  coschi = 1 -2.0*frand(); /* Isotropic */
  
  sinchi= sqrt(fabs(1.0 -coschi*coschi));
  
  phi1  = 2*M_PI*base2();
  cosphi= cos(phi1);
  sinphi= sin(phi1);
  
  if(e_flag) vel *= sqrt(1.0 -2.0*m[ecolsp]*(1-coschi)/m[ionsp]);
  
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
/*  e + Ne -> e + Ne  Mom. Trans.  */

SCALAR nsigma0(SCALAR energy)
{
  return(interpolateTableLinear(&xSectTable[0], energy));
}

/************************************/
/*  e + Ne -> e + Ne*   Excitation  */

SCALAR nsigma1(SCALAR energy)
{
  return(interpolateTableLinear(&xSectTable[1], energy));
}

/************************************/
/*  e + Ne -> e + Ne**  Excitation  */

SCALAR nsigma2(SCALAR energy)
{
  return(interpolateTableLinear(&xSectTable[2], energy));
}

/************************************/
/*  e + Ne -> e + Ne 18.54 E Loss   */

SCALAR nsigma3(SCALAR energy)
{
  return(interpolateTableLinear(&xSectTable[3], energy));
}

/************************************/
/*  e + Ne -> e + Ne 18.95 E Loss   */

SCALAR nsigma4(SCALAR energy)
{
  return(interpolateTableLinear(&xSectTable[4], energy));
}

/************************************/
/*  e + Ne -> e + e + Ne+  Ion.     */

SCALAR nsigma5(SCALAR energy)
{
  return(interpolateTableLinear(&xSectTable[5], energy));
}

/************************************/
/*  Ne + Ne+ -> Ne+ + Ne  Charge X  */

SCALAR nsigma6(SCALAR energy)
{
  return(6.321e-20/(sqrt(energy)+1e-30));
}

/************************************/
/*  Ne + Ne+ -> Ne+ + Ne  Charge X  */

SCALAR old_sigma6(SCALAR energy)
{
  SCALAR b2, loge, loge2, loge3, loge4, loge5;
  
  if(energy<6.3) return(3.134464e-22);

	loge  = log(sqrt(energy/27.2));
	loge2 = loge*loge;
	loge3 = loge*loge2;
	loge4 = loge*loge3;
	loge5 = loge*loge4;
	b2 = 43.1281 +21.887*loge -35.4268*loge2 +16.9951*loge3 -3.67456*loge4 +0.294592*loge5; 
	
	return(4.39867e-21*b2);
  
}

/**************************************************************/

void nmakethefile(void)
{
  SCALAR e;
  FILE *DMPFile;
  
  DMPFile = fopen("neonxsects.dat", "w"); 
  
  for(e=0.1; e<100; e+=0.1)
    fprintf(DMPFile, "%e %e %e %e %e %e %e %e\n", e+1e-30, nsigma0(e)+1e-30,
	    nsigma1(e)+1e-30, nsigma2(e)+1e-30, nsigma3(e)+1e-30,
	    nsigma4(e)+1e-30, nsigma5(e)+1e-30, nsigma6(e)+1e-30);
  
  fclose(DMPFile);
}

/**************************************************************/
