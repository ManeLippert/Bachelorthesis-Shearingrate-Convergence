#include "pdp1.h"
#include "xgrafix.h"

#define   NEMAX    500
SCALAR sigma1(SCALAR), sigma2(SCALAR), sigma3(SCALAR), sigma4(SCALAR), sigma5(SCALAR);
SCALAR gden, vgth;
int ecolsp, icolsp, ionsp;
void newvel(SCALAR, SCALAR, SCALAR *, SCALAR *, SCALAR *, int);

/***************************************************************/
/* Monte Carlo scheme for electron and ion-neutrals collisions */
/* easly convert old mcc cross sections */

#define selsmax 1.0e-19 
#define elsengy0 0.0
#define elsengy1 0.0
#define elsengy2 10.0
#define sextmax 1.0e-20 
#define extengy0 12.0    
#define extengy1 50.0
#define extengy2 100.0 
#define sionmax 1.0e-20
#define ionengy0 13.6 
#define ionengy1 60.0
#define ionengy2 110.0
#define achrgx 3.0e-19
#define bchrgx 0.0  
#define ascat 2.0e-19     
#define bscat 0.0

void mcc(int isp)
{
  static int init_flag=1;
  static SCALAR ecol_extra, icol_extra;
  static SCALAR max_sigmav_e, max_sigmav_i;
  static SCALAR col_prob_e, col_prob_i;
  
  int N, i, k, nnp, j, s, index;
  SCALAR random, dum, vel,sigma_total;
  SCALAR engy, rengy, del, phi1;
  SCALAR cosphi, sinphi, coschi, sinchi;
  SCALAR temp, vneutx, vneuty, vneutz;
  SCALAR tempsigma1, tempsigma2;
  SCALAR tempsigma5;
  
  if (init_flag)
		{
			gden = NperTORR*pressure/(gtemp+DBL_MIN);  /* Calculating the gas density */

			ecolsp= ecollisional-1;
			icolsp= icollisional-1;
			ionsp= ionspecies-1;          /* Fixing the indices into the array of species */
			
			if(ionsp>=0) vgth= sqrt(gtemp/Escale[ionsp]);

			/*****************************************/
			/* Calculating the null collision prob.  */
			
			if(ecollisional)
				{
					max_sigmav_e = 0.0;
					for (i=0; i< 10*NEMAX; i++)
						{
							engy = 0.1*i;
							max_sigmav_e= max(max_sigmav_e, sqrt(2.0*1.602e-19*engy/m[ecolsp])
																*(sigma1(engy)+sigma2(engy)+sigma3(engy)));
						}
					col_prob_e = 1 -exp(-max_sigmav_e*gden*sp_k[ecolsp]*dt);
				}
			
			if(icollisional)
				{
					max_sigmav_i = 0.0;
					for (i=0; i<10*NEMAX; i++)
						{
							engy = 0.1*i;
							max_sigmav_i= max(max_sigmav_i, sqrt(2.0*1.602e-19*engy/m[icolsp])
																*(sigma4(engy)+sigma5(engy)));
						}
					col_prob_i = 1 -exp(-max_sigmav_i*gden*sp_k[icolsp]*dt);
				}
			ecol_extra= icol_extra= 0.5;
			init_flag = 0;
		}
  
  /************************/
  /* Electron collisions  */
  
  if(ecollisional && isp==ecolsp)
		{
			/**** First Clear the diagnostics  ******/
			if(theRunWithXFlag) {
				for (j=0; j<ng; j++)
					mccrate[0][j] = mccrate[1][j] = mccrate[2][j] = 0.0;
			}

			ecol_extra += np[ecolsp]*col_prob_e;
			N = ecol_extra;
			ecol_extra -= N;
			
			nnp = np[ecolsp];
			for(j=0; j< N; j++)
				{
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
			
			for(j=nnp; j<nnp+N; j++)
				{
					/* determine if a collision occurs */
					dum = (vx[ecolsp][j]*vx[ecolsp][j] +vy[ecolsp][j]*vy[ecolsp][j]
								 +vz[ecolsp][j]*vz[ecolsp][j]);
					engy= Escale[ecolsp]*dum;
					vel = sqrt(dum);
//					sigma_total = sp_k[icolsp]*dt*max_sigmav_e/vel/dx;
					sigma_total = max_sigmav_e/(vel*vscale);
					random= frand();
					
					/*********************************************************************/
					/* determine the type of collision and calculate new velocities, etc */
					
					/*******************************/
					/* if the collision is elastic */
					
					if (random <= (tempsigma1 =sigma1(engy))/sigma_total)
						{
							/* first normalize vel */
							vx[ecolsp][j] /= vel;
							vy[ecolsp][j] /= vel;
							vz[ecolsp][j] /= vel;
							
							/* scatter the electron */
							newvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 1);
							/* collect diagnostics*/
							if(theRunWithXFlag) {
								s= x[ecolsp][j];
								del= x[ecolsp][j] - s;
								mccrate[0][s]  += (!s) ? 2*(1-del) : 1-del;
								mccrate[0][s+1]+= (s== nc-1) ? 2*del : del;
							}
							
						}
					/**********************************/
					/* if the collision is excitation */
					
					else if (random <= (tempsigma1 +(tempsigma2 =sigma2(engy)))/sigma_total)
						{
							/* first normalize vel */
							vx[ecolsp][j] /= vel;
							vy[ecolsp][j] /= vel;
							vz[ecolsp][j] /= vel;
							
							engy -= extengy0;
							vel   = sqrt(fabs(engy)/Escale[ecolsp]);
							
							/* scatter the electron */
							newvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
							
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
					
					else if(random <= (tempsigma1 + tempsigma2 +sigma3(engy))/sigma_total)
						{
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
							newvel(rengy, vel, &vx[ecolsp][k], &vy[ecolsp][k], &vz[ecolsp][k], 0);
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
							newvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);
							
							if(++np[ionsp] >maxnp[ionsp] || ++np[ecolsp] >maxnp[ecolsp])
								puts("mcc(Ionization): too many particles. MUST EXIT!");

							/* collect diagnostics */
							if(theRunWithXFlag) {
								mccrate[2][s]  += (!s) ? 2*(1-del) : 1-del;
								mccrate[2][s+1]+= (s== nc-1) ? 2*del : del;
							}
						}
				}
		}

  /**************************************/

  if(icollisional && isp==icolsp){
		/**** First Clear the diagnostics  ******/
		if(theRunWithXFlag) {
			for (j=0; j<ng; j++) mccrate[3][j] = mccrate[4][j] =0.0;
		}
		/**** Now do the collisions *****/
		icol_extra += np[icolsp]*col_prob_i;
		N = icol_extra;
		icol_extra -= N;
		
		nnp = np[icolsp];
		for(j=0; j< N; j++)
			{
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
		
		for(j=nnp; j<nnp+N; j++)
			{
				maxwellv(&vneutx, &vneuty, &vneutz, vgth);
				vx[icolsp][j] -= vneutx;
				vy[icolsp][j] -= vneuty;
				vz[icolsp][j] -= vneutz;
				dum = (vx[icolsp][j]*vx[icolsp][j] +vy[icolsp][j]*vy[icolsp][j]
							 +vz[icolsp][j]*vz[icolsp][j]);
				if (dum){
					engy= Escale[icolsp]*dum;
					vel = sqrt(dum);
					//				sigma_total = dt*max_sigmav_i/vel/dx;
					sigma_total = max_sigmav_i/(vel*vscale);
					random= frand();
				
				/**********************************/
				/* if the collision is scattering */
				
				if (random <= (tempsigma5 =sigma5(engy))/sigma_total)
					{
						newivel(engy, vel, &vx[icolsp][j], &vy[icolsp][j], &vz[icolsp][j]);
						/* collect diagnostics */
						if(theRunWithXFlag) {
							s= x[icolsp][j];
							del= x[icolsp][j] - s;
							mccrate[3][s]  += (!s) ? 2*(1-del) : 1-del;
							mccrate[3][s+1]+= (s== nc-1) ? 2*del : del;
						}
					}
				/***************************************/
				/* if the collision is charge exchange */
				
				else if (random <= (tempsigma5 +sigma4(engy))/sigma_total)
					{
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

void newvel(SCALAR energy, SCALAR vel, SCALAR *vx, SCALAR *vy, SCALAR *vz, int e_flag)
{
  SCALAR phi1, cosphi, sinphi, coschi, sinchi, up1, up2, up3;
  SCALAR mag, r11, r12, r13, r21, r22, r23, r31, r32, r33;

  if(energy < 1e-30)  coschi = 1;
  else  coschi = (energy +2 -2*pow(energy +1,frand()))/energy;
  sinchi= sqrt(1. - coschi*coschi);
  
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
/* XSection for Elastic Collisions */

SCALAR sigma1(SCALAR energy)
{
  int i;
  static SCALAR *sels;
  static int init_flag=1;
  
  /******************  Initialization   ********************/
  if(init_flag)
		{
			sels= (SCALAR *)malloc(NEMAX*sizeof(SCALAR));
			for (i=0; i<NEMAX; i++)
				{
					if (i <= elsengy0) sels[i] = 0;
					else if (i > elsengy0 && i <= elsengy1)
						sels[i] = selsmax*(i - elsengy0)/(elsengy1 - elsengy0);
					else if (i > elsengy1 && i <= elsengy2)
						sels[i] = selsmax;
					else					
						sels[i] = selsmax*elsengy2/log(elsengy2)*log((SCALAR)i)/i;
				}
			init_flag =0;
		}
  /**************************************/
  
  i= energy +.5;
  if(NEMAX <= i) i= NEMAX-1;
  return(sels[i]);
}

/**************************************/
/* XSection for Excitation Collisions */

SCALAR sigma2(SCALAR energy)
{
  int i;
  static SCALAR *sext;
  static int init_flag=1;
  
  /******************  Initialization   ********************/
  if(init_flag)
		{
			sext= (SCALAR *)malloc(NEMAX*sizeof(SCALAR));
			for (i=0; i<NEMAX; i++)
				{
					if (i <= extengy0) sext[i] = 0;
					else if (i > extengy0 && i <= extengy1)
						sext[i] = sextmax*(i - extengy0)/(extengy1 - extengy0);
					else if (i > extengy1 && i <= extengy2)
						sext[i] = sextmax;
					else					
						sext[i] = sextmax*extengy2/log(extengy2)*log((SCALAR)i)/i;
				}
			init_flag =0;
		}
  /**************************************/

  i= energy +.5;
  if(NEMAX <= i) i= NEMAX-1;
  return(sext[i]);
}

/**************************************/
/* XSection for Ionization Collisions */

SCALAR sigma3(SCALAR energy)
{
  int i;
  static SCALAR *sion;
  static int init_flag=1;
  
  /******************  Initialization   ********************/
  if(init_flag)
		{
			sion= (SCALAR *)malloc(NEMAX*sizeof(SCALAR));
			for (i=0; i<NEMAX; i++)
				{
					if (i <= ionengy0) sion[i] = 0;
					else if (i > ionengy0 && i <= ionengy1)
						sion[i] = sionmax*(i - ionengy0)/(ionengy1 - ionengy0);
					else if (i > ionengy1 && i <= ionengy2)
						sion[i] = sionmax;
					else					
						sion[i] = sionmax*ionengy2/log(ionengy2)*log((SCALAR)i)/i;
				}
			init_flag =0;
		}
  /**************************************/
  
  i= energy +.5;
  if(NEMAX <= i) i= NEMAX-1;
  return(sion[i]);
}

/*******************************************/
/* XSection for Charge Exchange Collisions */

SCALAR sigma4(SCALAR energy)
{
  return(achrgx + bchrgx/(sqrt(energy) +1e-30));
}

/**************************************/
/* XSection for Scattering Collisions */

SCALAR sigma5(SCALAR energy)
{
  return(ascat + bscat/(sqrt(energy) +1e-30));
}

/**************************************************************/

void mmakethefile(void)
{
  int i;
  SCALAR e;
  FILE *DMPFile;

  DMPFile = fopen("xsections", "w"); 

  for(i=0; i<100; i++)
		{
			e = .01*i;
			fprintf(DMPFile, "%f %e %e %e %e %e\n", e, sigma1(e),
							sigma2(e), sigma3(e), sigma4(e), sigma5(e));
		}
  for(i=10; i<200; i++)
		{
			e = .1*i;
			fprintf(DMPFile, "%f %e %e %e %e %e\n", e, sigma1(e),
							sigma2(e), sigma3(e), sigma4(e), sigma5(e));
		}
  for(i=20; i<600; i++)
		{
			e = i;
			fprintf(DMPFile, "%f %e %e %e %e %e\n", e, sigma1(e),
							sigma2(e), sigma3(e), sigma4(e), sigma5(e));
		}
  fclose(DMPFile);
}

/* mindgame: common function to update ion scattered velocity which is used in mcc.c, argonmcc.c, and heliummcc.c */
void newivel(SCALAR energy, SCALAR vel, SCALAR *vx, SCALAR *vy, SCALAR *vz)
{
	SCALAR up1, up2, up3, mag;
	SCALAR r11, r12, r13, r21, r22, r23, r31, r32, r33;
	SCALAR phi1, cosphi, sinphi, coschi, sinchi;

	coschi= sqrt(frand());
	sinchi= sqrt(fabs(1.0 -coschi*coschi));
	
	phi1  = 2*M_PI*frand();
	cosphi= cos(phi1);
	sinphi= sin(phi1);
	
	r13 = *vx/vel;
	r23 = *vy/vel;
	r33 = *vz/vel;
	
	if(r33 == 1.0) { up1= 0;  up2= 1;  up3= 0; }
	else           { up1= 0;  up2= 0;  up3= 1; }
	
	r12 = r23*up3 -r33*up2;
	r22 = r33*up1 -r13*up3;
	r32 = r13*up2 -r23*up1;
	mag = sqrt(r12*r12 + r22*r22 + r32*r32) + DBL_MIN;
	
	r12/= mag;
	r22/= mag;
	r32/= mag;
	
	r11 = r22*r33 -r32*r23;
	r21 = r32*r13 -r12*r33;
	r31 = r12*r23 -r22*r13;
	
	vel*= coschi;
	/* do not multiply coschi twice */
	*vx= vel*(r11*sinchi*cosphi +r12*sinchi*sinphi +r13*coschi);
	*vy= vel*(r21*sinchi*cosphi +r22*sinchi*sinphi +r23*coschi); 
	*vz= vel*(r31*sinchi*cosphi +r32*sinchi*sinphi +r33*coschi);
}

