#include "pdp1.h"

void realft(SCALAR data[], int, int);
void velocity(void), freqanalysis(void);
SCALAR bit_rever_maxwell(), maxwell(void);

/**********************************************************************/

double base2()
{
  static unsigned int counter2=0;
  return(revers_base(counter2++, 2));
}

/**********************************************************************/

double revers_base(int num, int n)
{ 
  double power, rev;
  unsigned int inum, iquot, irem;
  
  rev = 0.;
  inum = num;
  power = 1.;
  
  do {
    iquot = inum/n;
    irem = inum - n*iquot;
    power /= n;
    rev += irem*power;
    inum = iquot;
  } while (inum > 0);
  
  return (rev);
}

/************************************************************************/
/* irand(iseed) returns values 1 through 2**31-2.                       */
/* From "Random number generators: good ones are hard to find", S. Park */
/* and K. Miller, Communications of ACM, October 1988, pp 1192-1201.    */
/* This is from page 1195, and is to work on any system for which       */
/* maxint is 2**31-1 or larger. Due earlier to Schrage, as cited by P&M.*/
/*                                                                      */
/* Note: OK values for iseed are 1 through 2**31-2. Give it 0 or 2*31-1 */
/* and it will return the same values thereafter!                       */
/*                                                                      */
/* C version 6/91, Bruce Langdon.                                       */
/*                                                                      */
/* Algorithm replaces seed by mod(a*seed,m). First represent            */
/* seed = q*hi + lo.  Then                                              */
/* a*seed = a*q*hi + lo = (m - r)*hi + a*lo = (a*lo - r*hi) + m*hi,     */
/* and new seed = a*lo - r*hi unless negative; if so, then add m.       */

SCALAR frand(void)
{
  long af = 16807, mf = 2147483647, qf = 127773, r = 2836;
  long hi, lo;
  SCALAR fnumb;
  
  hi = seed/qf;
  lo = seed - qf*hi;
  seed = af*lo - r*hi;
  /* "seed" will always be a legal integer of 32 bits (including sign). */
  if(seed <= 0) seed = seed + mf;
  fnumb = seed/2147483646.0;

  return(fnumb);
}

/**********************************************************************/

SCALAR bit_rever_maxwell()
{
	static int iset=0;
	static int count=0;
	static SCALAR gset;
	SCALAR fac,r,v1,v2;

	if(iset==0){
		do{
			v1=2.0*revers_base(count++,2)-1.0;
			v2=2.0*revers_base(count,5)-1.0;
			r = v1*v1+v2*v2;
		} while (r>=1.0);
		fac=sqrt(-2*log(r)/r);
		gset=v1*fac;
		iset = 1;
		return v2*fac;
	}
	else{
		iset=0;
		return gset;
	}
}

/**********************************************************************/

SCALAR maxwell()
{
	static int iset=0;
	static SCALAR gset;
	SCALAR fac,r2,v1,v2;

	if(iset==0){
		do{
			v1=2.0*frand()-1.0;
			v2=2.0*frand()-1.0;
			r2 = v1*v1+v2*v2;
		} while (r2>=1.0);
		fac=sqrt(-2*log(r2)/r2);
		gset=v1*fac;
		iset = 1;
		return v2*fac;
	}
	else{
		iset=0;
		return gset;
	}
}

/**********************************************************************/

SCALAR maxwellian(int type)
{
	switch (type){
	case RANDOM:
		return maxwell();
		break;
	case BIT_REVERSE:
		return bit_rever_maxwell();
		break;
	default:
		puts("LOAD: Bad value for loader flag");
		exit(-1);
		break;
	}
}

/**************************************************************/

#define  NVTC    3.3 // number of thermal velocity for upper cutoff

SCALAR distribution(int k, int isp, int type)
{
	static int init=1;
	static int cutoff=1;
	static int init_sp[NSMAX][2];
	static int nvel[NSMAX][2];
	static SCALAR *v[NSMAX][2];
	SCALAR erfvl, diff_erf, vlower, vupper, flower, fupper, f, fmid;
	SCALAR dv, dvi;
	int i, n, ksign;
	SCALAR vnew;

	if(init){
		for(i=0; i<NSMAX; i++)
			for(ksign=0; ksign<2; ksign++)
				init_sp[i][ksign] = 1;
		init = 0;
	}
	
	if (init_sp[isp][k]) {
		init_sp[isp][k] = 0;
		if (vc[isp][k]!=0){
			erfvl = erf(vc[isp][k]);
			diff_erf = erf((SCALAR)NVTC)-erfvl;
			nvel[isp][k] = 1.0/(1-diff_erf/(1-erfvl));
			v[isp][k] = (SCALAR *) malloc(nvel[isp][k]*sizeof(SCALAR));

			v[isp][k][0] = vc[isp][k];
			v[isp][k][nvel[isp][k]-1] = NVTC;
			dvi = (NVTC-vc[isp][k])/((SCALAR)nvel[isp][k]);
			dv = dvi;
			fmid = 0;
			vlower = vc[isp][k];
			f = 0;
			flower = erf(vlower)-f*(diff_erf)-erfvl;
			vupper = vlower + dv;
			for (n=1; n<(nvel[isp][k]-1); n++){
				f = (SCALAR)n/(SCALAR)(nvel[isp][k]-1);
				flower = erf(vlower) - f*diff_erf - erfvl;
				fupper = erf(vupper) - f*diff_erf - erfvl;
				while ((flower*fupper)>0.0){
					dv*=2;
					vupper = vlower + dv;
					fupper = erf(vupper) - f*diff_erf - erfvl;
				}
				while (dv>1e-8){
					dv /=2;
					fmid = erf(vlower+dv) - f*diff_erf - erfvl; 
					if (fmid<=0) 
						vlower+=dv;
					else
						vupper-=dv;
				}
				v[isp][k][n] = v0[isp][k]+vt[isp][k]*vlower;
			}
		}
		else
			cutoff = 0.0;
	}
	
	if (vt[isp][k]==0)
		vnew = v0[isp][k];
	else if (cutoff){
		switch (type){
		case RANDOM:
			vnew = v[isp][k][(int)(nvel[isp][k]*frand())];
			break;
		case BIT_REVERSE:
			vnew = v[isp][k][(int)(nvel[isp][k]*base2())];
			break;
    default:
      puts("LOAD: Bad value for loader flag");
      exit(-1);
      break;
		}
	}
	else
		vnew = v0[isp][k]+vt[isp][k]*fabs(maxwellian(type));

	return vnew;			
		
}

/**************************************************************/
/*
SCALAR maxwellian_flux(void)
{
	return sqrt(-log(frand()));
}
*/

SCALAR distribution_flux(int k, int isp, int type)
{
	static int init=1;
	static int cutoff=1;
	static int init_sp[NSMAX][2];
	static SCALAR vc2[NSMAX][2];
	SCALAR v;
	int i,ksign;

	if(init){
		for(i=0; i<NSMAX; i++)
			for(ksign=0; ksign<2; ksign++)
				init_sp[i][ksign] = 1;
		init = 0;
	}
	if (init_sp[isp][k]){
		if (vc[isp][k]==0)
			cutoff = 0;
		else
			vc2[isp][k] = vc[isp][k]*vc[isp][k];
		init_sp[isp][k] = 0;
	}
	
	if (!cutoff)
		v = v0[isp][k]+vt[isp][k]*sqrt(-2*log(frand()));
	else
		v = v0[isp][k]+vt[isp][k]*sqrt(vc2[isp][k]-log(frand()));
	
	return v;
}
/************************************************************************/
/* Returns bit reversed num in base 2 */

double revers(unsigned int num)
{
  double f=0.5, sum=0.0;
  
  while(num)
		{
			if (num & 1) sum += f;     /* is 1st bit set? */
			f *= 0.5;
			num >>= 1;		       /* fast divide by 2 */
		}
  return (sum);
}

/**************************************************/
/* bulk plasma Efield */
SCALAR E_p(int N)
{
  int i, j;
  SCALAR v;

  j = nc/2;
  j++;
  N = 30;
  v = e[j];

  for(i=1; i <= N; i++)
    {
      v += (e[j+i] + e[j-i]);
    }

  v /= (SCALAR)(2*N + 1);
  return v;
}

/**************************************************/
/* left sheath width */
SCALAR s_1(void)
{
  SCALAR eprime, E;
  int i;
  int ncoeff = 3;
  int coeff[3];
  coeff[0] = -3.0;
  coeff[1] = 4.0;
  coeff[2] = -1.0;

  for(i=0, eprime=0.0; i < ncoeff; i++)
    {
      eprime += coeff[i]*e[i];
    }
  eprime /= (2.0*dx);

  /*  fprintf(stderr, "\neprime = %f, Ep=%f", eprime, E_p()); */

  E = E_p((int)(nc/4));

  eprime =  (E - e[0])/eprime;
  if(eprime < 0.0) { eprime = 0.0; }

  return eprime;
}

/**************************************************/
/* right sheath width */
SCALAR s_2(void)
{
  SCALAR eprime, E;
  int i;
  int ncoeff = 3;
  int coeff[3];
  coeff[0] = -3.0;
  coeff[1] = 4.0;
  coeff[2] = -1.0;

  for(i=0, eprime=0.0; i < ncoeff; i++)
    {
      eprime -= coeff[i]*e[ng-1-i];
    }
  eprime /= (2.0*dx);
  
  E = E_p((int)(nc/4));

  eprime = (e[ng-1] - E)/eprime;
  if(eprime < 0.0) { eprime = 0.0; }

  return eprime;
  
}

int sheath0(void)
{
  int ii;
  for (ii = 0; ii < nc/2; ii++) {
    if (sp_n[1][ii] && sp_n[1][ii+1])
     if (sp_n[0][ii]/sp_n[1][ii] + sp_n[0][ii+1]/sp_n[1][ii+1] >= 1.2) break;
  }
  return ii;
}

int sheath1(void)
{
  int ii;
  for (ii = nc; ii > nc/2; ii--) {
    if (sp_n[1][ii] && sp_n[1][ii-1])
      if (sp_n[0][ii]/sp_n[1][ii] + sp_n[0][ii-1]/sp_n[1][ii-1] >= 1.2) break;
  }
  return nc-ii;
}

/* sheath diagnostic for which E field reverses sign */
SCALAR sheath_E0(int left, SCALAR *vc)
{
  int ii, inc, end, frac, start;
  SCALAR E0, Ej, Ek, phij, phik, sdist;
  if(left) {start=ii=0; inc = 1; end=nc; }
  else { start=ii=nc; inc = -1; end=0; }
  E0 = e[ii];
  if(!E0) { return 0.; }
  for(;ii != end; ii += inc)
    {
      if(e[ii+inc]*e[ii] <= 0)
	{
	  if(inc > 0)
	    {
	      Ej = e[ii]; Ek = e[ii+inc];	      
	      phij = phi[ii]; phik = phi[ii+inc];
	    }
	  else
	    {
	      Ej = e[ii+inc]; Ek = e[ii];
	      phij = phi[ii+inc]; phik = phi[ii];
	    }
	  frac = Ej/(Ej - Ek);
	  *vc = frac*phik + (1. - frac)*phij;
	  sdist = dx*(inc*frac + (SCALAR)ii); 
	  if(inc < 0) { sdist = length - sdist; }
	  return sdist;
	}
    }
  *vc = phi[end];
  return length;
}


/************************************************************************/
/* time history accumulator; calculates and stores all history values,  */
/* and performs combing on history values when low on memory            */

void history(void)
{
  SCALAR sle0, sre0, vctemp;
  register int i, j, isp, k;
  static SCALAR jtemp[NSMAX];
  static int count=1, nrftcount=0, nrfcount=0;
  
  /******   MCC rates  ******/
  if(it[0]) {
    for(i=0; i<ndiag; i++)
      for (j=0; j<ng; j++)
				rate_show[i][j] = (rate_show[i][j]*(it[0]-1) +mccrate[i][j])/((SCALAR)it[0]);
  }
  
  /****** sheath widths ****/
  s1 = s_1();
  s2 = s_2();  
  
  if(fabs(s1) > length) { s1 = 0.0; }
  if(fabs(s2) > length) { s2 = 0.0; }
  s1s2 = s1 + s2;

  if(nsp > 1)
    {
      s1n = dx*sheath0();
      s2n = dx*sheath1();
    }
  else
    {
      s1n = s2n = 0;
    }
  s1s2n = s1n + s2n;
  sle0 = sheath_E0(1, &vctemp);
  sre0 = sheath_E0(0, &vctemp);

  /*****************************************/
  /*  Fixing the diagnostic arrays         */

	if (vel_dist_accum) velocity();

  if(nsp && it[0]) {
    for (isp=0; isp<nsp; isp++) {
      jtemp[isp] += dt*jwall[isp];
      
      for (i=0; i<nbin_mid[isp]; i++) 
				fe_mid_show[isp][i] = fe_mid[isp][i]/(it[isp]*sqrt(e_mid_array[isp][i]+DBL_MIN));
      
      for (i=0; i< ng; i++) {
				if(sp_n_k[isp][i] <1e-10) {
					sp_u_x_show[isp][i] = sp_u_y_show[isp][i] = sp_u_z_show[isp][i] = 0.0;
					sp_ke_show[isp][i] = sp_ke_x_show[isp][i] = sp_ke_y_show[isp][i] = sp_ke_z_show[isp][i] = 0.0;
				}
				else {
					sp_u_x_show[isp][i] = sp_j_x[isp][i]/sp_n_k[isp][i];
					sp_u_y_show[isp][i] = sp_j_y[isp][i]/sp_n_k[isp][i];
					sp_u_z_show[isp][i] = sp_j_z[isp][i]/sp_n_k[isp][i];
					sp_ke_x_show[isp][i] = Escale[isp]*sp_ke_x[isp][i]/sp_n_k[isp][i];
					sp_ke_y_show[isp][i] = Escale[isp]*sp_ke_y[isp][i]/sp_n_k[isp][i];
					sp_ke_z_show[isp][i] = Escale[isp]*sp_ke_z[isp][i]/sp_n_k[isp][i];
					sp_ke_show[isp][i]=sp_ke_x_show[isp][i]+sp_ke_y_show[isp][i]+sp_ke_z_show[isp][i];
				}
				jdote_show[isp][i]= jdote[isp][i]*jdote_scale[isp]/it[isp];
				sp_j_x_show[isp][i] = sp_j_x[isp][i]*j_scale[isp];
				sp_j_y_show[isp][i] = sp_j_y[isp][i]*j_scale[isp];
				sp_j_z_show[isp][i] = sp_j_z[isp][i]*j_scale[isp];
      }
    }
  }
  
  /*****************************************/
  
  if(nfft) {
    if (thist_hi >= nfft) {
      freqanalysis();
      thist_hi=0;
    }	
    cur_hist[thist_hi]= exti;
    pow_hist[thist_hi]= exti*phi[0];
    phi_hist[0][thist_hi] = phi[0];
    phi_hist[1][thist_hi] = phi[nc/2];
    Local_t_array[thist_hi] = t;
    thist_hi++;
    
    /********************************************/
    /* Calculating the time ave. moments        */
	}
	if (n_ave){
		if(nrftcount == n_ave) {
			nrfcount++;
			
			for(isp=0; isp<nsp; isp++){
				for(i=0; i<ng; i++) {
					if (sp_n_ave[isp][i]>DEN_MIN){
						sp_u_x_ave_show[isp][i] = sp_j_x_ave[isp][i]/sp_n_ave[isp][i];
						sp_u_y_ave_show[isp][i] = sp_j_y_ave[isp][i]/sp_n_ave[isp][i];
						sp_u_z_ave_show[isp][i] = sp_j_z_ave[isp][i]/sp_n_ave[isp][i];

						sp_ke_x_ave_show[isp][i] = sp_ke_x_ave[isp][i]/sp_n_ave[isp][i];
						sp_ke_y_ave_show[isp][i] = sp_ke_y_ave[isp][i]/sp_n_ave[isp][i];
						sp_ke_z_ave_show[isp][i] = sp_ke_z_ave[isp][i]/sp_n_ave[isp][i];

						Tx_ave_show[isp][i] = sp_ke_x_ave_show[isp][i]-sqr(sp_u_x_ave_show[isp][i]);
						Ty_ave_show[isp][i] = sp_ke_y_ave_show[isp][i]-sqr(sp_u_y_ave_show[isp][i]);
						Tz_ave_show[isp][i] = sp_ke_z_ave_show[isp][i]-sqr(sp_u_z_ave_show[isp][i]);

						sp_ke_x_ave_show[isp][i] *= Escale[isp];
						sp_ke_y_ave_show[isp][i] *= Escale[isp];
						sp_ke_z_ave_show[isp][i] *= Escale[isp];

						Tx_ave_show[isp][i] *= Escale[isp];
						Ty_ave_show[isp][i] *= Escale[isp];
						Tz_ave_show[isp][i] *= Escale[isp];
					}
					else{
						sp_u_x_ave_show[isp][i] = 0.0;
						sp_u_y_ave_show[isp][i] = 0.0;
						sp_u_z_ave_show[isp][i] = 0.0;

						sp_ke_x_ave_show[isp][i] = 0.0;
						sp_ke_y_ave_show[isp][i] = 0.0;
						sp_ke_z_ave_show[isp][i] = 0.0;

						Tx_ave_show[isp][i] = 0.0;
						Ty_ave_show[isp][i] = 0.0;
						Tz_ave_show[isp][i] = 0.0;
					}

					T_ave_show[isp][i] = Tx_ave_show[isp][i]+Ty_ave_show[isp][i]+Tz_ave_show[isp][i]/3;
					sp_ke_ave_show[isp][i] = sp_ke_x_ave_show[isp][i]+sp_ke_y_ave_show[isp][i]+sp_ke_z_ave_show[isp][i];

					sp_ke_x_ave[isp][i] = sp_ke_x[isp][i];
					sp_ke_y_ave[isp][i] = sp_ke_y[isp][i];
					sp_ke_z_ave[isp][i] = sp_ke_z[isp][i];

					sp_n_ave_show[isp][i] = sp_n_ave[isp][i]/n_ave; sp_n_ave[isp][i] = sp_n[isp][i];
					sp_j_x_ave_show[isp][i] = sp_j_x_ave[isp][i]/n_ave; sp_j_x_ave[isp][i] = sp_j_x[isp][i];
					sp_j_y_ave_show[isp][i] = sp_j_y_ave[isp][i]/n_ave; sp_j_y_ave[isp][i] = sp_j_y[isp][i];
					sp_j_z_ave_show[isp][i] = sp_j_z_ave[isp][i]/n_ave; sp_j_z_ave[isp][i] = sp_j_z[isp][i];
				}
				for (i=0; i<nbin_mid[isp]; i++){
					sp_fe_show[isp][i] = sp_fe_ave[isp][i]/(n_ave*sqrt(e_mid_array[isp][i]+DBL_MIN)); 
					sp_fe_ave[isp][i] = sp_fe[isp][i];
				}
			}
			for (i=0; i<ng; i++) {
				e_ave_show[i]= e_ave[i]/n_ave;      e_ave[i]  = e[i];
				phi_ave_show[i]= phi_ave[i]/n_ave;  phi_ave[i]= phi[i];
			}
			
			nrftcount = 1;
		}
		else {
			for(isp=0; isp<nsp; isp++){
				for(i=0; i<ng; i++) {
					sp_n_ave[isp][i] += sp_n[isp][i];
					sp_j_x_ave[isp][i] += sp_j_x[isp][i];
					sp_j_y_ave[isp][i] += sp_j_y[isp][i];
					sp_j_z_ave[isp][i] += sp_j_z[isp][i];
					sp_ke_x_ave[isp][i]+= sp_ke_x[isp][i];
					sp_ke_y_ave[isp][i]+= sp_ke_y[isp][i];
					sp_ke_z_ave[isp][i]+= sp_ke_z[isp][i];
				}
				for (i=0; i<nbin_mid[isp]; i++)
					sp_fe_ave[isp][i] += sp_fe[isp][i];
			}
			for (i=0; i<ng; i++) {
				e_ave[i]   += e[i];
				phi_ave[i] += phi[i];
			}
			nrftcount++;
		}
	}
	
	for(isp=0; isp<nsp; isp++)
		for (i=0; i<nbin_mid[isp]; i++)
			sp_fe[isp][i]=0.0;
	
	
	/********************************************/
	
	if (--count) return;		/* only accum every interval steps */
	if (hist_hi >= HISTMAX) {	/* comb time histories */
		for (isp=0; isp<nsp; isp++) {
			for (i=1, k=4; i<HISTMAX/4; i++, k+=4) {
        			  vel_hist[isp][i] = vel_hist[isp][k];
				np_hist[isp][i] = np_hist[isp][k];
				np_trapped[isp][i] = np_trapped[isp][k];
				np_untrapped[isp][i] = np_untrapped[isp][k];
				TE_trapped[isp][i] = TE_trapped[isp][k];
				TE_untrapped[isp][i] = TE_untrapped[isp][k];
				TE_particle[isp][i] = TE_particle[isp][k];
				kes_hist[isp][i] = kes_hist[isp][k];
				kes_x_hist[isp][i] = kes_x_hist[isp][k];
				kes_y_hist[isp][i] = kes_y_hist[isp][k];
				kes_z_hist[isp][i] = kes_z_hist[isp][k];
				jwall_hist[isp][i]= jwall_hist[isp][k];
			}
		}
		for (i=1, k=4; i<HISTMAX/4; i++, k+=4) {
			com_phi_hist[0][i] = com_phi_hist[0][k];
			com_phi_hist[1][i] = com_phi_hist[1][k];
			com_pow_hist[i] = com_pow_hist[k];
			com_cur_hist[i] = com_cur_hist[k];
			wall_sigma_hist[i] = wall_sigma_hist[k];
			ese_hist[i] = ese_hist[k];
			s1_hist[i] = s1_hist[k];
                        s2_hist[i] = s2_hist[k];
                        s1s2_hist[i] = s1s2_hist[k];
			s1n_hist[i] = s1n_hist[k];
                        s2n_hist[i] = s2n_hist[k];
                        s1s2n_hist[i] = s1s2n_hist[k];
			s1e0[i] = s1e0[k];
			s2e0[i] = s2e0[k];
			v2e0[i] = v2e0[k];
			t_array[i] = t_array[k];
		}
		hist_hi = i;
		interval *= 4;	
	}
	for (isp=0; isp<nsp; isp++) {		/* accumulate histories */
	  /****** bulk velocity -- for mobility calculation *****/
	  for(i=0, vel_hist[isp][hist_hi]=0.0; i<np[isp];i++)
	    { 
	      vel_hist[isp][hist_hi] += vx[isp][i]; 
	    }
	  vel_hist[isp][hist_hi] *= vscale;
	  if(np[isp] > 0) { vel_hist[isp][hist_hi] /= np[isp]; }

		jwall_hist[isp][hist_hi]= jtemp[isp];
		jtemp[isp]= 0.0;
		if(hist_hi)  jwall_hist[isp][hist_hi] +=jwall_hist[isp][hist_hi-1]; 
		np_hist[isp][hist_hi] = np[isp];
		np_trapped[isp][hist_hi] = N_trapped[isp];
		np_untrapped[isp][hist_hi] = N_untrapped[isp];
		TE_trapped[isp][hist_hi] = E_trapped[isp];
		TE_untrapped[isp][hist_hi] = E_untrapped[isp];
		TE_particle[isp][hist_hi] = E_particles[isp];
		kes_hist[isp][hist_hi] = 0;
		kes_x_hist[isp][hist_hi] = 0;
		kes_y_hist[isp][hist_hi] = 0;
		kes_z_hist[isp][hist_hi] = 0;

		for (j=0; j<ng; j++){
			kes_hist[isp][hist_hi]+=sp_ke_x[isp][j]+sp_ke_y[isp][j]+sp_ke_z[isp][j];
			kes_x_hist[isp][hist_hi]+=sp_ke_x[isp][j];
			kes_y_hist[isp][hist_hi]+=sp_ke_y[isp][j];
			kes_z_hist[isp][hist_hi]+=sp_ke_z[isp][j];
		}
		
		if(np[isp]) {
			kes_hist[isp][hist_hi] *= Escale[isp]/np[isp];  // note: 1/4 scaling in initwin.c
			kes_x_hist[isp][hist_hi] *= Escale[isp]/np[isp];
			kes_y_hist[isp][hist_hi] *= Escale[isp]/np[isp];
			kes_z_hist[isp][hist_hi] *= Escale[isp]/np[isp];
		}
	}
	t_array[hist_hi] = t;
	wall_sigma_hist[hist_hi] = sigma;  
	com_cur_hist[hist_hi] = exti;
	com_phi_hist[0][hist_hi] = phi[0];
	com_phi_hist[1][hist_hi] = phi[nc/2];
	com_pow_hist[hist_hi] = exti*phi[0];
	
	s1_hist[hist_hi] = s1;
        s2_hist[hist_hi] = s2;
        s1s2_hist[hist_hi] = s1s2;
	s1n_hist[hist_hi] = s1n;
        s2n_hist[hist_hi] = s2n;
        s1s2n_hist[hist_hi] = s1s2n;
	s1e0[hist_hi] = sle0;
	s2e0[hist_hi] = sre0;
	v2e0[hist_hi] = vctemp;

	ese_hist[hist_hi] = 0.5*(rho[0]*phi[0] + rho[nc]*phi[nc]);
	for (i=1; i<nc; i++) ese_hist[hist_hi] += rho[i]*phi[i];
	ese_hist[hist_hi] *= 0.5*area*dx;
	ese_hist[hist_hi] -= epsilon*e[nc]*area*phi[nc]; // what is this here for, phi[nc]==0
	ese_hist[hist_hi] += .5*sigma*area*phi[0];  // add in induced wall charge.

	ese_hist[hist_hi] = fabs(ese_hist[hist_hi]+DBL_MIN);
	for (isp=0; isp<nsp; isp++)
		kes_hist[isp][hist_hi] = kes_hist[isp][hist_hi]+DBL_MIN;
	
	hist_hi++; 
	count = interval;
}

/***************************************************************/

void freqanalysis()
{
	register int i, j;
	static int init_freq_flag=1;
	static SCALAR *temp1, *temp2;
	
	if(init_freq_flag) {
		if(!(temp1= (SCALAR *)malloc((nfft/2)*sizeof(SCALAR)))
			 || !(temp2= (SCALAR *)malloc((nfft/2)*sizeof(SCALAR))))
			puts("Null ptr in freqanalysis()");
		init_freq_flag=0;
	}
	
	/******************************************************/
	for(i=0; i< nfft; i++) {
		cur_fft[i]= cur_hist[i];
		pow_fft[i]= pow_hist[i];
		phi_fft[i]= phi_hist[0][i];
		mphi_fft[i]= phi_hist[1][i];
	} 
	
	realft(phi_fft-1, freq_hi, 1);
	realft(mphi_fft-1,freq_hi, 1);
	realft(cur_fft-1, freq_hi, 1);
	realft(pow_fft-1, freq_hi, 1);
	
	/******************************************************/
	/**** Computing mag and phase of the current signal ***/

	temp1[0]= fabs(cur_fft[0])/freq_hi/2;
	temp2[0]= (cur_fft[0] > 0.0) ? 0.0 : 180.0;
	for(i=1, j=2; i< freq_hi; i++, j+=2) {
		temp1[i]= sqrt(cur_fft[j]*cur_fft[j]+ cur_fft[j+1]*cur_fft[j+1])/freq_hi;
		if(fabs(cur_fft[j+1]) < 1e-30 && fabs(cur_fft[j]) < 1e-30)
			temp2[i]= 0.0;
		else
			temp2[i]= (180./M_PI)*atan2(cur_fft[j], cur_fft[j+1]);
	}
	for(i=0; i< freq_hi; i++) {
		cur_fft[i]= temp1[i];
		cur_fft[freq_hi+i]= temp2[i];
	}
	
	/******************************************************/ 
	/**** Computing mag and phase of the power signal *****/

	temp1[0]= fabs(pow_fft[0])/freq_hi/2;
	temp2[0]= (pow_fft[0] > 0.0) ? 0.0 : 180.0;
	for(i=1, j=2; i< freq_hi; i++, j+=2) {
		temp1[i]= sqrt(pow_fft[j]*pow_fft[j]+ pow_fft[j+1]*pow_fft[j+1])/freq_hi;
		if(fabs(pow_fft[j+1]) < 1e-30 && fabs(pow_fft[j]) < 1e-30)
			temp2[i]= 0.0;
		else
			temp2[i]= (180./M_PI)*atan2(pow_fft[j], pow_fft[j+1]);
	}
	for(i=0; i< freq_hi; i++) {
		pow_fft[i]= temp1[i];
		pow_fft[freq_hi+i]= temp2[i];
	}

	/******************************************************/  
	/**** Computing mag and phase of the voltage signal ***/

	temp1[0]= fabs(phi_fft[0])/freq_hi/2;
	temp2[0]= (phi_fft[0] > 0.0) ? 0.0 : 180.0;
	for(i=1, j=2; i< freq_hi; i++, j+=2) {
		temp1[i]= sqrt(phi_fft[j]*phi_fft[j]+ phi_fft[j+1]*phi_fft[j+1])/freq_hi;
		if(fabs(phi_fft[j+1]) < 1e-30 && fabs(phi_fft[j]) < 1e-30)
			temp2[i]= 0.0;
		else
			temp2[i]= (180./M_PI)*atan2(phi_fft[j], phi_fft[j+1]);
	}
	for(i=0; i< freq_hi; i++) {
		phi_fft[i]= temp1[i];
		phi_fft[freq_hi+i]= temp2[i];
	}

	/************************************************************/
	/**** Computing mag and phase of the mid-potential signal ***/

	temp1[0]= fabs(mphi_fft[0])/freq_hi/2;
	temp2[0]= (mphi_fft[0] > 0.0) ? 0.0 : 180.0;
	for(i=1, j=2; i< freq_hi; i++, j+=2) {
		temp1[i]= sqrt(mphi_fft[j]*mphi_fft[j]+ mphi_fft[j+1]*mphi_fft[j+1])/freq_hi;
		if(fabs(mphi_fft[j+1]) < 1e-30 && fabs(mphi_fft[j]) < 1e-30)
			temp2[i]= 0.0;
		else
			temp2[i]= (180./M_PI)*atan2(mphi_fft[j], mphi_fft[j+1]);
	}
	for(i=0; i< freq_hi; i++) {
		mphi_fft[i]= temp1[i];
		mphi_fft[freq_hi+i]= temp2[i];
	}
}

/***************************************************************/

void velocity(void)
{
	int i, isp, index, dum1;
	static int initflag=1;
	static SCALAR dvx[NSMAX], xoffset[NSMAX],dvy[NSMAX], yoffset[NSMAX],
	dvz[NSMAX], zoffset[NSMAX], normal;
	
	if(initflag) {
		normal = 1;
		initflag =0;
		for(isp=0; isp<nsp; isp++) {
			dvx[isp] = (vxu[isp] - vxl[isp])/((nvxbin[isp]-1)*vscale);
			xoffset[isp]= vxl[isp]/vscale/dvx[isp];
			dvy[isp] = (vyu[isp] - vyl[isp])/((nvybin[isp]-1)*vscale);
			yoffset[isp]= vyl[isp]/vscale/dvy[isp];
			dvz[isp] = (vzu[isp] - vzl[isp])/((nvzbin[isp]-1)*vscale);
			zoffset[isp]= vzl[isp]/vscale/dvz[isp];
		}
		for(isp=0; isp<nsp;isp++)
			{
				for(i=0; i<nvxbin[isp]; i++)
					vx_array[isp][i] = vxl[isp]+i*dvx[isp]*vscale;
				for(i=0; i<nvybin[isp]; i++)
					vy_array[isp][i] = vyl[isp]+i*dvy[isp]*vscale;
				for(i=0; i<nvzbin[isp]; i++)
					vz_array[isp][i] = vzl[isp]+i*dvz[isp]*vscale;			
			}
	}

	for(isp=0;isp<nsp;isp++){ 
		if (nvxbin[isp])
			for(i=0; i<np[isp];i++){ 
				index = (int)(vx[isp][i]/dvx[isp] - xoffset[isp] +0.5);
				if((index>=0) && (index < nvxbin[isp])){
					dum1 = (int) (x[isp][i] + 0.499);
					vx_dist[isp][dum1][index] += normal;
				}
			}
		if (nvybin[isp])
			for(i=0; i<np[isp];i++){ 
				index = (int)(vy[isp][i]/dvy[isp] - yoffset[isp] +0.5);
				if((index>=0) && (index < nvybin[isp])){
					dum1 = (int) (x[isp][i] + 0.499);
					vy_dist[isp][dum1][index] += normal;
				}
			}
		if (nvzbin[isp])
			for(i=0; i<np[isp];i++){ 
				index = (int)(vz[isp][i]/dvz[isp] - zoffset[isp] +0.5);
				if((index>=0) && (index < nvzbin[isp])){
					dum1 = (int) (x[isp][i] + 0.499);
					vz_dist[isp][dum1][index] += normal;
				}
				/* note temporary hokey normalization*/
			}
	}
}


