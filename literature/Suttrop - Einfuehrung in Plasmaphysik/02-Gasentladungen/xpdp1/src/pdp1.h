#include <math.h>
#include <stdio.h>
//#include <unistd.h>
#include <stdlib.h>
#include "xgscalar.h"
#define EPS0		8.8542e-12			/* (F/m)  */
#define NperTORR	8.3221e20
#define NSMAX 		7
#define HISTMAX 	32768 				/* upper bound on histories */

#ifndef max
#define max(x, y)       (((x) > (y)) ? (x) : (y))
#endif

#ifndef min
#define min(x, y)       (((x) < (y)) ? (x) : (y))
#endif

#ifndef DBL_MIN
#define DBL_MIN 	1E-200
#endif

#ifndef FLOAT_MIN
#define FLOAT_MIN 	1e-4
#endif

#ifndef DEN_MIN
#define DEN_MIN 	1
#endif

#ifndef True
#define True            1
#endif

#ifndef False
#define False           0
#endif

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

#ifndef onesixth
#define onesixth 0.1666666667
#endif

#ifndef one24
#define one24 0.041666667
#endif

#ifndef one12
#define one12 0.083333333
#endif

#ifndef onethird
#define onethird 0.333333333
#endif

#ifndef twothirds
#define twothirds 0.666666667
#endif

#ifdef linux
#define cosd(d) (cos(M_PI*d/180))
#define sind(d) (sin(M_PI*d/180))
#endif

#define BIT_REVERSE     1
#define RANDOM          0

#define sqr(a)                           ((a)*(a)) 

/*  PG, attempt to read the gas type from input deck */
#define HELIUM 1
#define ARGON 2
#define NEON 3
#define OXYGEN 4
#define MCC 5


/*******************************************************************/

SCALAR nc2p, length, area, rhoback, backj, dde, epsilon, b, psi, extr,
             extl, extc, w0, dcbias, acbias, extq, extq_1, extq_2, extq_3,
             exti, sigma, oldsigma, dx, dt, vxscale, vscale, xnc, pressure,
             gtemp, ramp, theta0, risetime, seec[NSMAX], jwall[NSMAX], jnorm[NSMAX],
             q[NSMAX], m[NSMAX], qm[NSMAX], jj0[NSMAX][2], v0[NSMAX][2],
             vt[NSMAX][2], vc[NSMAX][2], v0y[NSMAX], vty[NSMAX], v0z[NSMAX], 
						 vtz[NSMAX],
             tx[NSMAX], tz[NSMAX], sx[NSMAX], sz[NSMAX], emin[NSMAX], de[NSMAX],
             dtheta[NSMAX], enter[NSMAX][2], emin_mid[NSMAX], de_mid[NSMAX],
             xs_mid[NSMAX], xf_mid[NSMAX], a_scale[NSMAX], sp_chi_scale[NSMAX],
             Escale[NSMAX], jdote_scale[NSMAX], j_scale[NSMAX];

int   nsp, nc, ng, secondary, ionspecies, ecollisional, icollisional,
             hist_hi, thist_hi, freq_hi, interval, nsmoothing, ntimestep,
             nfft, n_ave, dcramped, reflux, np[NSMAX], nbin[NSMAX],
             inject[NSMAX], nbin_mid[NSMAX], sp_k[NSMAX], it[NSMAX], maxnp[NSMAX],
             k_count[NSMAX], ndiag, gas, psource, nstrt, vel_dist_accum, 
						 vxloader[NSMAX][2], vyloader[NSMAX], vzloader[NSMAX], N_trapped[NSMAX],
						 N_untrapped[NSMAX], E_trapped[NSMAX], E_untrapped[NSMAX], E_particles[NSMAX];

long int seed;

double t, ttemp;
SCALAR s1, s2, s1s2, s1n, s2n, s1s2n;

char src, **rate_title;

int wraparound;

SCALAR **x, **vx, **vy, **vz, **sp_n, *rho, *e, *phi, *a, *x_grid,
              **fe, **ftheta, **e_array, **th_array, **fe_mid, **fe_mid_show,
						 **sp_fe_ave, **sp_fe_show, **sp_fe,
             **e_mid_array, *e_ave, *e_ave_show, *phi_ave, *phi_ave_show,
             **sp_n_ave, **sp_n_ave_show, 
             **sp_n_0, **sp_n_k, **sp_n_mcc, *chi,
             **jdote, **jdote_show,
             **mccrate, **rate_show, **np_trapped, **np_untrapped, **kes_x_hist,
						 **TE_trapped, **TE_untrapped, **TE_particle;

/* stuff for history diagnostics */

SCALAR *t_array, **np_hist, **jwall_hist, **phi_hist, **vel_hist,
      *wall_sigma_hist, *ese_hist, 
			**kes_hist, **kes_x_hist,  **kes_y_hist,  **kes_z_hist,
			*cur_hist,
      *com_cur_hist, *pow_hist, *com_pow_hist, **com_phi_hist, *f_array,
  *s1_hist, *s2_hist, *s1s2_hist,
  *s1n_hist, *s2n_hist, *s1s2n_hist, 
      *cur_fft, *phi_fft, *pow_fft, *mphi_fft, *Local_t_array;
SCALAR *s1e0, *s2e0, *v2e0;

/* stuff for velocity moments diagnostics */

SCALAR **sp_ke_x, **sp_ke_y, **sp_ke_z,
      **sp_ke_x_show, **sp_ke_y_show, **sp_ke_z_show,
			**sp_ke_x_ave, **sp_ke_y_ave, **sp_ke_z_ave,
			**sp_ke_x_ave_show, **sp_ke_y_ave_show, **sp_ke_z_ave_show,
			**sp_ke_show, **sp_ke_ave_show,
			**sp_j_x, **sp_j_y, **sp_j_z,
			**sp_j_x_show, **sp_j_y_show, **sp_j_z_show,
			**sp_j_x_ave, **sp_j_y_ave, **sp_j_z_ave,
			**sp_j_x_ave_show, **sp_j_y_ave_show, **sp_j_z_ave_show,
			**sp_u_x_show, **sp_u_y_show, **sp_u_z_show,
			**sp_u_x_ave_show, **sp_u_y_ave_show, **sp_u_z_ave_show,
			**Tx_ave, **Ty_ave, **Tz_ave,
			**Tx_ave_show, **Ty_ave_show, **Tz_ave_show, **T_ave_show;
			
			

/*stuff for velocity distribution diagnostics*/
SCALAR ***vx_dist, **vx_array,***vy_dist, **vy_array,***vz_dist, **vz_array;
SCALAR  vxu[NSMAX], vxl[NSMAX], vyu[NSMAX], vyl[NSMAX], vzu[NSMAX], vzl[NSMAX];
int    nvxbin[NSMAX], nvybin[NSMAX], nvzbin[NSMAX];

/*stuff for the volume source */
SCALAR endpts[2], vol_source, ionization_energy;

/*stuff for injection*/
SCALAR eold[2], W[NSMAX], sin4W[NSMAX], sin22W[NSMAX], cos22W[NSMAX], cos_psi, sin_psi;

SCALAR frand(void), tstrt;
//SCALAR bit_rever_maxwellian(void), maxwellian(void), maxwellian_flux(void);
SCALAR maxwellian(int);
// maxwellian_flux(int);
SCALAR distribution(int,int,int), distribution_flux(int,int,int);
void maxwellv(SCALAR *, SCALAR *, SCALAR *, SCALAR);
double revers_base(int,int), base2(void), revers(unsigned int);
void history(void), gather(int), adjust(int);
int start(void);
void fields(void), setrho(void);
void imp_move(int isp, const int EnergyFlag), exp_move(int isp, const int EnergyFlag), (*moveptr)(int isp, const int EnergyFlag);
void (*mccptr)(int isp);
void mccdiag_init(void);
void heliummcc(int isp), argonmcc(int isp), neonmcc(int isp), oxygenmcc(int isp), mcc(int isp);
void newivel(SCALAR, SCALAR, SCALAR *, SCALAR *, SCALAR *);

