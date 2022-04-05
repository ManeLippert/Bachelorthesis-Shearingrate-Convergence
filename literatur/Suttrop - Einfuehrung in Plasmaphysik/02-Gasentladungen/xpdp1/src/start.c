#include <unistd.h> 
#include "pdp1.h"
#include "xgrafix.h"

FILE *InputDeck;


/***********************************************************/
void species(int initnp[][2], int isp);
void Restore(char *filename);
void load(int [][2]);
void DiagArray(void);
void SpeciesDiagArray(void);

int start(void)
{

  SCALAR t, s;
  char aachar[200];
  int i, j, isp, imp_flag, initnp[NSMAX][2];
  
  /***********************************************/
  /* Open files and read in "general" parameters */
  
  if (!WasInputFileGiven) {
    puts("Syntax:");
    puts("\tPDP2 <inputdeck> [dumpfile]");
    exit(1);
  }
  InputDeck = fopen(theInputFile, "r");
  
  /*****************************************/
  /* Read lines until we get to numbers    */
  
  while (fscanf(InputDeck,"%d %d %lf %lf %lf %lf %lf %lf %lf",
		&nsp, &nc, &nc2p, &dt, &length, &area, &epsilon, &b, &psi) <9)
    fscanf(InputDeck, "%s", aachar);
  
  while (fscanf(InputDeck, "%lf %lf %lf %lf %lf %lf %lf",
		&rhoback, &backj, &dde, &extr, &extl, &extc, &extq) <7)
    fscanf(InputDeck, "%s", aachar);
  
  while (fscanf(InputDeck, "%d %c %lf %lf %lf %lf %lf", &dcramped,
		&src, &dcbias, &ramp, &acbias, &w0, &theta0) <7)
    fscanf(InputDeck, "%s", aachar);
  
  while (fscanf(InputDeck, "%d %d %d %d %d %d %d %d", &secondary, &ecollisional,
		&icollisional, &reflux, &nfft, &n_ave, &nsmoothing, &ntimestep) <8)
    fscanf(InputDeck, "%s", aachar);

  while (fscanf(InputDeck, "%lf %lf %d %lf %lf %d", &seec[0], &seec[1],
		&ionspecies, &pressure, &gtemp, &imp_flag) <6)
    fscanf(InputDeck, "%s", aachar);
  
  while (fscanf(InputDeck, "%d %d %d %lf %lf %lf %lf", &gas, &psource, &nstrt, &vol_source, &endpts[0], &endpts[1], &ionization_energy) <7)
    fscanf(InputDeck, "%s", aachar);

  wraparound = 0;
  if((src == 'W') || (src == 'w'))
    {
      wraparound = 1;
      src = 'v';
    }

  /******************************************************/
  /* Check for errors in the "general" input parameters */
  
  if (extl < 0. || extr < 0. || extc < 0.) {
    puts("START: dont be cute with R,L,C < 0\n");
    exit(1);
  }
  if (dt<= 0.0) {
    puts("START: dt <= 0!!\n");
    exit(1);
  }
  if (nsp > NSMAX) {
    puts("START: nsp > NSMAX.\n");
    exit(1);
  }
  
  /* make sure that oxygen has 3 species */

  if ( (gas == OXYGEN) && (nsp != 3) ){
    puts("START: OXYGEN, nsp must be 3\n");
    exit(1);
  }


  /* 
     by default the gas is argon argonmcc eturns without any action if
     collisions are off
     */
  mccptr =argonmcc;

  if (ecollisional || icollisional){

    if ( (gas != OXYGEN) && (nsp == 3) ){
      puts("START: if you want nsp = 3 AND collisions on, gas must be OXYGEN\n");
      exit(1);
    }

    switch (gas){

    case HELIUM:
      mccptr =heliummcc;
      break;

    case ARGON:
      mccptr =argonmcc;
      break;

    case NEON:
      mccptr =neonmcc;
      break;

    case OXYGEN:
      mccptr =oxygenmcc;
      break;

    case MCC:
      mccptr =mcc;
      break;
    }
  }
	

  /* Check to see if nfft is an integer power of 2 */
  if(nfft) {
    for(i=0; i<= 31; i++) {
      if(nfft == (1<<i)) {
	i=0;
	break;
      }
    }
    if(i) {
      puts("START: nfft is not an integer power of 2");
      exit(1);
    }
  }

  if(imp_flag) {
    if(b>0.0) {
      puts("SORRY: The implicit version of PDP1 does not have a B field!");
      exit(1);
    }
    else moveptr = imp_move;
  }
  else moveptr = exp_move;
  
  xnc= nc;
  ng = nc+1;
  dx = length/xnc; 
  
  /***************************************/
  /* Allocate space for field parameters */
  
  x_grid= (SCALAR *)malloc(ng*sizeof(SCALAR));
  for(j=0; j<ng; j++) x_grid[j]= (SCALAR)j;
  
  sp_n    = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
  sp_n_0  = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
  sp_n_k  = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
  sp_n_mcc= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
  
  for(isp=0; isp<nsp; isp++) {
    sp_n[isp]    = (SCALAR *)calloc(ng,sizeof(SCALAR));
    sp_n_0[isp]  = (SCALAR *)calloc(ng,sizeof(SCALAR));
    sp_n_k[isp]  = (SCALAR *)calloc(ng,sizeof(SCALAR));
    sp_n_mcc[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
  }
  a  = (SCALAR *)calloc(ng,sizeof(SCALAR));
  rho= (SCALAR *)calloc(ng,sizeof(SCALAR));
  phi= (SCALAR *)calloc(ng,sizeof(SCALAR));
  e  = (SCALAR *)calloc(ng,sizeof(SCALAR));
  chi= (SCALAR *)calloc(ng,sizeof(SCALAR));

  if (theRunWithXFlag&&n_ave) /**** Allocate Dianostic arrays  *******/
    DiagArray();
    

  /**************************************/

  vscale = dx/dt;
  theta0 *= 2*M_PI/360.0;
  w0 *= 2*M_PI;
  epsilon *= EPS0;
  extq_1 = extq_2 = extq_3 = extq;
  interval = 1;

  seed = getpid(); /* This is the seed used in the random number generator */

  if(dcramped) {
    if(ramp*dcbias < 0.0) {
      fprintf(stderr, "\nWARNING: Ramp and dcbias have opposite signs.");
      fprintf(stderr, "\nThe external source will be discontineous!\n\n");
    }
    if(ramp) risetime= fabs(dcbias)/fabs(ramp);
    else     risetime= M_PI/w0;
  }

  /********************************************/
  /* Set up the parameters for each "species" */

  for (isp=0; isp<nsp; ++isp) {
    species(initnp, isp);

    vxscale= max(vxscale, v0[isp][0]+vt[isp][0]);
    vxscale= max(vxscale, -v0[isp][1]-vt[isp][1]);

    /* Scaling factors to convert from physical to code units and vis versa */
    qm[isp] = q[isp]/m[isp];
    q[isp] *= nc2p/area;
    Escale[isp] = .5*vscale*vscale*m[isp]/(1.602e-19);
    jdote_scale[isp]= q[isp]/dt;
    j_scale[isp] = q[isp]/dt;
    //		u_scale[isp] = vscale;
    if (sp_k[isp])
      jnorm[isp] = q[isp]/sp_k[isp]/dt;
    a_scale[isp] = qm[isp]*sp_k[isp]*dt*dt/dx;
    if(b>0.0 || imp_flag) a_scale[isp] *= 0.5;
		
    if(imp_flag)
      sp_chi_scale[isp] = 0.25*q[isp]*a_scale[isp]*sp_k[isp]/epsilon;
    else {
      sp_chi_scale[isp] = 0.0;
      jdote_scale[isp] *= 0.5;
      j_scale[isp] *= 0.5;
      //			u_scale[isp] *= 0.5;
    }
    /* Normalize velocities and energies */
    v0[isp][0] /= vscale;
    v0[isp][1] /= vscale;
    vc[isp][0] /= vscale;
    vc[isp][1] /= vscale;
    vt[isp][0] /= vscale;
    vt[isp][1] /= vscale;
    v0y[isp]/= vscale;
    vty[isp]/= vscale;
    v0z[isp]/= vscale;
    vtz[isp]/= vscale;
    emin[isp]  /= Escale[isp];
    de[isp]    /= Escale[isp];
    xs_mid[isp]/= dx;
    xf_mid[isp]/= dx;
    emin_mid[isp]/= Escale[isp];
    de_mid[isp]  /= Escale[isp];

    /* Normalize "enter" to = no. of particles injected PER dt */

    enter[isp][0] = jj0[isp][0]*dt*sp_k[isp]/(fabs(q[isp])+DBL_MIN);
    enter[isp][1] = jj0[isp][1]*dt*sp_k[isp]/(fabs(q[isp])+DBL_MIN);
    if(enter[isp][0]> 0.0 || enter[isp][1]> 0.0) inject[isp] = 1;
  }

  fclose(InputDeck);

  if(!vxscale) vxscale = 1e4;
  vxscale *= 5;

  /*
    psource =1, starts a uniform ionization source which compensates every
    ion leaving with an ion-electron pair
    */
  /*
    vol_source is a uniform ionization source the units are m-3 t-1
    */

  vol_source *= dt*area*(endpts[1]-endpts[0])/nc2p;//number of ionization events in a single dt																					
						       endpts[0] /= dx;
  endpts[1] /= dx;
  /* have to call here since at this point, svel has all it needs */

  if (psource||vol_source) 
    {
      //			svel();
      tstrt = 10*dt; /* default, start source after 10 dT */
      if (nstrt) tstrt = nstrt*dt;
    }


  /******************************/
  /** Allocate particle arrays **/

  x  = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
  vx = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
  vy = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
  vz = (SCALAR **)malloc(nsp*sizeof(SCALAR *));

  for(isp=0; isp<nsp; isp++) {
    x[isp]  = (SCALAR *)calloc(maxnp[isp],sizeof(SCALAR));
    vx[isp] = (SCALAR *)calloc(maxnp[isp],sizeof(SCALAR));
    vy[isp] = (SCALAR *)calloc(maxnp[isp],sizeof(SCALAR));
    vz[isp] = (SCALAR *)calloc(maxnp[isp],sizeof(SCALAR));
  }

  if (nsp && (!x[nsp-1] || !vx[nsp-1] || !vy[nsp-1] || !vz[nsp-1])) {
    puts("START: x, v malloc failed");
    exit(-1);
  }

  if (theRunWithXFlag)
    SpeciesDiagArray();

  /*****************************************************************/
  /* Calculating the coeff. for the mover in presence of a B field */

  cos_psi=cosd(psi);
  sin_psi=sind(psi);

  for (isp=0; isp<nsp; ++isp) {
    t= tan(.5*b*qm[isp]*dt*sp_k[isp]);
    s= 2*t/(1+t*t);

    /*this is for the injection mover*/
    W[isp]=0.25*b*qm[isp]*dt*sp_k[isp];
    sin4W[isp]=sin(4*W[isp]);
    sin22W[isp]=sqr(sin(2*W[isp]));
    cos22W[isp]=sqr(cos(2*W[isp]));
		
    tx[isp]= t*cos_psi;
    tz[isp]= t*sin_psi;
    sx[isp]= s*cos_psi;
    sz[isp]= s*sin_psi;
  }

  /**********************************************************/
  /* Set up array of bit-reversed indices. "revers(num)"    */
  /*  reverses the order of the base-i representation of    */
  /* the number "num", yielding a fraction between 0 and 1. */
  /* Index[] is also used in padjus.c.                      */

  /* for (i=0; i<1024; i++) Index[i] = 1024*revers(i); */
	/* not used anymore */

  /*********************************************************/
  /* Load in all species' particles, properly distributed. */
  /* Start from a dump file if given otherwise load the    */
  /* system with the initial distribution                  */

  if (WasDumpFileGiven)
    Restore(theDumpFile);
  else
    load(initnp);

  return(1);

}

/***************************************************************/

void species(int initnp[][2], int isp)
{
  char aachar[512];
  SCALAR iinitn, j0l, j0r, v0l, v0r, vcl, vcr, vtl, vtr, temin, temax,
    temin_mid, temax_mid;
	
  /* Read in SPECIES parameters... */
  while (fscanf(InputDeck, "%lf %lf %lf %lf %lf %d", &q[isp], &m[isp], &j0l,
		&j0r, &iinitn, &sp_k[isp]) <5)
    fscanf(InputDeck, "%s", aachar);

  while (fscanf(InputDeck, "%lf %lf %lf %d", 
		&v0l, &vtl, &vcl, &vxloader[isp][0])<4)
    fscanf(InputDeck, "%s", aachar);

  while (fscanf(InputDeck, "%lf %lf %lf %d", 
		&v0r, &vtr, &vcr, &vxloader[isp][1])<4)
    fscanf(InputDeck, "%s", aachar);
	
  //	while (fscanf(InputDeck, "%lf %lf %lf %lf %lf %lf %d %d", &v0l, &v0r, &vtl,
		      //								&vtr, &vcl, &vcr, &vxloader[isp][0], &vxloader[isp][1]) <8)
    //		fscanf(InputDeck, "%s", aachar);

  while (fscanf(InputDeck, "%lf %lf %d %lf %lf %d", &v0y[isp], &vty[isp], &vyloader[isp], &v0z[isp],
		&vtz[isp], &vzloader[isp]) <6)
    fscanf(InputDeck, "%s", aachar);

  while (fscanf(InputDeck, "%d %lf %lf %d",
		&nbin[isp], &temin, &temax, &maxnp[isp]) <4)
    fscanf(InputDeck, "%s", aachar);
	
  while (fscanf(InputDeck, "%d %lf %lf %lf %lf", &nbin_mid[isp], &temin_mid,
		&temax_mid, &xs_mid[isp], &xf_mid[isp]) <5)
    fscanf(InputDeck, "%s", aachar);

  while (fscanf(InputDeck, "%lf %lf %d %lf %lf %d %lf %lf %d", &vxl[isp], &vxu[isp],
		&nvxbin[isp], &vyl[isp], &vyu[isp],
		&nvybin[isp], &vzl[isp], &vzu[isp],
		&nvzbin[isp]) <9)
    fscanf(InputDeck, "%s", aachar);

  if (nvxbin[isp]||nvybin[isp]||nvzbin[isp])
    vel_dist_accum = 1;
	
  /* Assign input parameters to arguments */
  jj0[isp][0] =  fabs(j0l);
  jj0[isp][1] =  fabs(j0r);
  v0[isp][0]  =  fabs(v0l);
  v0[isp][1]  = -fabs(v0r);
  vc[isp][0]  =  fabs(vcl);
  vc[isp][1]  = -fabs(vcr);
  vt[isp][0]  =  fabs(vtl);
  vt[isp][1]  =  -fabs(vtr);
  vty[isp]    = fabs(vty[isp]);
  vtz[isp]    = fabs(vtz[isp]);
  emin[isp]   =  temin;
  de[isp]     =  (temax-temin)/(nbin[isp]-1);
  dtheta[isp] =  0.5*M_PI/(nbin[isp]-1);
	
  emin_mid[isp]= temin_mid;
  de_mid[isp] = (temax_mid-temin_mid)/(nbin_mid[isp]);
	
  if((vt[isp][0] || v0[isp][0]) && (vt[isp][1] ==0.0 && v0[isp][1] ==0.0)) {
    initnp[isp][0] = iinitn*length*area/nc2p +.5;
    initnp[isp][1] = 0.0;
  }
  else if((vt[isp][1] || v0[isp][1]) && (vt[isp][0] ==0.0 && v0[isp][0] ==0.0)) {
    initnp[isp][0] = 0.0;
    initnp[isp][1] = iinitn*length*area/nc2p +.5;
  }
  else if(!vt[isp][0] && !vt[isp][1] && !v0[isp][0] && !v0[isp][1]) {
    initnp[isp][0] = iinitn*length*area/nc2p +.5;
    initnp[isp][1] = 0.0;
  }
  else {
    initnp[isp][0] = .5*iinitn*length*area/nc2p +.5;
    initnp[isp][1] = .5*iinitn*length*area/nc2p +.5;
  }
}   /* End SPECIES */

/***************************************************************/

void DiagArray(void)
{
  /**** Allocate Diagnostic arrays  *******/
  /****** For the frgd version also allocate space for diagnostics *****/


  SCALAR df;
  int i, isp;
  
  if (n_ave){
    sp_j_x  = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_j_y  = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_j_z  = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_x = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_y = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_z = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    jdote = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_u_x_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_u_y_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_u_z_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_j_x_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_j_y_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_j_z_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_x_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_y_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_z_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_show= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    jdote_show= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    
    for(isp=0; isp<nsp; isp++) {
      sp_j_x[isp]  = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_j_y[isp]  = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_j_z[isp]  = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_x[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_y[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_z[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      jdote[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_u_x_show[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_u_y_show[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_u_z_show[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_j_x_show[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_j_y_show[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_j_z_show[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_show[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_x_show[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_y_show[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_z_show[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      jdote_show[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
    }
  
    /* Allocate space for time history diagnostics */
    t_array = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
    ese_hist= (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
    wall_sigma_hist= (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
    com_pow_hist   = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
    com_cur_hist   = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
    com_phi_hist   = (SCALAR **)malloc(2*sizeof(SCALAR *));
    com_phi_hist[0]= (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
    com_phi_hist[1]= (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
		
    np_hist= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    vel_hist= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    np_trapped = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    np_untrapped = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    TE_trapped = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    TE_untrapped = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    TE_particle = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    kes_hist= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    kes_x_hist= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    kes_y_hist= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    kes_z_hist= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    jwall_hist= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    for(i=0; i<nsp; i++) {
      np_hist[i]= (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
      vel_hist[i]= (SCALAR *)calloc(HISTMAX, sizeof(SCALAR));
      np_trapped[i] = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR *));
      np_untrapped[i] = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR *));
      TE_trapped[i] = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR *));
      TE_untrapped[i] = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR *));
      TE_particle[i] = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR *)); 
      kes_hist[i]= (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
      kes_x_hist[i]= (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
      kes_y_hist[i]= (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
      kes_z_hist[i]= (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
      jwall_hist[i]= (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
    }
  
    /* Allocate space for ave diagnostics */
    sp_n_ave = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_j_x_ave = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_j_y_ave = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_j_z_ave = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_x_ave= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_y_ave= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_z_ave= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    s1_hist = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
    s2_hist = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
    s1s2_hist = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
    s1n_hist = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
    s2n_hist = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
    s1s2n_hist = (SCALAR *)calloc(HISTMAX,sizeof(SCALAR));
    s1e0 = (SCALAR *)calloc(HISTMAX, sizeof(SCALAR));
    s2e0 = (SCALAR *)calloc(HISTMAX, sizeof(SCALAR));
    v2e0 = (SCALAR *)calloc(HISTMAX, sizeof(SCALAR));
    Tx_ave = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    Ty_ave = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    Tz_ave = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    Tx_ave_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    Ty_ave_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    Tz_ave_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    T_ave_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_n_ave_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_u_x_ave_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_u_y_ave_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_u_z_ave_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_j_x_ave_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_j_y_ave_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_j_z_ave_show = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_x_ave_show= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_y_ave_show= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_z_ave_show= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_ke_ave_show= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    for(isp=0; isp<nsp; isp++) {
      sp_n_ave[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_n_ave_show[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_j_x_ave[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_j_y_ave[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_j_z_ave[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_x_ave[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_y_ave[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_z_ave[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_u_x_ave_show[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_u_y_ave_show[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_u_z_ave_show[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_j_x_ave_show[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_j_y_ave_show[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_j_z_ave_show[isp] = (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_x_ave_show[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_y_ave_show[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_z_ave_show[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      sp_ke_ave_show[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      Tx_ave[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      Ty_ave[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      Tz_ave[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      Tx_ave_show[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      Ty_ave_show[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      Tz_ave_show[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
      T_ave_show[isp]= (SCALAR *)calloc(ng,sizeof(SCALAR));
    }
 
    phi_ave= (SCALAR *)calloc(ng,sizeof(SCALAR));
    phi_ave_show= (SCALAR *)calloc(ng,sizeof(SCALAR));
  	
    e_ave= (SCALAR *)calloc(ng,sizeof(SCALAR));
    e_ave_show= (SCALAR *)calloc(ng,sizeof(SCALAR));
  }

  if (nfft){
    freq_hi= nfft/2;
    df= 1/(nfft*dt);
    f_array = (SCALAR  *)malloc(nfft*sizeof(SCALAR));
    cur_hist= (SCALAR  *)malloc(nfft*sizeof(SCALAR));
    pow_hist= (SCALAR  *)malloc(nfft*sizeof(SCALAR));
    phi_hist= (SCALAR **)malloc(2*sizeof(SCALAR *));
    phi_hist[0]= (SCALAR *)malloc(nfft*sizeof(SCALAR));
    phi_hist[1]= (SCALAR *)malloc(nfft*sizeof(SCALAR));
			
    Local_t_array = (SCALAR *) malloc(nfft*sizeof(SCALAR));
			
    pow_fft= (SCALAR *)malloc(nfft*sizeof(SCALAR));
    cur_fft= (SCALAR *)malloc(nfft*sizeof(SCALAR));
    phi_fft= (SCALAR *)malloc(nfft*sizeof(SCALAR));
    mphi_fft=(SCALAR *)malloc(nfft*sizeof(SCALAR));
    for(i=0; i< nfft; i++) {
      f_array[i]= i*df;
      pow_fft[i] = cur_fft[i] = phi_fft[i] = mphi_fft[i] = 0.0;
    }
  }
}

void SpeciesDiagArray()
{	

  int i, isp;

  /* allocate space for v_dist */
  if (vel_dist_accum)
    {
      vx_dist = (SCALAR ***)malloc(nsp*sizeof(SCALAR **));
      vy_dist = (SCALAR ***)malloc(nsp*sizeof(SCALAR **));
      vz_dist = (SCALAR ***)malloc(nsp*sizeof(SCALAR **));
      for(isp=0; isp<nsp; isp++){
	  vx_dist[isp] = (SCALAR **)malloc(ng*sizeof(SCALAR*));
	  vy_dist[isp] = (SCALAR **)malloc(ng*sizeof(SCALAR*));
	  vz_dist[isp] = (SCALAR **)malloc(ng*sizeof(SCALAR*));
	  for (i=0; i<ng; i++){ 
	      vx_dist[isp][i] = (SCALAR *)calloc(nvxbin[isp],sizeof(SCALAR));
	      vy_dist[isp][i] = (SCALAR *)calloc(nvybin[isp],sizeof(SCALAR));
	      vz_dist[isp][i] = (SCALAR *)calloc(nvzbin[isp],sizeof(SCALAR));
	    } 
	}

      vx_array = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
      vy_array = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
      vz_array = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
      for(isp=0; isp<nsp; isp++){
	  vx_array[isp] = (SCALAR *)malloc(nvxbin[isp]*sizeof(SCALAR));
	  vy_array[isp] = (SCALAR *)malloc(nvybin[isp]*sizeof(SCALAR));
	  vz_array[isp] = (SCALAR *)malloc(nvzbin[isp]*sizeof(SCALAR));
	}
    }

  /*************************************************/
  /* Allocate space for the dist. functions arrays */

  if(theRunWithXFlag) {
    fe      = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    e_array = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    ftheta  = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    th_array= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    fe_mid     = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    fe_mid_show= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_fe    = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_fe_ave    = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    sp_fe_show    = (SCALAR **)malloc(nsp*sizeof(SCALAR *));
    e_mid_array= (SCALAR **)malloc(nsp*sizeof(SCALAR *));
		
    for(isp=0; isp<nsp; isp++) {
      fe[isp]      = (SCALAR *)malloc(nbin[isp]*sizeof(SCALAR));
      ftheta[isp]  = (SCALAR *)malloc(nbin[isp]*sizeof(SCALAR));
      e_array[isp] = (SCALAR *)malloc(nbin[isp]*sizeof(SCALAR));
      th_array[isp]= (SCALAR *)malloc(nbin[isp]*sizeof(SCALAR));
      fe_mid[isp]  = (SCALAR *)malloc(nbin_mid[isp]*sizeof(SCALAR));
      fe_mid_show[isp]= (SCALAR *)malloc(nbin_mid[isp]*sizeof(SCALAR));
      sp_fe[isp]  = (SCALAR *)malloc(nbin_mid[isp]*sizeof(SCALAR));
      sp_fe_ave[isp]  = (SCALAR *)malloc(nbin_mid[isp]*sizeof(SCALAR));
      sp_fe_show[isp]  = (SCALAR *)malloc(nbin_mid[isp]*sizeof(SCALAR));
      e_mid_array[isp]= (SCALAR *)malloc(nbin_mid[isp]*sizeof(SCALAR));
    }
		
    for(isp=0; isp<nsp; isp++) {
      for(i=0; i<nbin[isp]; i++) {
	fe[isp][i]      = 0.0;
	ftheta[isp][i]  = 0.0;
	e_array[isp][i] = emin[isp] +i*de[isp];
	th_array[isp][i]= i*dtheta[isp];
      }
      for(i=0; i<nbin_mid[isp]; i++) {
	fe_mid[isp][i] = fe_mid_show[isp][i] = 0.0;
	sp_fe[isp][i] = sp_fe_ave[isp][i] = sp_fe_show[isp][i] = 0.0;
	e_mid_array[isp][i]= emin_mid[isp] +(i+.5)*de_mid[isp];
      }
    }
  }
}
