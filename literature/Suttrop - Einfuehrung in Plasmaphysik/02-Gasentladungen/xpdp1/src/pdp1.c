#include "pdp1.h"
#include "xgrafix.h"

void display_title(void);
void InitWindows(void);

void main(int argc, char *argv[])
{
  display_title();        /* Display XPDP1 title           */
  XGInit(argc, argv, &t); /* Initialize XGrafix stuff      */
  
  start();                /* Allocate arrays and initailze */
  if(theRunWithXFlag) {
    mccdiag_init();       /* Initialize MCC rate diagnostic*/
    InitWindows();        /* Initialize diagnostic windows */
  }
  setrho();               /* Set initial charge density    */
  fields();               /* Initialize field arrays       */
  if(theRunWithXFlag)
    history();            /* Initialize history arrays     */

  XGStart();              /* Start XGrafix main loop       */
}

/***********************************************************/
/*  The main physics loop                                  */

void XGMainLoop(void)
{
  register int j, isp;
  int DiagFlag;
  SCALAR frac;
  
  t += dt; 

  if (theRunWithXFlag&&n_ave)
    DiagFlag=1;
  else
    DiagFlag=0;


	/*no subcycling loop*/
  for(isp=0; isp< nsp; isp++) {
    it[isp]++;                      /* Advance time for the species isp         */
    (*moveptr)(isp, DiagFlag);    /* Advance position and velocity            */
    adjust(isp);                    /* Remove particles that cross boundaries   */
		(*mccptr) (isp);                /* Monte Carlo collisions for species isp   */
	}
	for(isp=0; isp< nsp; isp++) {
		gather(isp);                    /* Assign charge densities to the grid      */
		for (j=0; j<ng; j++)
      sp_n[isp][j]=sp_n_k[isp][j];
	}
	/* gather needs to be in a seperate loop because new particles might
		 be created in mcc and adjust that might not have been weighted. */
	
	
	/*subcycling loop*/
	/*subcycling has bugs when used with mcc and secondaries*/
	/*the mcc is mostly right except when smoothing is used */
	/*the subcycling method is first order when using the explicit push*/
	/*ie the subcycling is for the implicit push */
	
#ifdef OLD_CODE  
  for(isp=0; isp< nsp; isp++) {
    if(!(k_count[isp]%sp_k[isp])) {
      it[isp]++;          /* Advance time for the species isp         */
      (*moveptr)(isp, DiagFlag);    /* Advance position and velocity            */
      adjust(isp);        /* Remove particles that cross boundaries   */
      (*mccptr) (isp);    /* Monte Carlo collisions for species isp   */
      gather(isp);        /* Assign charge densities to the grid      */
      k_count[isp]=0;     /* Reset species counter for subcycling     */
    }
    k_count[isp]++;
    frac = ((SCALAR)k_count[isp])/sp_k[isp];
    for (j=0; j<ng; j++)
      sp_n[isp][j]= (1- frac)*sp_n_0[isp][j]
				+frac*sp_n_k[isp][j] +sp_n_mcc[isp][j];
  }
#endif

  fields();
  
  if(theRunWithXFlag) history();
	
}

/***********************************************************/

void display_title(void)
{
  puts("\n\nXPDP1 Version 4.11");
  puts("Copyright (C) 1988-2002");
  puts("Regents of the University of California");
  puts("Plasma Theory and Simulation Group");
  puts("University of California - Berkeley\n");
}

/***************************************************************/
