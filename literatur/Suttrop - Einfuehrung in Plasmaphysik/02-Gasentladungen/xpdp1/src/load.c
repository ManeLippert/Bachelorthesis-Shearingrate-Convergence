#include "pdp1.h"

/***************************************************************/

void load(int initnp[][2])
{
  int i, ix, k, isp;
  SCALAR ddx;
  
  /************************************************/
  /* Load one species and one direction at a time */
  
  for (isp=0; isp<nsp; isp++) {
    np[isp]= 0;
    
    for (k=0; k<2; k++) {
      if (!initnp[isp][k]) continue;
      
      ddx= xnc/initnp[isp][k];
      if (np[isp] + initnp[isp][k] > maxnp[isp]) {
	puts("LOAD: too many particles, species ");
	exit(1);
      }
      if (vt[isp][k]==0. ) {    /* Loader for a COLD beam   */
	for(i=0; i<initnp[isp][k]; i++) {
	  ix = i +np[isp];
	  if(!i)  x[isp][ix]= .5*ddx;
	  else    x[isp][ix] = x[isp][ix-1] +ddx;
	  vx[isp][ix] = v0[isp][k];
	  vy[isp][ix] = v0y[isp]+vty[isp]*maxwellian(vyloader[isp]);
	  vz[isp][ix] = v0z[isp]+vtz[isp]*maxwellian(vzloader[isp]);
	}
	if(dde) {
	  for(i=0; i<initnp[isp][k]; i++) {
	    ix = np[isp] +i;
	    x[isp][ix] += dde*sin(2*M_PI*x[isp][ix]/xnc);
	  }
	}
      }
      else {  /* Loader for THERMAL distribution */
	for(i=0; i<initnp[isp][k]; i++)	{
	  ix = i +np[isp];
	  if(!i)  x[isp][ix]= .5*ddx;
	  else    x[isp][ix] = x[isp][ix-1] + ddx;
	  vx[isp][ix] = distribution(k,isp,vxloader[isp][k]);
	  vy[isp][ix] = v0y[isp]+vty[isp]*maxwellian(vyloader[isp]);
	  vz[isp][ix] = v0z[isp]+vtz[isp]*maxwellian(vzloader[isp]);
	}
	if(dde)	{
	  for(i=0; i<initnp[isp][k]; i++) {
	    ix = np[isp] +i;
	    x[isp][ix] += dde*sin(2*M_PI*x[isp][ix]/xnc);
	  }
	}
      }	/* end of THERMAL loader  "if( VT>0.." */
      np[isp] += initnp[isp][k];
      
    }	/* end  "for( k=0,1...)" loop */
  }	/* end for(isp=0 thru nsp-1..)" loop */
}	/* end LOAD  */

/***************************************************************/

