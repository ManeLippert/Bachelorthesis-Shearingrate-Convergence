#include "pdp1.h"

/***************************************************************/
void One_2_One(SCALAR *ary, int nmax);

void gather(int isp)
{
  register int i, j;
  register SCALAR s;
  
  for (j=0; j< ng; j++) {
    sp_n_0[isp][j]= sp_n_k[isp][j] +sp_n_mcc[isp][j];
    sp_n_mcc[isp][j]= sp_n_k[isp][j]= 0.0;
  }

  for (i=np[isp]-1; i>=0; i--) {
    j = x[isp][i];
    s = x[isp][i] - j;
    sp_n_k[isp][j]  += 1. - s;
    sp_n_k[isp][j+1]+= s;
  }
	sp_n_k[isp][0]  *= 2.;
	sp_n_k[isp][nc] *= 2.;

  /************************************************/
  /* Smoothing the charge density of each species */
  
  for(i=0; i< nsmoothing; i++) One_2_One(sp_n_k[isp], nc);
}

/***************************************************************/

void setrho(void)
{
  int j, isp;

  for(isp=0; isp<nsp; isp++) {
    k_count[isp] = 0;
    gather(isp);
    for (j=0; j<ng; j++) {
      sp_n[isp][j]= sp_n_k[isp][j];
      sp_n_mcc[isp][j]= 0.0;
    }
  }
}

/***************************************************************/
/* Smoothing the array using the 1-2-1 method with the proper   */
/* boundary conditions to conserve charge.                     */

void One_2_One(SCALAR *ary, int nmax)
{
  register int i;
  static int nlocal=0;
  static SCALAR *temp=0;
  
  if(nlocal < nmax) {
    if(temp) free(temp);
    temp = (SCALAR *)malloc((nmax+1)*sizeof(SCALAR));
    nlocal = nmax;
  }
  
	temp[0]= (ary[0] +ary[1])/2;
	for(i=1; i< nmax; i++)  temp[i]= (ary[i-1] +2*ary[i] +ary[i+1])/4;
	temp[nmax]= (ary[nmax] +ary[nmax-1])/2;

  for(i=0; i<= nmax; i++) ary[i]= temp[i];
}

/***************************************************************/
