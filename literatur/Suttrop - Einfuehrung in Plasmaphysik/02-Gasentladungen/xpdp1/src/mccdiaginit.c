#include "pdp1.h"

void mccdiag_init(void)
{
  int i, j;
  
  if(ecollisional==0 && icollisional==0) {
    ndiag = 0;
    return;
  }
  
	/* default gas = argon */

	if ( (gas < HELIUM) || (gas > OXYGEN) ) gas = ARGON;

	switch (gas) {

	case HELIUM:
		ndiag = 5;
		
		rate_title = (char **)malloc(ndiag*sizeof(char *));
		rate_title[0]= "MCC Elastic";
		rate_title[1]= "MCC Excitation";
		rate_title[2]= "MCC Ionization";
		rate_title[3]= "MCC Ion Scattering";
		rate_title[4]= "MCC Charge X";
		break;
		
	case ARGON:
		ndiag = 5;
		
		rate_title = (char **)malloc(ndiag*sizeof(char *));
		rate_title[0]= "MCC Elastic";
		rate_title[1]= "MCC Excitation";
		rate_title[2]= "MCC Ionization";
		rate_title[3]= "MCC Ion Scattering";
		rate_title[4]= "MCC Charge X";
		break;

	case NEON:
		ndiag = 7;
		
		rate_title = (char **)malloc(ndiag*sizeof(char *));
		rate_title[0]= "MCC Elastic";
		rate_title[1]= "MCC Excitation1";
		rate_title[2]= "MCC Excitation2";
		rate_title[3]= "MCC Excitation3";
		rate_title[4]= "MCC Excitation4";
		rate_title[5]= "MCC Ionization";
		rate_title[6]= "MCC Charge X";
		break;
		
	case OXYGEN:
		ndiag = 21;
		
		rate_title = (char **)malloc(ndiag*sizeof(char *));
		rate_title[0]=  "MCC Elastic";
		rate_title[1]=  "MCC Rotational";
		rate_title[2]=  "MCC Vibrational1";
		rate_title[3]=  "MCC Vibrational2";
		rate_title[4]=  "MCC Vibrational3";
		rate_title[5]=  "MCC Vibrational4";
		rate_title[6]=  "MCC O2 Sing Delta";
		rate_title[7]=  "MCC O2 Sing Sigma";
		rate_title[8]=  "MCC Diss Attach";
		rate_title[9]=  "MCC E Loss1";
		rate_title[10]= "MCC E Loss2";
		rate_title[11]= "MCC E Loss3";
		rate_title[12]= "MCC E Loss4";
		rate_title[13]= "MCC Ionization";
		rate_title[14]= "MCC E Loss5";
		rate_title[15]= "MCC Diss Recombine";
		rate_title[16]= "MCC E Impact Detach";
		rate_title[17]= "MCC O- Elastic";
		rate_title[18]= "MCC O- Detachment";
		rate_title[19]= "MCC Mutual Neutral";
		rate_title[20]= "MCC Charge X";
		break;

	case MCC:
		ndiag = 5;
		
		rate_title = (char **)malloc(ndiag*sizeof(char *));
		rate_title[0]= "MCC Elastic";
		rate_title[1]= "MCC Excitation";
		rate_title[2]= "MCC Ionization";
		rate_title[3]= "MCC Scattering";
		rate_title[4]= "MCC Charge X";
		break;
	}
  
  /*  Allocate arrays for diagnostic rates */
  mccrate  = (SCALAR **)malloc(ndiag*sizeof(SCALAR *));
  rate_show= (SCALAR **)malloc(ndiag*sizeof(SCALAR *));
  for(i=0; i<ndiag; i++) {
    mccrate[i]  = (SCALAR *)malloc(ng*sizeof(SCALAR));
    rate_show[i]= (SCALAR *)malloc(ng*sizeof(SCALAR));
    for (j=0; j<ng; j++) mccrate[i][j] = rate_show[i][j] = 0.0;
  }
}
