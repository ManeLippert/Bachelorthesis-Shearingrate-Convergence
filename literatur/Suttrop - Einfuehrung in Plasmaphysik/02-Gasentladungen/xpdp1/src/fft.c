#include <math.h>
#include "pdp1.h"

#define SWAP(a, b) tempr=(a); (a)=(b); (b)=tempr

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif



/***************************************************************/
/* Replace `data' by its discrete Fourier transform, if `isign'is
	input as 1, or by its inverse discrete Fourier transform, if 
	`isign' is input as -1. `data' is a complex array of length `nn',
	input as a real array data[1..2*nn]. `nn' MUST be an integer
	power of 2 (this is not checked for!?)  */

void four1(SCALAR data[], int nn, int isign)
{
  register int i, j, m, n;
  int  mmax, istep;
  SCALAR wtemp, wr, wpr, wpi, wi, theta;
  SCALAR tempr, tempi;
  
  n= nn<<1;
  j=1;
  
  for(i=1; i<n; i+=2)
  {
    if(j >i)
    {
      SWAP(data[j], data[i]);
      SWAP(data[j+1], data[i+1]);
    }
    m= n>>1;
    while(m >=2 &&j >m)
    {
      j -=m;
      m >>=1;
    }
    j +=m;
  }
  mmax =2;
  while(n> mmax)
  {
    istep= 2*mmax;
    theta= 6.28318530717959/(isign*mmax);
    wtemp= sin(.5*theta);
    wpr= -2.0*wtemp*wtemp;
    wpi= sin(theta);
    wr= 1.0;
    wi= 0.0;
    for(m=1; m<mmax; m+=2)
    {
      for(i=m; i<n; i+=istep)
      {
	j=i+mmax;
	tempr= wr*data[j]-wi*data[j+1];
	tempi= wr*data[j+1]+wi*data[j];
	data[j]= data[i]- tempr;
	data[j+1]= data[i+1]- tempi;
	data[i] +=tempr;
	data[i+1] +=tempi;
      }
      wr= (wtemp=wr)*wpr - wi*wpi+wr;
      wi= wi*wpr + wtemp*wpi + wi;
    }
    mmax= istep;
  }
}

/***************************************************************/
/* Calculates the Fourier transform of a set of 2n real-valued
	data points. Replaces `data' by the positive frequency half of
	its complex Fourier transform. The real-valued first and last
	components of the complex transform are returnedas elements
	data[1] and data[2] respectively. `n' MUST be a power of 2.
	This routine also calculates the inverse transform of a complex
	data array if it is the transform of real data. (Result in
	this case MUST be divided by `n'.)  */

void realft(SCALAR data[], int n, int isign)
{
  register int i, i1, i2, i3, i4, n2p3;
  SCALAR c1=0.5, c2, h1r, h1i, h2r, h2i;
  SCALAR wr, wi, wpr, wpi, wtemp, theta;
  
  theta= M_PI/(SCALAR) n;
  if(isign ==1)
  {
    c2= -0.5;
    four1(data, n, 1);
  }
  else
  {
    c2= 0.5;
    theta= -theta;
  }
  wtemp= sin(0.5*theta);
  wpr= -2.0*wtemp*wtemp;
  wpi= sin(theta);
  wr= 1.0+wpr;
  wi= wpi;
  n2p3= 2*n+3;
  
  for(i=2; i<= n/2; i++)
  {
    i4= 1+(i3=n2p3 -(i2=1 +(i1=i+i-1)));
    h1r= c1*(data[i1] +data[i3]);
    h1i= c1*(data[i2] -data[i4]);
    h2r= -c2*(data[i2] +data[i4]);
    h2i= c2*(data[i1] -data[i3]);
    
    data[i1]= h1r +wr*h2r -wi*h2i;
    data[i2]= h1i +wr*h2i +wi*h2r;
    data[i3]= h1r -wr*h2r +wi*h2i;
    data[i4]= -h1i +wr*h2i +wi*h2r;
    
    wr= (wtemp=wr)*wpr -wi*wpi +wr;
    wi= wi*wpr +wtemp*wpi +wi;
  }
  
  if(isign ==1)
  {
    data[1]= (h1r= data[1]) +data[2];
    data[2]= h1r - data[2];
  }
  else
  {
    data[1]= c1*((h1r=data[1]) +data[2]);
    data[2]= c1*(h1r -data[2]);
    four1(data, n, -1);
  }

	  /**djc**/
  /* Is it too much to ask to get the correctly normalized value of the forward 
     or backward tranforms? */
  for (i=1; i<=2*n; i++)
    {
      if (isign== 1) data[i] /= (2*n);
      if (isign==-1) data[i] *= 2;
    }
}

/***************************************************************/
