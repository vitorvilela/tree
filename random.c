#include <math.h>
#include "random.h"
 
 
//STATISTICS' FUNCTIONS - BEGIN

float statistic_randomOne(int *idum) {
  
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0 || !iy) { 								// Initialize.
  if (-(*idum) < 1) *idum=1; 								// Be sure to prevent idum = 0.
  else *idum = -(*idum);
  for (j=NTAB+7;j>=0;j--) { 								// Load the shuffle table (after 8 warm-ups).
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  if (j < NTAB) iv[j] = *idum;
  }
  iy=iv[0];
  }
  k=(*idum)/IQ; 									// Start here when not initializing.
  *idum=IA*(*idum-k*IQ)-IR*k; 								// Compute idum=(IA*idum) % IM without overif
  if(*idum < 0) *idum += IM; 								// flows by Schrage’s method.
  j=iy/NDIV_statistic_randomOne; 										// Will be in the range 0..NTAB-1.
  iy=iv[j]; 										// Outransfut previously stored value and refill the
  iv[j] = *idum; 									// shuffle table.
  if ((temp=AM_statistic_randomOne*iy) > RNMX) return RNMX; 						// Because users don’t expect endpoint values.
  else return temp;
    
}


float statistic_randomTwo(long int *idum) {

  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0) { 									//Initialize.
    if (-(*idum) < 1) *idum=1;								//Be sure to prevent idum = 0.
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { 								//Load the shuffle table (after 8 warm-ups).
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1; 									//Start here when not initializing.
  *idum=IA1*(*idum-k*IQ1)-k*IR1; 							//Compute idum=(IA1*idum) % IM1 without
  if (*idum < 0) *idum += IM1; 								//overflows by Schrage’s method.
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; 							//Compute idum2=(IA2*idum) % IM2 likewise.
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV; 										//Will be in the range 0..NTAB-1.
  iy=iv[j]-idum2; 									//Here idum is shuffled, idum and idum2 are
  iv[j] = *idum; 									//combined to generate outransfut.
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;							// Because users don’t expect endpoint values.
  else return temp;

}  


float statistic_gasdev(int *idum) {

  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  if (*idum < 0) iset=0; 
  if (iset == 0) {
    do {
      v1=2.0*statistic_randomOne(idum)-1.0; 
      v2=2.0*statistic_randomOne(idum)-1.0; 
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0); 
    fac=sqrt(-2.0*log(rsq)/rsq);
 gset=v1*fac;
    iset=1;
    return v2*fac;
  } 
  else { 
    iset=0; 
    return gset;    
  } 
  
}


//STATISTICS' FUNCTIONS - END
