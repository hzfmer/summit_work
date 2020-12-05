#include "upsamp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void *upsamp(float *S1, float dT1, float dT2, int npts1, int npts2, float *S2)
/* upsampling of seismograms by simple linear interpolation.
   not suitable for downsampling, since no anti-aliasing filter is
   implemented */
{
   int n, n_l, n_h;
   float t_n, t_l, t_h;

   if (dT2 > dT1)  
      {
      fprintf(stderr, "anti-aliasing filter not (yet) implemented.\n");
      fprintf(stderr, "Error: downsampling is not supported!\n");
      return(NULL);
      }

   for (n=0;n<npts2;n++)
     {
     t_n=dT2*n;
     n_l=floorf(t_n/dT1);
     n_h=ceilf(t_n/dT1);
     if (n_l == n_h)
      /* data is in NR type array, starting at index 1 ... */
        S2[n]=S1[n_l];
     else if (n_h > npts1-1)
        S2[n]=S1[n_l];
     else
        {
        t_l=n_l*dT1;
        t_h=n_h*dT1;
        S2[n]=S1[n_l] + (S1[n_h]-S1[n_l])/(t_h-t_l) * (t_n-t_l);
        }
     }
   return(0);
}
