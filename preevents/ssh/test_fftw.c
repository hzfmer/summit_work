#include <stdio.h>
#include <fftw3.h>

int main () {
     fftw_complex in[50], out[50];
     fftw_plan p;

     p = fftw_create_plan(50, FFTW_FORWARD, FFTW_ESTIMATE);
     return(0);
}
