#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

double dist_sphere(double lat1, double lon1, double lat2, double lon2, 
                   double R);

float mindist(long int np, float *lon, float *lat, float tlon, float tlat, long int *mdi){
   int k;
   float dx, dy, adist;
   float md=1.e13;
   float R=6378137; /*Earth's radius im m*/

   for (k=0; k<np; k++){
      adist = dist_sphere(lat[k], lon[k], tlat, tlon, R);
      if (adist < md){
         md=adist;
         *mdi=k;
      }
   }
   return(md);
}

int main(){
   long int np, k, n, idx=-1;
   int npt=1, m;
   double *buff;
   float *lat, *lon;
   FILE *fid, *fid2;
   float tlat, tlon, md;
   int xi, yi;
   int nx=9000, ny=6750; 
   char cname[10];

   np = (long int) nx * ny;
   buff=(double*) calloc(np*3, sizeof(double));
   lat=(float*) calloc(np, sizeof(float));
   lon=(float*) calloc(np, sizeof(float));

   fprintf(stdout, "Reading mesh...");
   fflush(stdout);
   fid=fopen("la_habra_large.grid", "r");
   fread(buff, np*3, sizeof(double), fid);
   for (k=0; k<np; k++) {
      lon[k] = (float) buff[k*3];
      lat[k] = (float) buff[k*3+1];
   }
   free(buff);
   fclose(fid);
   fprintf(stdout, "%f %f\n", lon[np-1], lat[np-1]);
   fprintf(stdout, " ok.\n");

   fid=fopen("stat.ll", "r");
   fid2=fopen("stat.idx", "w");
   for (m=0; m<npt; m++){
      fprintf(stdout, "\rProcessing station %d of %d", m+1, npt);
      fflush(stdout);
      fscanf(fid, "%f %f %s\n", &tlon, &tlat, cname);
   
      md=mindist(np, lon, lat, tlon, tlat, &idx);
      xi = idx % nx;
      yi = idx / nx;
      fprintf(fid2, "%d %d %s\n", xi, yi, cname);
      fflush(fid2);
   }
   fclose(fid);
   fclose(fid2);
   fprintf(stdout, " - finished.\n");
   return(0);
}
