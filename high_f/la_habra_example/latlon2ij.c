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
   float md=1.e16;
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
   long int np=1960000, k, n, idx=-1;
   int npt=35, m;
   float *lat, *lon;
   FILE *fid, *fid2;
   char newline[200], newline2[200];
   char *sep="|", *tmp, *latchar, *lonchar;
   float tlat, tlon, md;
   int xi, yi;
   float dh=20.;
   int nx=1400, ny=1400; 
   time_t t0, t1;
//   char cname[10];

   lat=(float*) calloc(np, sizeof(float));
   lon=(float*) calloc(np, sizeof(float));

   fprintf(stdout, "Reading mesh...");
   fflush(stdout);
   fid=fopen("mesh_20m.ll", "r");
   for (k=0; k<np; k++) fscanf(fid, "%f %f\n", lon+k, lat+k);
   fclose(fid);
   fprintf(stdout, " ok.\n");

   fid=fopen("merged_lahabra_ll.txt", "r");
   fid2=fopen("stat.idx", "w");
   for (m=0; m<npt; m++){
      fprintf(stdout, "\rProcessing station %d of %d", m+1, npt);
      fflush(stdout);
      fscanf(fid, "%f %f\n", &tlon, &tlat);
   
      md=mindist(np, lon, lat, tlon, tlat, &idx);
      fprintf(stdout, "station %d: %f, %f: min_dist = %f\n", m, tlon, tlat, md);
      if (md < dh){
         xi = idx % nx;
         yi = idx / nx;   
         fprintf(fid2, "%d %d\n", xi, yi);
         fflush(fid2);
      }
   }
   fclose(fid);
   fclose(fid2);
   fprintf(stdout, " - finished.\n");
   return(0);
}
