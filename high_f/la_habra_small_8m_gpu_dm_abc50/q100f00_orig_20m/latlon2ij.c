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
   int npt=351, m;
   double *buff;
   float *lat, *lon;
   FILE *fid, *fid2;
   float tlat, tlon, md;
   int xi, yi;
   int nx=1400, ny=1400; 
   char cname[10];

   np = (long int) nx * ny;
   buff=(double*) calloc(np*3, sizeof(double));
   lat=(float*) calloc(np, sizeof(float));
   lon=(float*) calloc(np, sizeof(float));

   fprintf(stdout, "Reading mesh...");
   fflush(stdout);
   fid=fopen("surf.grid", "r");
   fread(buff, np*3, sizeof(double), fid);
   for (k=0; k<np; k++) {
      lon[k] = (float) buff[k*3];
      lat[k] = (float) buff[k*3+1];
   }
   free(buff);
   fclose(fid);
   fprintf(stdout, "%f %f\n", lon[np-1], lat[np-1]);
   fprintf(stdout, " ok.\n");

   fid=fopen("../la_habra_large_statlist.txt", "r");
   fid2=fopen("la_habra_small_statlist_20m.idx", "w");
   for (m=0; m<npt; m++){
      fprintf(stdout, "\rProcessing station %d of %d", m+1, npt);
      fflush(stdout);
      fscanf(fid, "%s %f %f\n", cname, &tlon, &tlat);
   
      md=mindist(np, lon, lat, tlon, tlat, &idx);
      xi = idx % nx;
      yi = idx / nx;
      if (0 < xi && xi < nx - 1 && 0 < yi && yi < ny - 1) {
        fprintf(fid2, "%s %d %d\n", cname, xi, yi);
        fflush(fid2);
      }
   }
   fclose(fid);
   fclose(fid2);
   fprintf(stdout, " - finished.\n");
   return(0);
}