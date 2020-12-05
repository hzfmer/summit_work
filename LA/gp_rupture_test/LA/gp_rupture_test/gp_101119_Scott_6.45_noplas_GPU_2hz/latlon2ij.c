#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

float dist_sphere(float lat1, float lon1, float lat2, float lon2, 
                   float R);

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
   int npt=2, m;
   float *buff;
   float *lat, *lon;
   FILE *fid, *fid2;
   float tlat, tlon, md;
   int xi, yi;
   int nx=6320, ny=4200; 
   char cname[10];

   np = (long int) nx * ny;
   buff=(float*) calloc(np*2, sizeof(float));
   lat=(float*) calloc(np, sizeof(float));
   lon=(float*) calloc(np, sizeof(float));

   fprintf(stdout, "Reading mesh...");
   fflush(stdout);
   fid=fopen("surf.grid", "r");
   fread(buff, np*2, sizeof(float), fid);
   fprintf(stdout, "Mesh read");
   for (k=0; k<np; k++) {
      lon[k] = (float) buff[k*2];
      lat[k] = (float) buff[k*2+1];
   }
   fprintf(stdout, "latlon read");
   free(buff);
   fclose(fid);
   fprintf(stdout, "%f %f\n", lon[np-1], lat[np-1]);
   fprintf(stdout, " ok.\n");

   fid=fopen("fault_loc.txt", "r");
   fid2=fopen("fault_loc.idx", "w");
   for (m=0; m<npt; m++){
      fprintf(stdout, "\rProcessing fault %d / %d", m+1, npt);
      fflush(stdout);
      fscanf(fid, "%f %f\n", &tlon, &tlat);
      fprintf(stdout, "lon, lat = %f, %f\n", tlon, tlat);
   
      md=mindist(np, lon, lat, tlon, tlat, &idx);
      xi = idx % nx;
      yi = idx / nx;
      fprintf(fid2, "%d %d\n", xi, yi);
      fflush(fid2);
   }
   fclose(fid);
   fclose(fid2);
   fprintf(stdout, " - finished.\n");
   return(0);
}
