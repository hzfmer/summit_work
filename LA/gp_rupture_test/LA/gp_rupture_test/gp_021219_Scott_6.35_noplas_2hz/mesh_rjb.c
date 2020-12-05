#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

float dist_sphere(float lat1, float lon1, float lat2, float lon2, 
                   float R);

float mindist(int np, float *lon, float *lat, float tlon, float tlat){
    int k;
    float dx, dy, adist;
    float md=1.e13;
    float R=6378.137; /*Earth's radius im m*/
    for (k=0; k<np; k++) {
        adist = dist_sphere(lat[k], lon[k], tlat, tlon, R);
        if (adist < md) {
            md=adist;
        }
    }
    return(md);
}

int main(){
   long int npt, m, k, n;
   int np, idx=-1;
   float *buff;
   double *buff2;
   float *lat, *lon, *tlat, *tlon;
   FILE *fid, *fid2, *fid_mesh;
   float *md;
   int xi, yi;
   int nx=6320, ny=4200; 
   char cname[10];

   np = 354;  /* Twice the length of fault */
   npt = (long int) nx * ny;
   buff=(float*) calloc(np*2, sizeof(float));
   buff2 = (double*) calloc(npt * 2, sizeof(double));
   lat=(float*) calloc(np, sizeof(float));
   lon=(float*) calloc(np, sizeof(float));
   tlat = (float*) calloc(npt, sizeof(float));
   tlon = (float*) calloc(npt, sizeof(float));
   md = (float*) calloc(npt, sizeof(float));

   fprintf(stdout, "Reading mesh...");
   fflush(stdout);
   fid_mesh = fopen("surf.grid", "r");
   fread(buff2, npt * 2, sizeof(double), fid_mesh);
   for (k=0; k<npt; k++) {
       tlon[k] = (float) buff2[k * 2];
       tlat[k] = (float) buff2[k * 2 + 1];
   }
   free(buff2);
   fclose(fid_mesh);
   
   fid=fopen("fault_full_loc.bin", "r");
   fread(buff, np * 2, sizeof(float), fid);
   for (k=0; k<np; k++) {
      lon[k] = (float) buff[k*2];
      lat[k] = (float) buff[k*2+1];
   }
   free(buff);
   fclose(fid);
   fprintf(stdout, "The last point of the fault: %f %f\n", lon[np-1], lat[np-1]);
   fprintf(stdout, " ok.\n");

   for (m=0; m<npt; m++) {
      md[m]=mindist(np, lon, lat, tlon[m], tlat[m]);
   }
   fid2=fopen("fault_distance.bin", "wb");
   fprintf(stdout, "sizeof(md): %d", sizeof(md));
   fwrite(md, npt * sizeof(float), 1, fid2); 
   fclose(fid2);
   fprintf(stdout, " - finished.\n");
   return(0);
}
