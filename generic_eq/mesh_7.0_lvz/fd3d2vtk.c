#include <stdio.h>
#include <math.h>

int main () {
   int nx=740;
   int ny=450;
   int nz=250;
   int step=1, stepx=20;

   int nx1, ny1, nz1;

   int mod_tot;
  
   int n, m, l;

   float dh=200.;
 
   char errmsg[200];

   char *vname="mojave_100m_7.0_lvz2.vtk";
   char *fname ="mojave_100m_7.0_lvz2";

   FILE *vfid;
   FILE *fid;

   nx1=(int) ceil( (double) nx/stepx);
   ny1=(int) ceil( (double) ny/step);
   nz1=(int) ceil( (double) nz/step);

   vfid=fopen(vname, "w");
   if (vfid==NULL) {
      sprintf(errmsg, "could not open %s", vname);
      perror(errmsg);
      return(-1);
   }

   fprintf(vfid, "# vtk DataFile Version 3.0\n");
   fprintf(vfid, "created from %s\n", fname);
   fprintf(vfid, "ASCII\n");
   fprintf(vfid, "DATASET STRUCTURED_POINTS\n");
   fprintf(vfid, "DIMENSIONS %d %d %d\n", nx1, ny1, nz1);
   fprintf(vfid, "ORIGIN 0 0 0 \n");
   fprintf(vfid, "SPACING %f %f %f\n", dh*stepx, dh*step, -dh*step);
   fprintf(vfid, "POINT_DATA %d\n", nx1 * ny1 * nz1);
   fprintf(vfid, "SCALARS Beta float\n");
   fprintf(vfid, "LOOKUP_TABLE default\n");

   struct {
      int xi, yi, zi;
      float alpha, beta, rho, Qp, Qs;
   } node_par;

   fid=fopen(fname, "r");
   if (fid==NULL) {
      sprintf(errmsg, "could not open %s", fname);
      perror(errmsg);
      return(-1);
   }
   
   for (l=0; l<nz; l++){
     for (m=0; m<ny; m++){
        for (n=0; n<nx; n++){
           fread(&node_par, sizeof(node_par), 1, fid);
           mod_tot = (l % step ) + (m % step) + (n % stepx);
           if (mod_tot == 0){ 
               //fwrite(&node_par.beta, sizeof(float), 1, vfid);
               fprintf(vfid, "%e\n", node_par.beta);
           }
        }
      }
   }

   fclose(vfid);
   fclose(fid);
   return(0);
}
