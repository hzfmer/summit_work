#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void errhandle(int ierr, char *where){
   int errlen;
   char *errstr;

   if (ierr != 0) {
      fprintf(stderr, "error in %s\n", where);
      errstr=calloc(500, sizeof(char));
      MPI_Error_string(ierr, errstr, &errlen);
      fprintf(stderr, "%s", errstr);
      MPI_Finalize();
   }
}

void readlayer(int nx, int ny, int zi, int nvar, MPI_File fh, float *buf){
     MPI_Offset off; 
     int np, ierr;

     off = (MPI_Offset) nx * (MPI_Offset) ny * (MPI_Offset) zi * (MPI_Offset) sizeof(MPI_FLOAT) * 
           (MPI_Offset) nvar ;

     np = nx * ny * nvar;

     ierr=MPI_File_read_at_all(fh, off, buf, np, MPI_FLOAT, MPI_STATUS_IGNORE);
     errhandle(ierr, "read_at_all");
}

void writelayer(int nx, int ny, int zi, int nvar, MPI_File fh, float *buf){
     MPI_Offset off; 
     int np, ierr;

     off = (MPI_Offset) nx * (MPI_Offset) ny * (MPI_Offset) zi * (MPI_Offset) sizeof(MPI_FLOAT) *
           (MPI_Offset) nvar ;

     np = nx * ny * nvar;

     ierr=MPI_File_write_at_all(fh, off, buf, np, MPI_FLOAT, MPI_STATUS_IGNORE);
     errhandle(ierr, "write_at_all");
}

void def_damage_zone(int nx, int ny, int zi, float dh, int nvar, int faultpos, float depth1, float depth2, 
     float dist1, float dist2, float vsfact, float vsmin, float *buf){

     int k,l;
     long int bp;
     float depth, fltdist, vs, rfact, zfact;
     int vspos;

     if (nvar > 3) vspos = 4;
     else vspos = 1;

     depth = dh * (float) zi;
     if (depth < depth1) zfact = vsfact;
     else if ((depth >= depth1) && (depth <= depth2))
        zfact = vsfact + (depth - depth1)/(depth2 - depth1) * (1.f - vsfact);
     else
        zfact = 1.0f;
   
     fprintf(stdout, "depth = %f, zfact = %f\n", depth, zfact);
 
     for (k=0; k<ny; k++){
        for (l=0; l<nx; l++){
           bp=(k*nx+l)*nvar + vspos;

           vs=buf[bp];

           fltdist = fabsf((float) k + 0.5f - (float) faultpos) * dh;
           if ((fltdist <= dist2) && (fltdist >= dist1))
              rfact = zfact + (fltdist - dist1) / (dist2 - dist1) * (1.f - zfact);

           else if (fltdist <= dist1) 
              rfact = zfact;
           else 
              rfact = 1.f;
           buf[bp]=fmaxf(vsmin, vs * rfact);
        }
     }
}

int main (int argc, char *argv[]) {

   int nx=3000, ny=1372, nz=800;
   int nvar=8;   
   float dh=100.f;

   float *buf;
   long int np;
   char *infile="SAF_dyn_100m+idx";
   char *outfile="SAF_dyn_100m+idx_dmz";
   int rank, ncpus;
   MPI_File infid, outfid;
   int ierr;
   int k;
   
   int faultpos=1172;
   float depth1=4000.f, depth2=6000.f;
   float dist1=225.f, dist2=750.f;    
   float vsfact=0.7f, vsmin=500.f;

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_size(MPI_COMM_WORLD,&ncpus);

   ierr=MPI_File_open(MPI_COMM_WORLD, infile, MPI_MODE_RDONLY, MPI_INFO_NULL, &infid);
   errhandle(ierr, "input file open");
   MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outfid);
   errhandle(ierr, "output file open");

   np = (long int) nx * ny * nvar;

   buf = (float*) calloc (np, sizeof(float));
 
   for (k=rank; k<nz; k+=ncpus) {
      fprintf(stdout, "rank=%d, k=%d\n", rank, k);
      readlayer(nx, ny, k, nvar, infid, buf);      
      def_damage_zone(nx, ny, k, dh, nvar, faultpos, depth1, depth2, dist1, dist2, vsfact, vsmin, buf);
      writelayer(nx, ny, k, nvar, outfid, buf);      
   }

   MPI_File_close(&infid);

   MPI_File_close(&outfid);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return(0);

} 
