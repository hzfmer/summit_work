#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

//#define FLTPOS 40
#define FLTPOS 80
#define DH 100.

void writelayer2(int nx, int ny, int zi, float vp, float vs, float rho, float qp, 
   float qs, MPI_File fh){
     int k, l, np;
     float *buf;
     int *idx;
     MPI_Offset off, bp, offidx, offpar; 
     MPI_Datatype ftype1, ftype2;
     float fltdist, rfact, depth, zfact;

     np=nx*ny;
     idx=(int*) calloc(np*3, sizeof(int));
     buf=(float*) calloc(np*5, sizeof(float));
     off = np*8*sizeof(float)*zi;

     depth = zi * DH;
     if (depth < 4000.f) zfact = 0.70f;
     else if ((depth >= 4000.f) && (depth <= 6000.f))
        zfact = 0.7f + (depth - 4000.f)/2000.f* 0.3f;
     else
        zfact = 1.0f;
   
     fprintf(stdout, "depth = %f, zfact = %f\n", depth, zfact);
  
     for (k=0; k<ny; k++){
        for (l=0; l<nx; l++){
           bp=(k*nx+l)*3;
           idx[bp]=l+1;
           idx[bp+1]=k+1;
           idx[bp+2]=zi+1;

           bp=(k*nx+l)*5;
           buf[bp]=vp;
           fltdist = fabsf((float) k + 0.5f - (float) FLTPOS) * (float) DH;
           if ((fltdist <= 750.f) && (fltdist >= 225.f))
              rfact = zfact + (fltdist - 225.f) / (750.f - 225.f) * (1.f - zfact);

           else if (fltdist <= 225.f) 
              rfact = zfact;
           else 
              rfact = 1.f;
           buf[bp+1]=fmaxf(450.f, vs * rfact);
           buf[bp+2]=rho;
           buf[bp+3]=qp;
           buf[bp+4]=qs;
        }
     }

     MPI_Type_vector(nx*ny, 3, 8, MPI_INT, &ftype1);
     MPI_Type_commit(&ftype1);
     MPI_Type_vector(nx*ny, 5, 8, MPI_FLOAT, &ftype2);
     MPI_Type_commit(&ftype2);

     offidx=0;
     offpar=3*sizeof(int);

     MPI_File_set_view(fh, off+offidx, MPI_INT, ftype1, "native", MPI_INFO_NULL);
     MPI_File_write_all(fh, idx, np*3, MPI_INT, MPI_STATUS_IGNORE);

     MPI_File_set_view(fh, off+offpar, MPI_FLOAT, ftype2, "native", MPI_INFO_NULL);
     MPI_File_write_all(fh, buf, np*5, MPI_FLOAT, MPI_STATUS_IGNORE);

     free(buf);     
}

int main(int argc, char **argv){
   int rank, ncpus;
   MPI_File fh;
   float *vp, *vs, *d1, *qp, *qs;
   int nx=200, ny=180, nz=180;
   float dh=100.;
   int k;
   FILE *fid;
   int p;
   float dpt;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&ncpus);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   vp=(float*) calloc(nz, sizeof(float));
   vs=(float*) calloc(nz, sizeof(float));
   d1=(float*) calloc(nz, sizeof(float));
   qp=(float*) calloc(nz, sizeof(float));
   qs=(float*) calloc(nz, sizeof(float));

   if (rank==0){
      fid=fopen("../velmods/mojave_100m.dat", "r");
      for (p=0; p<nz; p++){
          fscanf(fid, "%f %e %e %e %e %e\n", &dpt, vp+p, vs+p, d1+p, qp+p, qs+p);
      }
      fclose(fid);
   }
   MPI_Bcast(vp, nz, MPI_FLOAT, 0, MPI_COMM_WORLD);
   MPI_Bcast(vs, nz, MPI_FLOAT, 0, MPI_COMM_WORLD);
   MPI_Bcast(d1, nz, MPI_FLOAT, 0, MPI_COMM_WORLD);
   MPI_Bcast(qp, nz, MPI_FLOAT, 0, MPI_COMM_WORLD);
   MPI_Bcast(qs, nz, MPI_FLOAT, 0, MPI_COMM_WORLD);

   //MPI_File_open(MPI_COMM_WORLD, "mojave_100m_toy_lvz2", 
   MPI_File_open(MPI_COMM_WORLD, "mojave_100m_toy2_lvz2", 
       MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

   for (k=rank; k<nz; k+=ncpus){
      fprintf(stdout, "%d, %f, %f\n", k, qp[k], qs[k]);
      writelayer2(nx, ny, k, vp[k], vs[k], d1[k], qp[k], qs[k], fh);
   }
   MPI_File_close(&fh);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return(0);

}

