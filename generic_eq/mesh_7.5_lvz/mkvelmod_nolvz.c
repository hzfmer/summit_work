#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define FLTPOS 250

void writelayer(int nx, int ny, int zi, float vp, float vs, float rho, float qp, 
   float qs, MPI_File fh){
     int k, l, np;
     float *buf;
     MPI_Offset off, bp; 

     np=nx*ny;
     buf=(float*) calloc(np*8, sizeof(float));
     off = np*8*sizeof(float)*zi;
 
     for (k=0; k<ny; k++){
        for (l=0; l<nx; l++){
           bp=(k*nx+l)*8;
           buf[bp]=l+1;
           buf[bp+1]=k+1;
           buf[bp+2]=zi+1;
           buf[bp+3]=vp;
           if ((k >= (FLTPOS - 2)) & (k <= (FLTPOS + 1))){
              buf[bp+4]=vs * 0.75;
              if (vs < 450.) vs = 450.;
           }
           else {
              buf[bp+4]=vs;
           }
           buf[bp+5]=rho;
           buf[bp+6]=qp;
           buf[bp+7]=qs;
        }
     }
     MPI_File_write_at_all(fh, off, buf, np*8, MPI_FLOAT, MPI_STATUS_IGNORE);
     free(buf);     
}

void writelayer2(int nx, int ny, int zi, float vp, float vs, float rho, float qp, 
   float qs, MPI_File fh){
     int k, l, np;
     float *buf;
     int *idx;
     MPI_Offset off, bp, offidx, offpar; 
     MPI_Datatype ftype1, ftype2;

     np=nx*ny;
     idx=(int*) calloc(np*3, sizeof(int));
     buf=(float*) calloc(np*5, sizeof(float));
     off = np*8*sizeof(float)*zi;
 
     for (k=0; k<ny; k++){
        for (l=0; l<nx; l++){
           bp=(k*nx+l)*3;
           idx[bp]=l+1;
           idx[bp+1]=k+1;
           idx[bp+2]=zi+1;

           bp=(k*nx+l)*5;
           buf[bp]=vp;
           if ((k >= (FLTPOS - 2)) & (k <= (FLTPOS + 1))){
              buf[bp+1]=vs * 0.75;
              if (vs < 450.) vs = 450.;
           }
           else {
              buf[bp+1]=vs;
           }
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
   int nx=1920, ny=500, nz=300;
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

   MPI_File_open(MPI_COMM_WORLD, "mojave_100m_7.5_lvz", 
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

