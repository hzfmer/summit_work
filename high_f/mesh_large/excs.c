#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

MPI_Aint compute_odc_offset(int nx, int ny, int nz, int nvar, 
    int xi, int yi, int zi){
   MPI_Aint offset;
   offset = (MPI_Aint) ( (MPI_Aint) zi*nx*ny + yi*nx + xi) * nvar * sizeof(float);
   if (offset < 0) fprintf(stderr, "Error: offset is %ld.\n", offset);
   return offset;
}

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

void extract_vs_section(int rank, int nprocs, int nvar, int nx, int ny, int nz,
   int i0, int i1, int j0, int j1, int k0, int k1, char *infile, char *outfile){

   int i, j, k, np, p;
   MPI_File infid, outfid;
   int *blocklen;
   MPI_Aint *map;
   MPI_Datatype filetype;
   MPI_Offset disp=0, disp2;
   float *buf, *vs;
   int ierr;

   np=(i1-i0)*(j1-j0);

   if (((k1-k0) % nprocs) != 0) {
        fprintf(stderr, "NZ of output is not divisible by NCPUS. Quitting\n");
        MPI_Finalize();
   }
   
   blocklen=(int*) calloc(np, sizeof(int));
   map=(MPI_Aint*) calloc (np, sizeof(MPI_Aint));
   buf=(float*) calloc(np*nvar, sizeof(float));
   vs=(float*) calloc(np, sizeof(float));

   MPI_File_open(MPI_COMM_WORLD, infile, MPI_MODE_RDONLY, 
       MPI_INFO_NULL, &infid);
   MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &outfid);

   for (k=k0+rank; k<k1; k+=nprocs){
      if (rank==0) fprintf(stdout, "reading slice %d\n", k);
      p=0;
      for (j=j0; j<j1; j++){
         for (i=i0; i<i1; i++){
             blocklen[p]=nvar;
             map[p]=compute_odc_offset(nx, ny, nz, nvar, i, j, k);
             p++;
         }
      }
      //this won't work for large files
      //ierr=MPI_Type_create_indexed_block(np, nvar, map, MPI_FLOAT, &filetype);
      ierr=MPI_Type_create_hindexed(np, blocklen, map, MPI_FLOAT, &filetype);
      errhandle(ierr, "create_indexed_block");
      ierr=MPI_Type_commit(&filetype);
      errhandle(ierr, "commit");
      ierr=MPI_File_set_view(infid, disp, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
      errhandle(ierr, "set_view");
      ierr=MPI_File_read_all(infid, buf, np*nvar, MPI_FLOAT, MPI_STATUS_IGNORE);
      errhandle(ierr, "read_all");

      if (nvar == 8){
	  p=0;
	  for (j=j0; j<j1; j++){
	     for (i=i0; i<i1; i++){
		 buf[p*8] = i - i0 + 1;
		 buf[p*8+1] = j - j0 + 1;
		 p++;
	     }
	  }
      }
      disp2= (MPI_Offset) (i1-i0)* (MPI_Offset) (j1-j0)* (MPI_Offset) (k-k0)*nvar*sizeof(float);
      if (disp2 < 0) {
         fprintf(stderr, "%d: disp2 = %ld.  Quiting\n", rank, disp2);
         MPI_Finalize();
      }

      ierr=MPI_File_write_at_all(outfid, disp2, buf, np*nvar, MPI_FLOAT, MPI_STATUS_IGNORE);
      errhandle(ierr, "write_at_all");
      ierr=MPI_Type_free(&filetype);
      errhandle(ierr, "Type_free");
   }

   MPI_File_close(&infid);
   MPI_File_close(&outfid);
   MPI_Barrier(MPI_COMM_WORLD);
   }

int main (int argc, char *argv[] ) {
   int rank, nprocs;
   int nx=9000, ny=6750, nz=3072;
   int nvar=3;
   int i0=5569, i1=5570, j0=0, j1=6750, k0=0, k1=3000;
   char *infile="la_habra_large_hetdm.media";
   char *outfile="la_habra_large_hetdm_cs_epi.media";

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   extract_vs_section(rank, nprocs, nvar, nx, ny, nz,
      i0, i1, j0, j1, k0, k1, infile, outfile);

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   return(0);

}
