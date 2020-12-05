/* This program creates a kinematic source from results of a GP kinematic rupture 
   generator.

   Zhifeng Hu, Oct 2018 <zhh076@ucsd.edu>

   MPI-IO is used for both reading the srf file and writing the moment rate file.
   The moment rates are written to a single file, which must be partitioned using s
   srcpart
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "upsamp/upsamp.h"

void readpoints(int *xi, int *yi, int *zi, float *area, float *mu, 
        float *strike,float *dip, float *rake, int nflts) {
    FILE *fid;
    int k;
    fid = fopen("subfaults.idx", "r");
    if (fid==NULL) {
        perror("Could not open subfaults.idx");
        exit(-1);
    }

    for (k=0; k<nflts; k++){
        fscanf(fid, "%d %d %d %e %e %f %f %f\n", xi+k, yi+k, zi+k,
                area+k, mu+k, strike+k, dip+k, rake+k);
    }
    fclose(fid);
}

void errhandle(int ierr, char *where){
    int errlen;
    char *errstr;

    if (ierr!=0) {
        fprintf(stderr, "Error in %s\n", where);
        errstr=calloc(500, sizeof(char));
        MPI_Error_string(ierr, errstr, &errlen);
        fprintf(stderr, "%s", errstr);
        MPI_Finalize();
    }
}

/* Read STF from the file containing slip rates generated srf file */
void read_stf_mpiio(char *fname, int nx, int nz, int xi, int zi, int nt,
        float *x, int npoints){

    MPI_File fh;
    MPI_Datatype filetype;
    MPI_Offset disp;
    int ierr;

    ierr=MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    errhandle(ierr, "MPI_File_open()");

    ierr=MPI_Type_vector(nt, npoints, nx*nz, MPI_FLOAT, &filetype);
    errhandle(ierr, "MPI_Type_vector()");
    ierr=MPI_Type_commit(&filetype);
    errhandle(ierr, "MPI_Type_commi");

    disp = (MPI_Offset) sizeof(float) * ((MPI_Offset) nx*zi + (MPI_Offset) xi);
    ierr=MPI_File_set_view(fh, disp, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    errhandle(ierr, "MPI_File_set_view()");

    ierr=MPI_File_read_all(fh, x, nt*npoints, MPI_FLOAT, MPI_STATUS_IGNORE);
    errhandle(ierr, "MPI_File_read_all()");

    ierr=MPI_File_close(&fh);
    errhandle(ierr, "MPI_File_close()");

    ierr=MPI_Type_free(&filetype);
    errhandle(ierr, "MPI_Type_free()");
}

void write_momrate(int nst, int nx, int rank, int csize, int ion, int *xi, int *yi, int *zi,
        float **xx, float **yy, float **zz, float **xz, float **yz, float **xy){

    MPI_Offset offset;
    int *buf1, *blen1, *blen2;
    int ierr;
    int p, n, bpos;
    float *buf2;
    MPI_Aint extent;
    MPI_Aint *map1, *map2;
    MPI_File fh;
    MPI_Datatype filetype1, filetype2;
    
    MPI_File_open(MPI_COMM_WORLD, "momrate.dat", MPI_MODE_CREATE | MPI_MODE_WRONLY, 
            MPI_INFO_NULL, &fh);

    buf1 = (int*) calloc(3 * csize, sizeof(int));
    MPI_Type_extent(MPI_AINT, &extent);
    map1 = (MPI_Aint*) calloc(csize, extent);
    blen1 = (int*) calloc(csize, sizeof(int));

    for (p=0; p<csize; p++) {
        map1[p] = (MPI_Aint) (3*sizeof(int) + nst*6*sizeof(float)) * (MPI_Aint) + p;
        blen1[p] = 3;
        buf1[p*3] = xi[rank*nx+ion*csize+p];
        buf1[p*3+1] = yi[rank*nx+ion*csize+p];
        buf1[p*3+2] = zi[rank*nx+ion*csize+p];
    }

    offset = (MPI_Offset) (3*sizeof(int) + 6*nst*sizeof(float)) 
        * (MPI_Offset) (nx * rank + ion*csize);
    if (rank==0) fprintf(stdout, "offset = %lld\n", offset);
    if (rank==0) fprintf(stdout, "ion = %d, csize = %d\n", ion, csize);

    ierr=MPI_Type_create_hindexed(csize, blen1, map1, MPI_INT, &filetype1);
    errhandle(ierr, "MPI_Type_create_hindexed()");

    ierr=MPI_Type_commit(&filetype1);
    errhandle(ierr, "MPI_Type_commit()");

    ierr=MPI_File_set_view(fh, offset, MPI_INT, filetype1, "native", MPI_INFO_NULL);
    errhandle(ierr, "MPI_File_set_view()");

    ierr=MPI_File_write_all(fh, buf1, csize*3, MPI_INT, MPI_STATUS_IGNORE);
    errhandle(ierr, "MPI_File_write_all()");

    ierr=MPI_Type_free(&filetype1);
    errhandle(ierr, "MPI_Type_free()");

    free(buf1);
    free(map1);
    free(blen1);

    buf2 = (float*) calloc(csize*nst*6, sizeof(float));
    map2 = (MPI_Aint*) calloc(csize, extent);
    blen2 = (int*) calloc(csize, sizeof(int));

    for (p=0; p<csize; p++) {
        map2[p] = (MPI_Aint) (3*sizeof(int) + 6*nst*sizeof(float)) * (MPI_Aint) p + 3 * sizeof(int);
        blen2[p] = nst * 6;

        for (n=0; n<nst; n++) {
            bpos = p*nst*6 + n*6;
            buf2[bpos] = xx[p][n];
            buf2[bpos+1] = yy[p][n];
            buf2[bpos+2] = zz[p][n];
            buf2[bpos+3] = xz[p][n];
            buf2[bpos+4] = yz[p][n];
            buf2[bpos+5] = xy[p][n];
        }
    }

    ierr=MPI_Type_create_hindexed(csize, blen2, map2, MPI_FLOAT, &filetype2);
    errhandle(ierr, "MPI_Type_create_hindexed()");

    ierr=MPI_Type_commit(&filetype2);
    errhandle(ierr, "MPI_Type_commit()");

    ierr=MPI_File_set_view(fh, offset, MPI_FLOAT, filetype2, "native", MPI_INFO_NULL);
    errhandle(ierr, "MPI_File_set_view()");

    ierr=MPI_File_write_all(fh, buf2, csize*nst*6, MPI_FLOAT, MPI_STATUS_IGNORE);
    errhandle(ierr, "MPI_File_write_all()");

    MPI_Barrier(MPI_COMM_WORLD);

    ierr=MPI_Type_free(&filetype2);
    errhandle(ierr, "MPI_Type_free()");

    ierr=MPI_File_close(&fh);
    errhandle(ierr, "MPI_File_close()");

    free(buf2);
    free(map2);
    free(blen2);
}

float deg2rad(float deg) {
    float rad; 
    rad = deg / 180. * M_PI;
    return(rad);
}

void compmt(float str, float dip, float rake,
        float *xx, float *yy, float *zz, float *xz, float *yz, float *xy ){

      str=deg2rad(str);
      dip=deg2rad(dip);
      rake=deg2rad(rake);

      *yy= -(sinf(dip)*cosf(rake)*sinf(2.*str)+
           sinf(2.*dip)*sinf(rake)*sinf(str)*sinf(str));

      *xy= sinf(dip)*cosf(rake)*cosf(2.*str)+
           0.5*(sinf(2.*dip)*sinf(rake)*sinf(2.*str));

      *yz= (cosf(dip)*cosf(rake)*cosf(str)+
           cosf(2.*dip)*sinf(rake)*sinf(str));

      *xx= sinf(dip)*cosf(rake)*sinf(2.*str)-
           sinf(2.*dip)*sinf(rake)*cosf(str)*cosf(str);

      *xz= (cosf(dip)*cosf(rake)*sinf(str)-
           cosf(2.*dip)*sinf(rake)*cosf(str));

      *zz= sinf(2.*dip)*sinf(rake);
}

//#define MAG 6
int main(int argc, char *argv[]){
//*** M6.35 ***
    #if MAG == 6     // Mag 6.35
        int nx = 354, nz=212; /* size of fault */
        int nflts = 75048;
        int nio = 6; /* nio must divide nx */
//*** M7.35 ***
    #elif MAG == 7  // Mag 7.35
        int nx = 980, nz=240; /* size of fault */
        int nflts = 235200;
        int nio = 5; /* nio must divide nx */
//*** M8.45 ***
    #else           // Mag 8.45
        int nx = 5116, nz=220; /* size of fault */
        int nflts = 1125520;
        int nio = 4; /* nio must divide nx */
        fprintf(stdout, "Mag=8.45\n");
    #endif

    float dt1=0.05, dt2 = 0.005;
    int nt1 = 400, nt2 = 4000;
    
    int l, n, ki;
    float sla;

    float *buf, *x2;
    float **x1;
    
    float s_xx, s_yy, s_zz, s_xz, s_yz, s_xy;
    float **xx, **yy, **zz, **xz, **yz, **xy;

    int rank, ncpus;
    MPI_Status status;


    int *xi, *yi, *zi;
    float *area, *mu, *strike, *dip, *rake;

    int xip, zip;
    
    int p;
    int csize, c;

    csize = nx / nio;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ncpus);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    fprintf(stdout, "ncpus = %d\n", ncpus);
    
    if (nz!=ncpus) {
        if (rank==0) fprintf(stdout, "nz (%d) should equal ncpus (%d)\n", nz, ncpus);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    }

    x1=(float**) calloc(csize, sizeof(float*));

    xx=(float**) calloc(csize, sizeof(float*));
    yy=(float**) calloc(csize, sizeof(float*));
    zz=(float**) calloc(csize, sizeof(float*));
    xz=(float**) calloc(csize, sizeof(float*));
    yz=(float**) calloc(csize, sizeof(float*));
    xy=(float**) calloc(csize, sizeof(float*));

    for (p=0; p<csize; p++){
        x1[p] = (float*) calloc(nt1, sizeof(float));
        xx[p] = (float*) calloc(nt2, sizeof(float));
        yy[p] = (float*) calloc(nt2, sizeof(float));
        zz[p] = (float*) calloc(nt2, sizeof(float));
        xz[p] = (float*) calloc(nt2, sizeof(float));
        yz[p] = (float*) calloc(nt2, sizeof(float));
        xy[p] = (float*) calloc(nt2, sizeof(float));
    }

    buf = (float*) calloc(nt1*csize, sizeof(float));

    x2 = (float*) calloc(nt2, sizeof(float));

    MPI_Barrier(MPI_COMM_WORLD);

    xi = (int*) calloc(nflts, sizeof(int));
    yi = (int*) calloc(nflts, sizeof(int));
    zi = (int*) calloc(nflts, sizeof(int));

    area = (float*) calloc(nflts, sizeof(float));
    mu = (float*) calloc(nflts, sizeof(float));
    strike = (float*) calloc(nflts, sizeof(float));
    dip = (float*) calloc(nflts, sizeof(float));
    rake = (float*) calloc(nflts, sizeof(float));
    
    if (rank==0) {
        fprintf(stdout, "ion = %d, Done. \n Reading and broadcasting subfaults parameters ...", nio);
        fflush(stdout);
        readpoints(xi, yi, zi, area, mu, strike, dip, rake, nflts);
        #ifdef DEBUG
        for (l=0; l<nflts; l++) {
            fprintf(stdout, "%d %d %d\n", xi[l], yi[l], zi[l]);
        }
        #endif
    }

    MPI_Bcast(xi, nflts, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(yi, nflts, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(zi, nflts, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(area, nflts, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(mu, nflts, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(strike, nflts, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(dip, nflts, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(rake, nflts, MPI_FLOAT, 0, MPI_COMM_WORLD);

    zip = nz - rank -1;
    for (c=0; c<nio; c++) {
        xip = csize * c;
        fprintf(stdout, "%d, reading row %d, column %d\n", rank, zip, xip);
        // Care, it seems AWP stores fault properties upside-down
        // read_stf_mpiio("sliprate.bin", nx, nz, xip, zip, nt1, buf, csize);
        read_stf_mpiio("sliprate.bin", nx, nz, xip, rank, nt1, buf, csize);
        MPI_Barrier(MPI_COMM_WORLD);
      
        for (p=0; p<csize; p++) {
            for (n=0; n<nt1; n++) {
                x1[p][n] = buf[n*csize+p];
            }
        }

        for (p=0; p<csize; p++) {
            l = nx * rank + csize * c + p;
            upsamp(x1[p], dt1, dt2, nt1, nt2, x2);
            compmt(strike[l], dip[l], rake[l], &s_xx, &s_yy, &s_zz, &s_xz, &s_yz, &s_xy);

            for (n=0; n<nt2; n++) {
                sla=sqrt(pow(x2[n], 2.));
                sla *= area[l] * mu[l];
                xx[p][n] = s_xx * sla;
                yy[p][n] = s_yy * sla;
                zz[p][n] = s_zz * sla;
                xz[p][n] = s_xz * sla;
                yz[p][n] = s_yz * sla;
                xy[p][n] = s_xy * sla;
            }

        MPI_Barrier(MPI_COMM_WORLD);
        } /* end of p loop */
        MPI_Barrier(MPI_COMM_WORLD);

        if (rank==0) {
            fprintf(stdout, "Writing moment rates to disk ... \n");
            if (c==2) fprintf(stdout, "The 400th slip rate is %f\n", x1[0][400]);
        }
        write_momrate(nt2, nx, rank, csize, c, xi, yi, zi, xx, yy, zz, xz, yz, xy);
        MPI_Barrier(MPI_COMM_WORLD);
    } /* end of c loop */

MPI_Barrier(MPI_COMM_WORLD);
MPI_Finalize();
return(0);
}

















