#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

float dist_sphere(float lat1, float lon1, float lat2, float lon2, 
            float R);

float mindist(long int np, float *lon, float *lat, float tlon, float tlat, long int * mdi){
    long int k;
    float dx, dy, adist;
    float md=1.e13;
    float R=6378.137; /*Earth's radius im m*/

    // for (k=0; k<np; k++) {
    for (k=10744000; k<14536000; k++) {
        adist = dist_sphere(lat[k], lon[k], tlat, tlon, R);
        if (adist < md) {
            md=adist;
            *mdi = k;
        }
    }
    return(md);
}


int main (int argc, char *argv[]){
    long int npf, npt, m, k, n, idx=-1;
    int rank, nprocs;
    float *buff;
    float *buff2;
    float *lat, *lon, *tlat, *tlon;
    FILE *fid, *fid2, *fid_mesh;
    float md;
    int xi, yi;
    // int *xi, *yi;
    int ** ll;
    int nx=6320, ny=4200; 
    int err;
    int number_amount;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Status status;
    fprintf(stdout, "rank = %d, nprocs = %d\n", rank, nprocs); 
    npt = (long int) nx * ny;
    npf = (long int) 5116 * 220;

    buff=(float*) calloc(npf * 2, sizeof(float));
    tlat=(float*) calloc(npf, sizeof(float));
    tlon=(float*) calloc(npf, sizeof(float));
    ll = (int **) calloc(npf, sizeof(int *)); 
    // xi = (int*) calloc(npf, sizeof(int));
    // yi = (int*) calloc(npf, sizeof(int));
    buff2 = (float*) calloc(npt * 2, sizeof(float));
    lat = (float*) calloc(npt, sizeof(float));
    lon = (float*) calloc(npt, sizeof(float));
   
    // MPI_Bcast is too slow!
    // Just read in the mesh for each rank
    fid = fopen("fault_loc.bin", "r");
    size_t ret = fread(buff, sizeof(float), npf * 2, fid);
    fprintf(stdout, "read %zu bytes\n", ret*sizeof(float));
    if (rank==0) {
        fprintf(stdout, "Read fault location\n");
    }
    for (k=0; k<npf; k++) {
        tlon[k] = (float) buff[k*2];
        tlat[k] = (float) buff[k*2+1];
        /*if (tlat[k] < 1) {
            fprintf(stdout, "Error, tlat[%ld] < 1\n", k);
            getchar();
        }*/
    }
    free(buff);
    fclose(fid);
    
    fid_mesh = fopen("surf.grid", "r");
    if (rank==0) {
        fprintf(stdout, "Reading mesh...\n");
    }
    ret = fread(buff2, sizeof(float), npt * 2, fid_mesh);
    fprintf(stdout, "read %zu bytes\n", ret*sizeof(float));
    for (k=0; k<npt; k++) {
        if (k % 10000000 == 0) {
            fprintf(stdout, "Rank=%d, Read %ld / %ld\n", rank, k, npt);
            fflush(stdout);
        }
        lon[k] = buff2[k * 2];
        lat[k] = buff2[k * 2 + 1];
    }
    free(buff2);
    fclose(fid_mesh);
    fprintf(stdout, "&lon is %lu\n lon[2654000+rank]=%f\n", &lon, lon[2654000+rank]);
    fflush(stdout);

    for (m=rank; m<npf; m+=nprocs){
        md = mindist(npt, lon, lat, tlon[m], tlat[m], &idx);
        // xi[m] = idx % nx;
        // yi[m] = idx / nx;
        ll[m] = (int *) calloc(2, sizeof(int));
        ll[m][0] = idx % nx;
        ll[m][1] = idx / nx;
        if (rank == 0) {
            for (k=1; k<nprocs; k++){
                // MPI_Recv(&xi[m+k], 1, MPI_INT, k, 101, MPI_COMM_WORLD, &status);
                // MPI_Recv(&yi[m+k], 1, MPI_INT, k, 102, MPI_COMM_WORLD, &status);
                MPI_Recv(&ll[m][0], 2, MPI_INT, k, 102, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_FLOAT, &number_amount);
            }
            if (m % 5000 == 0) {
                fprintf(stdout, "Computed %ld / %ld, md=%f\n", m, npf, md);
                fflush(stdout);
            }
        } 
        else {
            // MPI_Send(&xi[m], 1, MPI_INT, 0, 101, MPI_COMM_WORLD);
            // MPI_Send(&yi[m], 1, MPI_INT, 0, 102, MPI_COMM_WORLD);
            MPI_Send(&ll[m][0], 2, MPI_INT, 0, 102, MPI_COMM_WORLD);
            if (m % 50001 == 0) {
                fprintf(stdout, "md[%ld] / %ld, md=%f\n", m, npf, md);
                fflush(stdout);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        fid2 = fopen("fault_loc.idx", "w");
        for (k=0; k<npf; k++) {
            // fprintf(fid2, "%d %d\n", xi[k], yi[k]);
            fwrite(ll[m], sizeof(int), 2, fid2);
        }
        fclose(fid2);
        fprintf(stdout, "Finished writing indecis of the whole fault\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return(0);
}
