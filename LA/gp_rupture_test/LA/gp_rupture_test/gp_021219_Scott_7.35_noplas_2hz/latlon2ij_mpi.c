#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define MIN(x, y) (((x)) < ((y)) ? (x):(y))
float dist_sphere(float lat1, float lon1, float lat2, float lon2, 
            float R);

float * mindist(long int npt, int nxf, float *lon, float *lat, float * tlon, float * tlat, int* x_idx, int * y_idx){
    long int k;
    long int k0, k1;
    int m, x1, y1, x2, y2, nx=6320;
    float adist;
    float md[nxf];
    float R=6378.137; /*Earth's radius im km*/
    FILE * fid;
    for (m=0; m<nxf; m++) {
        if (m==0) {
            fid = fopen("fault_loc.idx", "r");
            fscanf(fid, "%d %d\n", &x1, &y1);
            fclose(fid);
            k0 = (long int) nx * (y1 - 30) + x1 - 100;
            k1 = (long int) nx * (y1 + 30) + x1 + 100;
        }
        else {
            k0 = (long int) nx * (y_idx[m-1] - 10) + x_idx[m-1] - 10;
            k1 = (long int) nx * (y_idx[m-1] + 10) + x_idx[m-1] + 10;
        }
        md[m] = 1.e10;
        for (k=k0; k<k1; k++) {
            adist = dist_sphere(lat[k], lon[k], tlat[m], tlon[m], R);
            if (adist < md[m]) {
                md[m] = adist;
                x_idx[m] = k % nx;
                y_idx[m] = k / nx;
            }
        }
        if (m % 500 == 0) {
            fprintf(stdout, "%f %f %f %f\n", lat[k], lon[k], tlat[m], tlon[m]);
            fprintf(stdout, "%d / %d column computed\n", m, nxf);
        }
    }
    return md;
}


int main (int argc, char *argv[]){
    long int npf, npt, m, k, n, * idx;
    int rank, nprocs;
    float *buff;
    float *buff2;
    float *lat, *lon, ** tlat, ** tlon;
    FILE *fid, *fid2, *fid_mesh;
    float * md;
    int xi, yi;
    // int *xi, *yi;
    // int ** ll;
    int ** x_idx, ** y_idx;
    int nx = 6320, ny = 4200, nxf = 980, nzf = 240; 
    int err;
    int number_amount;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Status status;
    fprintf(stdout, "rank = %d, nprocs = %d\n", rank, nprocs); 
    npt = (long int) nx * ny;
    npf = (long int) nxf * nzf;

    buff=(float *) calloc(npf * 2, sizeof(float));
    tlat=(float **) calloc(nzf, sizeof(float *));
    tlon=(float **) calloc(nzf, sizeof(float *));
    x_idx = (int **) calloc(nzf, sizeof(int *));
    y_idx = (int **) calloc(nzf, sizeof(int *));
    // ll = (int **) calloc(npf, sizeof(int *)); 
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
    for (k=0; k<nzf; k++) {
        tlon[k] = (float *) calloc(nxf, sizeof(float));
        tlat[k] = (float *) calloc(nxf, sizeof(float));
        for (m=0; m<nxf; m++) {
            n = k * nxf + m;
            tlon[k][m] = (float) buff[n * 2];
            tlat[k][m] = (float) buff[n * 2 + 1];
        }
    }
    free(buff);
    fclose(fid);
    fprintf(stdout, "Fault location read!\n");
    fflush(stdout);
    
    fid_mesh = fopen("surf.grid", "r");
    if (rank==0) {
        fprintf(stdout, "Reading mesh...\n");
    }
    ret = fread(buff2, sizeof(float), npt * 2, fid_mesh);
    fprintf(stdout, "read %zu bytes\n", ret*sizeof(float));
    for (k=0; k<npt; k++) {
        lon[k] = buff2[k * 2];
        lat[k] = buff2[k * 2 + 1];
    }
    free(buff2);
    fclose(fid_mesh);
    fprintf(stdout, "Mesh read!\n&lon[0] is %lu\n lon[2654000+rank]=%f\n", &lon[0], lon[2654000+rank]);
    fflush(stdout);

    for (m=rank; m<nzf; m+=nprocs){
        md = (float *) calloc(nxf, sizeof(float));
        x_idx[m] = (int *) calloc(nxf, sizeof(int));
        y_idx[m] = (int *) calloc(nxf, sizeof(int));
        md = mindist(npt, nxf, lon, lat, tlon[m], tlat[m], x_idx[m], y_idx[m]);
        fprintf(stdout, "Computed %ld / %d, md=%f\n", m, nzf, md[0]);
        fflush(stdout);
        if (rank == 0) {
            for (k=1; k<nprocs; k++){
                x_idx[m+k] = (int *) calloc(nxf, sizeof(int));
                y_idx[m+k] = (int *) calloc(nxf, sizeof(int));
                MPI_Recv(&(x_idx[m+k][0]), nxf, MPI_INT, k, 101, MPI_COMM_WORLD, &status);
                MPI_Recv(&(y_idx[m+k][0]), nxf, MPI_INT, k, 102, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_FLOAT, &number_amount);
            }
        } 
        else {
            // x_idx[m] = (int *) calloc(nxf, sizeof(int));
            // y_idx[m] = (int *) calloc(nxf, sizeof(int));
            MPI_Send(&(x_idx[m][0]), nxf, MPI_INT, 0, 101, MPI_COMM_WORLD);
            MPI_Send(&(y_idx[m][0]), nxf, MPI_INT, 0, 102, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        fid2 = fopen("fault_idx.bin", "wb");
        for (m=0; m<nzf; m++) {
            for (k=0; k<nxf; k++) {
            // fprintf(fid2, "%d %d\n", xi[k], yi[k]);
                fwrite(&x_idx[m][k], sizeof(int), 1, fid2);
                fwrite(&y_idx[m][k], sizeof(int), 1, fid2);
            }
            free(x_idx[m]);
            free(y_idx[m]);
        }
        fclose(fid2);
        fprintf(stdout, "Finished writing indecis of the whole fault\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return(0);
}
