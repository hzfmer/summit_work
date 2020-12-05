#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

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


int main (int argc, char *argv[]){
    int npt, m, k, n;
    int np, idx=-1;
    int rank, nprocs;
    float *buff;
    float *buff2;
    float *lat, *lon, *tlat, *tlon;
    FILE *fid, *fid2, *fid_mesh;
    float *md;
    int xi, yi;
    int nx=6320, ny=4200; 
    int err;
    int number_amount;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Status status;
    fprintf(stdout, "rank = %d, nprocs = %d\n", rank, nprocs); 
    np = 5116;  /* Twice the length of fault */
    npt = (int) nx * ny;
     
    buff=(float*) calloc(np * 2, sizeof(float));
    buff2 = (float*) calloc(npt * 2, sizeof(float));
    lat=(float*) calloc(np, sizeof(float));
    lon=(float*) calloc(np, sizeof(float));
    tlat = (float*) calloc(npt, sizeof(float));
    tlon = (float*) calloc(npt, sizeof(float));
    md = (float*) calloc(npt, sizeof(float));
   
    // MPI_Bcast is too slow!
    // Just read in the mesh for each rank
    fid = fopen("fault_full_loc.bin", "r");
    fread(buff, sizeof(float), np * 2, fid);
    if (rank==0) {
        fprintf(stdout, "Read fault location\n");
    }
    for (k=0; k<np; k++) {
        lon[k] = (float) buff[k*2];
        lat[k] = (float) buff[k*2+1];
    }
    free(buff);
    fclose(fid);
    
    fid_mesh = fopen("surf.grid", "r");
    if (rank==0) {
        fprintf(stdout, "Reading mesh...\n");
    }
    size_t ret = fread(buff2, sizeof(float), npt * 2, fid_mesh);
    fprintf(stdout, "read %zu bytes\n", ret*sizeof(float));
    for (k=0; k<npt; k++) {
        if (k % 10000000 == 0) {
            fprintf(stdout, "Rank=%d, Read %d / %d\n", rank, k, npt);
            fflush(stdout);
        }
        tlon[k] = buff2[k * 2];
        tlat[k] = buff2[k * 2 + 1];
    }
    free(buff2);
    fclose(fid_mesh);
    fprintf(stdout, "&tlon is %lu\n tlon[2654000+rank]=%f\n", &tlon, tlon[2654000+rank]);
    fflush(stdout);

    if (rank == 100000) {
        fid = fopen("fault_full_loc.bin", "r");
        fread(buff, sizeof(float), np * 2, fid);
        fprintf(stdout, "Read fault location\n");
        for (k=0; k<np; k++) {
            lon[k] = (float) buff[k*2];
            lat[k] = (float) buff[k*2+1];
        }
        free(buff);
        fclose(fid);
        
        fid_mesh = fopen("surf.grid", "r");
        fprintf(stdout, "Reading mesh...\n");
        size_t ret = fread(buff2, sizeof(float), npt * 2, fid_mesh);
        fprintf(stdout, "read %zu bytes\n", ret*sizeof(float));
        for (k=0; k<npt; k++) {
            if (k % 1000000 == 0) {
                fprintf(stdout, "Read %d / %d\n", k, npt);
                fflush(stdout);
            }
            tlon[k] = buff2[k * 2];
            tlat[k] = buff2[k * 2 + 1];
        }
        free(buff2);
        fclose(fid_mesh);
        fprintf(stdout, "&tlon is %lu\n tlon[26543999]=%f\n", &tlon, tlon[26543999]);
        fprintf(stdout, "Done reading\n"); fflush(stdout);
       
        err = MPI_Bcast(lon, np, MPI_FLOAT, 0, MPI_COMM_WORLD);
        if (err != MPI_SUCCESS) {
            fprintf(stdout, "Failed to MPI_Bcast\n"); 
        }
        MPI_Bcast(lat, np, MPI_FLOAT, 0, MPI_COMM_WORLD);
        fprintf(stdout, "Done bcast\n");  fflush(stdout);

        MPI_Bcast(tlon, npt, MPI_FLOAT, 0, MPI_COMM_WORLD);
        err = MPI_Bcast(tlat, npt, MPI_FLOAT, 0, MPI_COMM_WORLD);
        if (err != MPI_SUCCESS) {
            fprintf(stdout, "Failed to MPI_Bcast!\n");
        }
        else {
            fprintf(stdout, "Finish bcast\n");
        }
        fflush(stdout);
    }

    // MPI_Barrier(MPI_COMM_WORLD);

    for (m=rank; m<npt; m+=nprocs){
        md[m] = mindist(np, lon, lat, tlon[m], tlat[m]);
        if (rank == 0) {
            for (k=1; k<nprocs; k++){
                MPI_Recv(&md[m+k], 1, MPI_FLOAT, k, 100, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_FLOAT, &number_amount);
            }
            if (m % 500000 == 0) {
                fprintf(stdout, "md[%d]=%f: Status.count=%d, msg source = %d\n", m, md[m], number_amount, status.MPI_SOURCE);
            fflush(stdout);
            }
        } 
        else {
            MPI_Send(&md[m], 1, MPI_FLOAT, 0, 100, MPI_COMM_WORLD);  
            if (m % 500001 == 0) {
                fprintf(stdout, "md[%d]=%f\n", m, md[m]);
                fflush(stdout);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        fid2=fopen("mesh_rjb.bin", "wb");
        fprintf(stdout, "sizeof(md): %d\n", sizeof(md));
        fwrite(md, npt * sizeof(float), 1, fid2); 
        fclose(fid2);
        fprintf(stdout, " - finished.\n");
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return(0);
}
