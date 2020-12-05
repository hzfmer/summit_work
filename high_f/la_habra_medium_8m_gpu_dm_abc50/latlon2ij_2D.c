#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <queue>
using namespace std;

double dist_sphere(double lat1, double lon1, double lat2, double lon2, 
                   double R);

float mindist(int npx, int npy, float **lon, float **lat, float tlon, float tlat, int & mdx, int & mdy) {
   float min_dist;
   float md=1.e13;
   float R=6378137; /*Earth's radius im m*/
   int d[8][2] = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,1},{1,-1},{1,0},{1,1}};

   queue<pair<int, int>> q;
   vector<vector<bool>> visited(npy, vector<bool>(npx, false));
   q.push(make_pair(mdx, mdy));
   visited[mdy][mdx] = true;

   while(!q.empty()) {
       int x = q.front().first;
       int y = q.front().second;
       min_dist = dist_sphere(lat[y][x], lon[y][x], tlat, tlon, R);
       int count = 0;
       for(int i=0; i<8; ++i) {
           int x1 = x + d[i][0];
           int y1 = y + d[i][1];
           if(!visited[y1][x1] && x1>=0 && x1<npx && y1 >= 0 && y1 <=npy) {
                md = dist_sphere(lat[y1][x1], lon[y1][x1], tlat, tlon, R);
                if(md <= min_dist) {
                    count++;
                    q.push(make_pair(x1, y1));
                }
                visited[y1][x1] = true;
           }
       }
       if(count==0) {
           mdx = x;
           mdy = y;
           return md;
       }
   }
}

/*
Input:
------
    NX, NY: int
        Number of grid
    grid_fname: string
        Name of grid FILE
    stat_fname: string
        Name of station ascii file with three columns [site_name, lat, lon]
    output_fname: string
        Name of output ascii file with three columns [site_name, idx, idy]
*/
int main(int argc, char * argv[]){
    if(argc < 6) {
        printf("Inputs: NX, NY, grid_fname, stat_fname, output_fname\nAborting!\n");
        return 0;
    }
    int NX = atoi(argv[1]);
    int NY = atoi(argv[2]);
    char *grid_fname = argv[3];
    char *stat_fname = argv[4];
    char *output_fname = argv[5];

    long int NP, k;
    int mdx=0, mdy=0;  // Initial trial 
 
    double *buff;
    float **lat, **lon;
    FILE *fid_grid, *fid_in, *fid_out;
    float tlat, tlon, md;
    char cname[15];
    int line;
 
    // Read the mesh grid
    NP = (long int) NX * NY;
    buff=(double*) calloc(NP*3, sizeof(double));
    lat=(float**) calloc(NP, sizeof(float));
    lon=(float**) calloc(NP, sizeof(float));
 
    fprintf(stdout, "Reading mesh...");
    fflush(stdout);
    fid_grid = fopen(grid_fname, "r");
    fread(buff, NP*3, sizeof(double), fid_grid);
    for (int i=0; i<NY; i++) {
        for (int j=0; j<NX; j++) {
            k = i * NX + j;
            lon[i][j] = (float) buff[k * 3];
            lat[i][j] = (float) buff[k * 3 + 1];
        }
    }
    free(buff);
    fclose(fid_grid);
    fprintf(stdout, "%f %f\n", lon[NY-1][NX-1], lat[NY-1][NX-1]);
    fprintf(stdout, " ok.\n");
 
    // Main loop
    fid_in = fopen(stat_fname, "r");
    fid_out = fopen(output_fname, "w");
    if(fid_in == NULL)
    {
        fprintf(stdout, "No input station file\n");
    }
    else {
        line = fscanf(fid_in, "%s %f %f\n", cname, &tlon, &tlat);
        while(line != EOF) {
            md = mindist(NX, NY, lon, lat, tlon, tlat, mdx, mdy);
            if(md > 10) {
                fprintf(stdout, "Nearest distance: %f > 10km, probably wrong?\n", md);
            }
            fprintf(fid_out, "%s %d %d", cname, mdx, mdy);
            fflush(fid_out);
            line = fscanf(fid_in, "%s %f %f\n", cname, &tlon, &tlat);
        }
    }
    free(lon);
    free(lat);
    fclose(fid_in);
    fclose(fid_out);
    fprintf(stdout, " - finished.\n");
    return(0);
 }
