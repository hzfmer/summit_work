#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(){
    float lon, lat;
    char cname[10];
    FILE *fid;
    fid=fopen("la_habra_large_statlist.txt", "r");
    fscanf(fid, "%s %f %f\n", cname, &lon, &lat);
    fprintf(stdout, "%s %f %f", cname, lon, lat);
    fflush(stdout);
    fclose(fid);
    return(0);
}
