#include<stdlib.h>
#include<stdio.h>

#define MIN(x, y) (((x) < (y)) ? (x):(y))

int main(int argc, char *argv[]) {
    FILE *fid;
    int x1, x2, y1, y2;
    fid = fopen("fault_loc.idx", "r");
    fscanf(fid, "%d %d\n", &x1, &y1);
    fscanf(fid, "%d %d\n", &x2, &y2);
    fprintf(stdout, "%d %d\n", x1, y1);
    fprintf(stdout, "%d %d\n", x2, y2);
    fprintf(stdout, "min(x): %d\n", MIN(x1, x2));
    fprintf(stdout, "min(y): %d\n", MIN(y1, y2));
    fclose(fid);
    return 0;
}
