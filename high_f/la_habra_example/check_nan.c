#include <stdio.h>
#include <stdlib.h>

int LEN = 4;
int count = 0;
int main(int argc, char *argv[])
{
    FILE *file;
    float num;
    file = fopen(argv[1],"rb");
    fread(&num, sizeof(num), 1, file);
    while(!feof(file))
    {
        count ++;
        if (num != num)
        {
          printf("NaN detected at: %d\n", count); 
          return 0;
        }
        fread(&num,sizeof(num),1,file);
// printf("count = %d, %f\n", count, num);
    }
    fclose(file);
    return 0;
}

