#include <stdio.h>
#include <stdlib.h>

int LEN = 4;
int count = 1;

int main(int argc, char *argv[])
{
    FILE *file;
    file = fopen(argv[1], "rb");
    float num;
    fread(&num, sizeof(num), 1, file);
    printf("size of num: %lu\n", sizeof(num));
    printf("filename %s\n", argv[1]);
    while(!feof(file))
    {
        count ++;
        if (num != num)
        {
            printf("NaN detected at : %d\n", count);
        }
        fread(&num, sizeof(num), 1, file);
//        printf("count = %d, %f\n", count, num);
    }
    fclose(file);
}
