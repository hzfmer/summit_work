#include <stdio.h>
#include <math.h>

int main(void)
{
    float tmax = 30.001;
    float dt = 0.001;
    int nt;
    nt = (int) floorf(tmax / dt);
    printf("%d", nt);
    return 0;
}
