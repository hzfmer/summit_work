void inverse(double *x, int xlen)
{
    int i;
    double temp;
    if (xlen % 2 == 1)
    {
        for (i=0; i<(xlen - 1)/2; i++)
        {
            temp = *(x + 1);
            *(x+1) = *(x + xlen - 1 - i);
            *(x + xlen - 1 - i) = temp;
        }
    }
    else
    {
        for (i=0; i < xlen / 2; i++)
        {
            temp = *(x + i);
            *(x + i) = *(x + xlen - 1 - i);
            *(x + xlen - 1 - i) = temp;
        }
    }
}

