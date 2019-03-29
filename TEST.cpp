# include "brownian_histogram.h"
# include <stdio.h>
# include <math.h>

int main()
{
    double param;

    param = 3.14159265;
    printf ("%f = %i \n", param, (int) floor( param ) );

    param = 4.73;
    printf ("%f = %i \n", param, (int) floor( param ) );

    param = 2.50;
    printf ("%f = %i \n", param, (int) floor( param ) );

    param = -1.5;
    printf ("%f = %i \n", param, (int) floor( param ) );

    param = -0.2;
    printf ("%f = %i \n", param, (int) floor( param ) );

    param = -2.8;
    printf ("%f = %i \n", param, (int) floor( param ) );

    return 0;
}