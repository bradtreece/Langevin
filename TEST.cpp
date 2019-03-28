# include "brownian_histogram.h"
# include <stdio.h>
# include <math.h>

int main()
{
    double param, fractpart, intpart;

    param = 3.14159265;
    fractpart = modf (param , &intpart);
    printf ("%f = %f + %f \n", param, intpart, fractpart);

    param = 4.73;
    fractpart = modf (param , &intpart);
    printf ("%f = %f + %f \n", param, intpart, fractpart);

    param = 2.50;
    fractpart = modf (param , &intpart);
    printf ("%f = %f + %f \n", param, intpart, fractpart);

    param = -1.5;
    fractpart = modf (param , &intpart);
    printf ("%f = %f + %f \n", param, intpart, fractpart);

    param = -0.2;
    fractpart = modf (param , &intpart);
    printf ("%f = %f + %f \n", param, intpart, fractpart);

    param = -2.8;
    fractpart = modf (param , &intpart);
    printf ("%f = %f + %f \n", param, intpart, fractpart);

    return 0;
}