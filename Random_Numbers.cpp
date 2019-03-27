#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Random_Numbers.h"

double rand_uniform_0_1()
{
    double p = rand() / nextafter((double) RAND_MAX, __DBL_MAX__);

    return p;
}

double stdnormal_inv(double p, double mu, double sigma)
{
    double a[6];
    double b[5];
    double c[6];
    double d[4];

    a[0] = -3.969683028665376E+01;
    a[1] =  2.209460984245205E+02;
    a[2] = -2.759285104469687E+02;
    a[3] =  1.383577518672690E+02;
    a[4] = -3.066479806614716E+01;
    a[5] =  2.506628277459239E+00;

    b[0] = -5.447609879822406E+01;
    b[1] =  1.615858368580409E+02;
    b[2] = -1.556989798598866E+02;
    b[3] =  6.680131188771972E+01;
    b[4] = -1.328068155288572E+01;

    c[0] = -7.784894002430293E-03;
    c[1] = -3.223964580411365E-01;
    c[2] = -2.400758277161838E+00;
    c[3] = -2.549732539343734E+00;
    c[4] =  4.374664141464968E+00;
    c[5] =  2.938163982698783E+00;

    d[0] =  7.784695709041462E-03;
    d[1] =  3.224671290700398E-01;
    d[2] =  2.445134137142996E+00;
    d[3] =  3.754408661907416E+00;

    double q, t, u;

    if (isnan(p) || p >= 1.0 || p <= 0.0)
    { return nan(""); }

    q = fmin(p, 1-p);
    if (q < 0.02425){
        // Rational approximation for central region
        u = q - 0.5;
        t = u*u;
        u = u*(((((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4])*t+a[5])
        /(((((b[0]*t+b[1])*t+b[2])*t+b[3])*t+b[4])*t+1);
    } else {
        // Rational approximation for tail region.
        t = sqrt(-2*log(q));
        u = (((((c[0]*t+c[1])*t+c[2])*t+c[3])*t+c[4])*t+c[5])
        /((((d[0]*t+d[1])*t+d[2])*t+d[3])*t+1);
    }
    // For refinement to full precision accuracy, review the following webpage:
    // https://web.archive.org/web/20150910091202/http://home.online.no:80/~pjacklam/notes/invnorm/impl/lea/lea.c
    
    if (p < 0.5)
    { return -sigma*u + mu;}
    else
    { return sigma*u + mu; }
    
}

void Rand_Normal_Array(double_t *arr, int length, double *mu, double *sigma)
{
    double MU[length], SIGMA[length], p;
    for (int i=0; i<length; i++)
    {
        MU[i] = 0.0;
        SIGMA[i] = 1.0;
    }

    if (mu != nullptr)
    {
        for (int i=0; i < length; i++)
        {MU[i] = mu[i];}
    }

    if (sigma != nullptr)
    {
        for (int i=0; i < length; i++)
        {SIGMA[i] = sigma[i];}
    }

    p = NAN;
    for (int i=0;i<length;i++){
        while (isnan(p))
        {
            p = stdnormal_inv(rand_uniform_0_1(), MU[i], SIGMA[i]);
        }
        arr[i] = p;
    }
    return;
}
