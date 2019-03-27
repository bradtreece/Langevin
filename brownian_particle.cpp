#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include "Random_Numbers.h"
#include "brownian_particle.h"
#include "stats_and_histogram.h"

void brownian_particle::calculate_trajectory(double variance_of_r)
{
    int tau = (int) round(dT_out/dt);
    double p[tau], sigma[tau], xtmp, vtmp, ftmp, atmp;

    if (isnan(variance_of_r))
    {
        for (int i = 0; i < tau; i++) {sigma[i]=1.0;}
    } else {
        for (int i = 0; i < tau; i++) {sigma[i]=pow(variance_of_r, 0.5);}
    }

    for (int i=1; i < num_steps; i++)
    {
        xtmp = x[i-1];
        vtmp = v[i-1];
        ftmp = 0.0;
        Rand_Normal_Array(p, tau, nullptr, sigma);
        for (int j=0; j<tau; j++)
        {
            atmp = ( -friction*vtmp + F*p[j] ) / m;
            xtmp += vtmp*dt + 0.5*atmp*dt*dt;
            vtmp += atmp*dt;
            ftmp = F*p[j];
        }
        x[i] = xtmp;
        v[i] = vtmp;
        f[i] = ftmp;
    }
    return;
}

void brownian_particle::repeat_trajectory_to_init_xsquared_histogram(int number_of_repeats, histogram *h)
{
    double set_of_x2m[number_of_repeats][num_steps];
    double observable[number_of_repeats];


    for (int i=0; i < number_of_repeats; i++)
    {
        calculate_trajectory();
        mean_square_displacement_from_reference(x, set_of_x2m[i], num_steps, x[0]);
    }
    for (int i=0; i<num_steps; i++)
    {
        for (int j=0; j < number_of_repeats; j++)
        {observable[j]=set_of_x2m[j][i];}


        h[i].init(observable, number_of_repeats);
    }
    return;
}

void brownian_particle::repeat_trajectory_to_update_xsquared_histogram(int number_of_repeats, histogram *h)
{
    double set_of_x2m[number_of_repeats][num_steps];
    double observable[number_of_repeats];
    int i;

    for (i=0; i < number_of_repeats; i++)
    {
        calculate_trajectory();
        mean_square_displacement_from_reference(x, set_of_x2m[i], num_steps, x[0]);
    }

    for (i=0; i<num_steps; i++)
    {
        for (int j=0; j < number_of_repeats; j++)
        {observable[j]=set_of_x2m[j][i];}

        h[i].add_observations(observable, number_of_repeats);

    }

    return;
}


void brownian_averages::add_samples(int number_of_samples, double variance_of_r)
{
    for (int i = 0; i < number_of_samples; i++)
    {
        b.calculate_trajectory(variance_of_r);
        num_samples += 1;
        for (int j = 0; j < b.num_steps; j++)
        {
            if (xbar)
            { xbar[j]  = ( xbar[j]*(num_samples-1)  + b.x[j] ) / num_samples; }
            if (x2bar)
            { x2bar[j] = ( x2bar[j]*(num_samples-1) + b.x[j]*b.x[j] ) / num_samples; }
            if (vbar)
            { vbar[j]  = ( vbar[j]*(num_samples-1)  + b.v[j] ) / num_samples; }
            if (v2bar)
            { v2bar[j] = ( v2bar[j]*(num_samples-1) + b.v[j]*b.v[j] ) / num_samples; }
            if (xvbar)
            { xvbar[j] = ( xvbar[j]*(num_samples-1) + b.x[j]*b.v[j] ) / num_samples; }
	}
    }
    return;
}

void brownian_averages::print_data_to_file(FILE *fp)
{
    fprintf(fp, "NEW_RUN\tdt\tdT_out\tm\teta\tkT\tF\tn_samples\n");
    fprintf(fp, "%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%i\n", b.dt, b.dT_out, b.m, b.friction, b.kT(), b.F, num_samples);

    if (xbar)
    {
        fprintf(fp, "<x>:");
        for (int i = 0; i < b.num_steps; i++)
        { fprintf(fp, "\t%.4e", xbar[i]);}
        fprintf(fp,"\n");
    }
    if (vbar)
    {
        fprintf(fp, "<v>:");
        for (int i = 0; i < b.num_steps; i++)
        { fprintf(fp, "\t%.4e", vbar[i]);}
        fprintf(fp,"\n");
    }
    if (x2bar)
    {
        fprintf(fp, "<xx>:");
        for (int i = 0; i < b.num_steps; i++)
        { fprintf(fp, "\t%.4e", x2bar[i]);}
        fprintf(fp,"\n");
    }
    if (v2bar)
    {
        fprintf(fp, "<vv>:");
        for (int i = 0; i < b.num_steps; i++)
        { fprintf(fp, "\t%.4e", v2bar[i]);}
        fprintf(fp,"\n");
    }
    if (xvbar)
    {
        fprintf(fp, "<xv>:");
        for (int i = 0; i < b.num_steps; i++)
        { fprintf(fp, "\t%.4e", xvbar[i]);}
        fprintf(fp,"\n");
    }
    return;
}
