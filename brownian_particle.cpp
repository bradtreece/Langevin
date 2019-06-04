#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include "Random_Numbers.h"
#include "brownian_particle.h"
#include "histogram.h"

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
            atmp = ( -friction*vtmp + F*p[j] + force_ext->force_calculation(xtmp) ) / m;
            //xtmp += vtmp*dt + 0.5*atmp*dt*dt;
	        xtmp += vtmp*dt;
            vtmp += atmp*dt;
            ftmp = F*p[j];
        }
        x[i] = xtmp;
        v[i] = vtmp;
        f[i] = ftmp;
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

void brownian_histogram::add_samples(int number_of_samples, double variance_of_r)
{   
    // This is not the fastest way to do it
    std::vector<int>::iterator it_i;
    std::vector<histogram>::iterator it_h;
    for (int i = 0; i < number_of_samples; i++)
    {
        //if (i == 1) {printf("Made It\n");}
        b.calculate_trajectory(variance_of_r);

        if (!x_h.empty())
        {
            it_h = x_h.begin();
            for (it_i = indices.begin(); it_i != indices.end(); it_i++)
            {
                (*it_h).increment_single_observation(b.x[*it_i]);
                if ( (*it_h).bins > 1000 ) { (*it_h).triple_bin_width(); }
                it_h++;
            }
        }
        if (!x2_h.empty())
        {
            it_h = x2_h.begin();
            for (it_i = indices.begin(); it_i != indices.end(); it_i++)
            {
                (*it_h).increment_single_observation(b.x[*it_i]*b.x[*it_i]);
                if ( (*it_h).bins > 1000 ) { (*it_h).triple_bin_width(); }
                it_h++;
            }
        }
        if (!v_h.empty())
        {
            it_h = v_h.begin();
            for (it_i = indices.begin(); it_i != indices.end(); it_i++)
            {
                (*it_h).increment_single_observation(b.v[*it_i]);
                if ( (*it_h).bins > 1000 ) { (*it_h).triple_bin_width(); }
                it_h++;
            }
        }
        if (!v2_h.empty())
        {
            it_h = v2_h.begin();
            for (it_i = indices.begin(); it_i != indices.end(); it_i++)
            {
                (*it_h).increment_single_observation(b.v[*it_i]*b.v[*it_i]);
                if ( (*it_h).bins > 1000 ) { (*it_h).triple_bin_width(); }
                it_h++;
            }
        }
        if (!xv_h.empty())
        {
            it_h = xv_h.begin();
            for (it_i = indices.begin(); it_i != indices.end(); it_i++)
            {
                (*it_h).increment_single_observation(b.x[*it_i]*b.v[*it_i]);
                if ( (*it_h).bins > 1000 ) { (*it_h).triple_bin_width(); }
                it_h++;
            }
        }
    }
    return;
}

void brownian_histogram::print_data_to_file(FILE *fp)
{
    fprintf(fp, "NEW_RUN\tdt\tdT_out\tm\teta\tkT\tF\tD\n");
    fprintf(fp, "%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n", b.dt, b.dT_out, b.m, b.friction, b.kT(), b.F, b.D());

    std::vector<int>::iterator it_indices;
    std::vector<histogram>::iterator it_h;
    std::vector<int>::iterator it_f;

    if (!x_h.empty())
    {
        fprintf(fp, "x\n");

        it_h = x_h.begin();
        for (it_indices = indices.begin(); it_indices != indices.end(); it_indices++)
        {
            fprintf(fp, "t\t%.3e\tsamples\t%i\txmin\t%.3e\tbin_w\t%.3e\n"
              , (*it_indices)*(b.dT_out), (*it_h).constructor_samples, (*it_h).xmin, (*it_h).bin_width);
            for (it_f = (*it_h).f.begin(); it_f != (*it_h).f.end(); it_f++)
            { fprintf(fp, "%i\t", *it_f);}
            fprintf(fp,"\n");
            it_h++;
        }
    }
    if (!v_h.empty())
    {
        fprintf(fp, "v\n");

        it_h = v_h.begin();
        for (it_indices = indices.begin(); it_indices != indices.end(); it_indices++)
        {
            fprintf(fp, "t\t%.3e\tsamples\t%i\txmin\t%.3e\tbin_w\t%.3e\n\t"
              , (*it_indices)*(b.dT_out), (*it_h).constructor_samples, (*it_h).xmin, (*it_h).bin_width);
            for (it_f = (*it_h).f.begin(); it_f != (*it_h).f.end(); it_f++)
            { fprintf(fp, "%i\t", *it_f);}
            fprintf(fp,"\n");
            it_h++;
        }
    }
    if (!x2_h.empty())
    {
        fprintf(fp, "x2\n");

        it_h = x2_h.begin();
        for (it_indices = indices.begin(); it_indices != indices.end(); it_indices++)
        {
            fprintf(fp, "t\t%.3e\tsamples\t%i\txmin\t%.3e\tbin_w\t%.3e\n\t"
              , (*it_indices)*(b.dT_out), (*it_h).constructor_samples, (*it_h).xmin, (*it_h).bin_width);
            for (it_f = (*it_h).f.begin(); it_f != (*it_h).f.end(); it_f++)
            { fprintf(fp, "%i\t", *it_f);}
            fprintf(fp,"\n");
            it_h++;
        }
}
    if (!v2_h.empty())
    {
        fprintf(fp, "v2\n");

        it_h = v2_h.begin();
        for (it_indices = indices.begin(); it_indices != indices.end(); it_indices++)
        {
            fprintf(fp, "t\t%.3e\tsamples\t%i\txmin\t%.3e\tbin_w\t%.3e\n\t"
              , (*it_indices)*(b.dT_out), (*it_h).constructor_samples, (*it_h).xmin, (*it_h).bin_width);
            for (it_f = (*it_h).f.begin(); it_f != (*it_h).f.end(); it_f++)
            { fprintf(fp, "%i\t", *it_f);}
            fprintf(fp,"\n");
            it_h++;
        }
}
    if (!xv_h.empty())
    {
        fprintf(fp, "xv\n");

        it_h = xv_h.begin();
        for (it_indices = indices.begin(); it_indices != indices.end(); it_indices++)
        {
            fprintf(fp, "t\t%.3e\tsamples\t%i\txmin\t%.3e\tbin_w\t%.3e\n\t"
              , (*it_indices)*(b.dT_out), (*it_h).constructor_samples, (*it_h).xmin, (*it_h).bin_width);
            for (it_f = (*it_h).f.begin(); it_f != (*it_h).f.end(); it_f++)
            { fprintf(fp, "%i\t", *it_f);}
            fprintf(fp,"\n");
            it_h++;
        }
}
    return;
}