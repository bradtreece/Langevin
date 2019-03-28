#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <algorithm>
#include "Random_Numbers.h"
#include "stats_and_histogram.h"

void histogram::init(double *observable, int samples, int bins)
{
    // bin_edges
    if (bins == 0)
    {
        bin_edges = floor(pow(4, log(samples) / log(20))) + 5;
    } else {
        bin_edges = bins + 1;
    }

    // constructor_samples
    constructor_samples = samples;

    // x
    x = new double[bin_edges];
    double smallest = *std::min_element(observable, observable + samples);
    double largest = *std::max_element(observable, observable + samples);
    if (smallest == largest) {largest = smallest + 0.001;}
    double bin_width = (largest - smallest)/(bin_edges - 5);
    for (int i = 0; i < bin_edges; i++)
    {
        x[i] = smallest + (i - 2)*bin_width;
    }

    // f
    f = new double[bin_edges];
    int counts, i, j;
    for (i = 0; i < bin_edges; i++)
    {
        counts = 0;
        for(j = 0; j < samples; j++)
        {
            if ((x[i] - 0.5*bin_width < observable[j]) and (observable[j] <= x[i] + 0.5*bin_width))
            {
                counts += 1;
            }
        }
        f[i] = counts/(bin_width*samples);
    }
    /*double sum = 0;
    for(i = 0; i < bin_edges - 1; i++)
    {
        sum += 0.5*(f[i] + f[i+1])*bin_width;
    }*/
    f_is_probability = true;
    
    return;
}

void histogram::extend_histogram_left(int numbins)
{
    int new_bin_edges = bin_edges + numbins;
    double temp_array[new_bin_edges] = {0.0};
    double bin_width = x[1] - x[0];

    for (int i=0; i<bin_edges; i++){temp_array[numbins + i] = f[i];}
    //delete [] f;
    f = new double[new_bin_edges];
    for (int i=0; i<new_bin_edges; i++){f[i] = temp_array[i];}

    for (int i=0; i<bin_edges; i++){temp_array[numbins + i] = x[i];}
    //delete [] x;
    x = new double[new_bin_edges];
    for (int i=0; i<new_bin_edges; i++){x[i] = temp_array[i];}
    for (int i=numbins - 1; i >= 0; i--){x[i] = x[i+1] - bin_width;}

    bin_edges = new_bin_edges;
    return;
}

void histogram::extend_histogram_right(int numbins)
{
    int new_bin_edges = bin_edges + numbins;
    double temp_array[new_bin_edges] = {0.0};
    double bin_width = x[1] - x[0];

    for (int i=0; i<bin_edges; i++){temp_array[i] = f[i];}
    //delete [] f;
    f = new double[new_bin_edges];
    for (int i=0; i<new_bin_edges; i++){f[i] = temp_array[i];}

    for (int i=0; i<bin_edges; i++){temp_array[i] = x[i];}
    //delete [] x;
    x = new double[new_bin_edges];
    for (int i=0; i<bin_edges;i++){x[i]=temp_array[i];}
    for (int i=bin_edges; i<new_bin_edges; i++){x[i] = x[i-1] + bin_width;}

    bin_edges = new_bin_edges;
    return;
}

void histogram::add_observations(double *observations, int num_of_observations)
{
    // Reconstruct Samples Per Bin
    double bin_width = x[1] - x[0];
    for (int i=0; i < bin_edges; i++)
    {
        f[i] = f[i]*constructor_samples*bin_width;
    }
    f_is_probability = false;

    // Update The Samples In The Bins
    double delta_x;
    int numbins, indx;
    for (int i=0; i < num_of_observations; i++)
    {
        if (bin_edges > 500)
        {
            /*consolidate_histogram(250, false);
            bin_width = x[1] - x[0];*/
            chop_histogram(250);
        }
        if ((x[0] - 0.5*bin_width <= observations[i]) 
        and (observations[i] <= x[bin_edges] + 0.5*bin_width))
        {
            indx = (int)(-1 + (observations[i] - x[0])/bin_width);
            if (indx < 0){indx = 0;}
            while (x[indx] - 0.5*bin_width <= observations[i])
            {
                indx += 1;
                if (indx > bin_edges)
                {
                    printf("PROBLEM CENTER\t");
                    return;
                }
            }
            f[indx-1] += 1;
        } else if (observations[i] < x[0] - 0.5*bin_width) 
        {
            delta_x = x[0] - 0.5*bin_width - observations[i];
            numbins = 3 + (delta_x / bin_width);
            extend_histogram_left(numbins);
            indx = (int)(-1 + (observations[i] - x[0])/bin_width);
            if (indx < 0){indx = 0;}
            while (x[indx] - 0.5*bin_width <= observations[i])
            {
                indx += 1;
                if (indx > bin_edges)
                {
                    printf("PROBLEM LEFT\t");
                    return;
                }
            }
            f[indx-1] += 1;
        } else if  (x[bin_edges] + 0.5*bin_width < observations[i])
        {
            delta_x = observations[i] - x[bin_edges] - 0.5*bin_width;
            numbins = 3 + (delta_x / bin_width);
            extend_histogram_right(numbins);
            indx = (int)(-1 + (observations[i] - x[0])/bin_width);
            if (indx < 0){indx = 0;}
            while (x[indx] - 0.5*bin_width <= observations[i])
            {
                indx += 1;
                if (indx > bin_edges)
                {
                    printf("PROBLEM RIGHT\t");
                    return;
                }
            }
            f[indx-1] += 1;
        } else
        {
            printf("\nSomething Went Wrong With Observation %i\n\n", i);
            return;
        }  
    }

    // Update the constructor_samples parameter in histogram
    constructor_samples += num_of_observations;

    // Reconstruct Frequency Per Bin
    for (int i=0; i < bin_edges; i++)
    {
        f[i] = f[i]/(constructor_samples*bin_width);
    }
    f_is_probability = true;
}

void histogram::consolidate_histogram(int reduction_amount, bool centered)
{
    int BIN_EDGES = bin_edges - reduction_amount;
    double BIN_WIDTH = (x[bin_edges-1] - x[0])/(BIN_EDGES - 2); // There will be one extra left and right of original
    double X[BIN_EDGES], F[BIN_EDGES];
    double bin_width = x[1] - x[0];
    int indx;

    if (centered)
    {
        X[0] = x[0] - BIN_WIDTH; // Center of bins aligned and add one to the left
    } else {
        X[0] = x[0] - 0.5*bin_width - 0.5*BIN_WIDTH; // Start the bins just slightly to the left of the original and add bin to the left
    }
    F[0] = 0.0;
    for (int i=1; i<BIN_EDGES; i++)
    {
        X[i] = X[i-1] + BIN_WIDTH;
        F[i] = 0.0;
    }

    if (!f_is_probability)
    {
        for (int i=0; i<bin_edges; i++)
        {f[i] = f[i]/bin_width;}
    }

    for (int i=0; i<bin_edges; i++)
    {
        indx = (int)(-1 + (x[i] - X[0])/BIN_WIDTH);
        if (indx < 0){indx = 0;}
        // While right side of new bin is less than the right side of the old bin, increment
        while (X[indx] + 0.5*BIN_WIDTH < x[i] + 0.5*bin_width) {indx+=1;}
        if (X[indx] - 0.5*BIN_WIDTH <= x[i] - 0.5*bin_width)
        {
            F[indx] += (f[i]*bin_width) / BIN_WIDTH;
        } else {
            F[indx] += (f[i]*( (x[i] + 0.5*bin_width) - (X[indx] - 0.5*BIN_WIDTH) )) / BIN_WIDTH;
            F[indx-1] += (f[i]*( (X[indx] - 0.5*BIN_WIDTH) - (x[i] - 0.5*bin_width) )) / BIN_WIDTH;
        }
    }

    if (!f_is_probability)
    {
        for (int i=0; i<BIN_EDGES; i++)
        {F[i] = F[i]*BIN_WIDTH;}
    }

    bin_edges = BIN_EDGES;
    x = new double[bin_edges];
    f = new double[bin_edges];
    for (int i=0; i<bin_edges; i++)
    {
        x[i] = X[i];
        f[i] = F[i];
    }

    return;
}

void histogram::chop_histogram(int leave_first_num)
{
    double sum1 = 0.0, sum2 = 0.0;
    double xnew[leave_first_num], fnew[leave_first_num];

    for (int i=0; i<bin_edges; i++)
    {
        sum1 += f[i];
        if (i < leave_first_num)
        {
            sum2 += f[i];
            xnew[i] = x[i];
            fnew[i] = f[i];
        }
    }

    x = new double[leave_first_num];
    f = new double[leave_first_num];

    for (int i=0; i < leave_first_num; i++)
    {
        x[i] = xnew[i];
        f[i] = fnew[i]*sum1/sum2;
    }

    bin_edges = leave_first_num;

    return;
}

double histogram::probability(double observation)
{
    if ((observation < x[0]) or (x[bin_edges-1] < observation))
    {return 0.0;}
    int i=0;
    while ( ( x[i] < observation ) and (i < bin_edges) )
    {i+=1;}
    if ((i == 0) or (i==bin_edges))
    {
        return 0.0;
    } else {
        double m = ( f[i] - f[i-1] )/( x[i] - x[i-1] );
        return m*( observation - x[i-1] ) + f[i-1];
    }
}