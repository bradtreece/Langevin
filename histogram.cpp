#include "histogram.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

std::vector<int>::iterator histogram::index_of_median()
{
    std::vector<int>::iterator it = f.begin();
    int samples = *it;
    while (samples < constructor_samples/2)
    {
        it++;
        samples += *it;
    }
    return it;
}

int histogram::get_index_of_observation(double observation)
{
    return (int) floor( ( observation - xmin ) / bin_width );
}

void histogram::extend_histogram_left(int num_bins_to_add)
{
    f.insert(f.begin(), num_bins_to_add, 0);
    bins += num_bins_to_add;
    xmin -= num_bins_to_add*bin_width;
    return;
}

void histogram::extend_histogram_right(int num_bins_to_add)
{
    bins += num_bins_to_add;
    f.resize(bins, 0);
    return;
}

void histogram::increment_single_observation(double observation)
{
    int indx;
    bool index_fails = true;
    
    while (index_fails)
    {
        indx = get_index_of_observation(observation);

        if ( indx < 0 )
        { /* extend left */
            extend_histogram_left(abs(indx));
        } else if ( indx >= bins )
        { /* extend right */
            extend_histogram_right( (indx - bins + 1) ); 
        } else
        {
            index_fails = false;
        }
    }

    f.at(indx)++;
    constructor_samples++;

    return;
}

void histogram::double_bin_width()
{
    if (bins % 2 != 0)
    {
        extend_histogram_right(1);
    }

    bins = bins/2;
    bin_width = 2.0*bin_width;

    std::vector<int> f_temp(bins, 0);
    std::vector<int>::iterator it_f = f.begin();
    std::vector<int>::iterator it_ft = f_temp.begin();

    for (int i=0; i<bins; i++)
    {
        // Put f_{2i} into f_temp_{i}
        *it_ft = *it_f;
        it_f++;
        // Add f_{2i+1} to f_temp_{i}
        *it_ft += *it_f;
        it_f++;
        it_ft++;
    }
    
    f.swap(f_temp);
    return;
}

void histogram::triple_bin_width()
{
    std::vector<int>::iterator it = index_of_median();
    //printf("\n\nMedian Bin: %li, Median x: %f\n\n", it - f.begin(), (it - f.begin() + 0.5)*bin_width + xmin);
    // We need one on either side of the median to make three for the central bin
    // and the number of bins left or right should be a multiple of three.
    // Hence, bins to the left and right need to be one more than a multiple of 3. 
    int L = 3 - ( (it - f.begin() - 1) % 3 );
    int R = 3 - ( (f.end() - it - 2) % 3 ); // f.end is one past the end

    if (L != 3)
    { extend_histogram_left(L); }
    if (R != 3)
    { extend_histogram_right(R); }

    if (bins % 3 != 0)
    {
        printf("SOMETHING WENT WRONG WITH TRIPLING THE BIN WIDTH\n");
        return;
    }

    bins = bins/3;
    bin_width = 3.0*bin_width;

    std::vector<int> f_temp(bins, 0);
    std::vector<int>::iterator it_f = f.begin();
    std::vector<int>::iterator it_ft = f_temp.begin();

    for (int i=0; i<bins; i++)
    {
        // Put f_{3i} into f_temp_{i}
        *it_ft = *it_f;
        it_f++;
        // Add f_{3i+1} to f_temp_{i}
        *it_ft += *it_f;
        it_f++;
        // Add f_{3i+2} to f_temp_{i}
        *it_ft += *it_f;
        it_f++;
        it_ft++;
    }
    
    f.swap(f_temp);
    //it = index_of_median();
    //printf("\n\nMedian Bin: %li, Median x: %f\n\n", it - f.begin(), (it - f.begin() + 0.5)*bin_width + xmin);
    return;
}

double histogram::probability(double observation)
{
    // The actual x positions for f_i are taken to be the middle of the bins.
    // Therefore, the histogram extends half a bin width in either direction,
    // giving a smooth decrease to zero at these points.

    // Figure out which bin. Again, the f_i are at the centers of bins.
    // This means the observation needs shifted left half a bin to 
    // correspond to the correct f bin (different from x bin).
    // i.e. - Anything between x_0 + 0.5 and x_1 + 0.5 falls
    //        between f_0 and f_1.
    double obs_shifted = observation - 0.5*bin_width;

    // Outside the bonds of the histogram, this happens (after shifting) when
    // obs is in (xmin - bin_width, xmax) since there is an extra f point left
    // and right of the histogram.
    // These extra f are zero by definition above.
    if ( (obs_shifted <= xmin - bin_width) or (obs_shifted >= xmin + bins*bin_width) )
    { return 0.0; }

    int indx = get_index_of_observation(obs_shifted);

    // Edge cases
    if (indx == -1)
    { /* left edge: f_{-1} = 0 */
        return  f.at(0) * (obs_shifted - xmin - indx*bin_width) / bin_width;
    } else if (indx == bins-1)
    { /* right edge: f_{bins} = 0 */
        return -f.at(indx) * (obs_shifted - xmin - indx*bin_width) / bin_width;
    } else 
    {
        return (f.at(indx+1) - f.at(indx)) * (obs_shifted - xmin - indx*bin_width) / bin_width;
    }
}