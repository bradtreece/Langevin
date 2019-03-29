#include "brownian_histogram.h"
#include <math.h>
#include <stdlib.h>

int brownian_histogram::get_index_of_observation(double observation)
{
    return (int) floor( ( observation - xmin ) / bin_width );
}

void brownian_histogram::extend_histogram_left(int num_bins_to_add)
{
    f.insert(f.begin(), num_bins_to_add, 0);
    bins += num_bins_to_add;
    xmin -= num_bins_to_add*bin_width;
    return;
}

void brownian_histogram::extend_histogram_right(int num_bins_to_add)
{
    bins += num_bins_to_add;
    f.resize(bins, 0);
    return;
}

void brownian_histogram::increment_single_observation(double observation)
{
    int indx;
    bool index_fails = true;
    
    while (index_fails)
    {
        indx = get_index_of_observation(observation);

        if ( indx < 0 )
        { /* extend left */
            extend_histogram_left(2*abs(indx));
        } else if ( indx >= bins )
        { /* extend right */
            extend_histogram_left( 2*(indx - bins + 1) ); 
        } else
        {
            index_fails = false;
        }
    }

    f.at(indx)++;

    return;
}