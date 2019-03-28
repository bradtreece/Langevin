#include "brownian_histogram.h"
#include "math.h"

void brownian_histogram::increment_single_observation(double observation)
{
    if ( observation < x.front() )
    { /* extend left */ } else if ( observation > x.back() )
    { /* extend rhgt */ } else
    {
        /* add observation */
    }
}