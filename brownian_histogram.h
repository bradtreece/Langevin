#ifndef BROWNIAN_HISTOGRAM_H
#define BROWNIAN_HISTOGRAM_H
#include <vector>
class brownian_histogram {
public:
    int bin_edges;
    int constructor_samples;
    double bin_width;
    std::vector<double> x;
    std::vector<int> f;
    brownian_histogram(double min_x, double max_x, int num_bins, bool add_padding = false)
    : constructor_samples(0)
    , bin_edges(num_bins+1)
    , bin_width( (max_x - min_x)/num_bins )
    {
//        double bin_width = ( max_x - min_x ) / num_bins;
        double x0 = min_x;
        if (add_padding)
        {
            // Add 4 bins (2 left, 2 right)
            bin_edges += 4;
            x0 -= 2*bin_width;
        }
        
        x.assign(bin_edges, 0.0);
        f.assign(bin_edges, 0.0);

        for (int i=0; i<bin_edges; i++)
        { x.at(i) = x0 + i*bin_width; }
    }

    ~brownian_histogram(){}
    
    void increment_single_observation(double observation);
    // void init(double *observable, int samples, int bins = 0);
    // void extend_histogram_left(int numbins);
    // void extend_histogram_right(int numbins);
    // void add_observations(double *observations, int num_of_observations);
    // void consolidate_histogram(int reduction_amount, bool centered);
    // void chop_histogram(int leave_first_num);
    // double probability(double observation);
};
#endif