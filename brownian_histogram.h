#include <vector>
#ifndef BROWNIAN_HISTOGRAM_H
#define BROWNIAN_HISTOGRAM_H
class brownian_histogram {
public:
    int bins;
    int constructor_samples;
    double bin_width;
    double xmin;
    std::vector<int> f;
    brownian_histogram(double min_x, double max_x, int num_bins, bool add_padding = false)
    : bins(num_bins+1)
    , constructor_samples(0)
    , bin_width( (max_x - min_x)/num_bins )
    , xmin(min_x)
    {
//        double bin_width = ( max_x - min_x ) / num_bins;
        if (add_padding)
        {
            // Add 4 bins (2 left, 2 right)
            bins += 4;
            xmin -= 2*bin_width;
        }
        f.assign(bins, 0);
    }

    ~brownian_histogram(){}
    
    int get_index_of_observation(double observation);
    void extend_histogram_left(int num_bins_to_add);
    void extend_histogram_right(int num_bins_to_add);
    void increment_single_observation(double observation);
    void double_bin_width();
    // void init(double *observable, int samples, int bins = 0);
    // void extend_histogram_left(int numbins);
    // void extend_histogram_right(int numbins);
    // void add_observations(double *observations, int num_of_observations);
    // void consolidate_histogram(int reduction_amount, bool centered);
    // void chop_histogram(int leave_first_num);
    // double probability(double observation);
};
#endif