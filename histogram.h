#include <vector>

#ifndef HISTOGRAM_H
#define HISTOGRAM_H
class histogram {
public:
    int bins;
    int constructor_samples;
    double bin_width;
    double xmin;
    std::vector<int> f;
    histogram(double min_x, double max_x, int num_bins, bool add_padding = false)
    : bins(num_bins)
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

    ~histogram(){}
    
    std::vector<int>::iterator index_of_median();
    int get_index_of_observation(double observation);
    void extend_histogram_left(int num_bins_to_add);
    void extend_histogram_right(int num_bins_to_add);
    void increment_single_observation(double observation);
    void double_bin_width();
    void triple_bin_width();
    double probability(double observation);
};
#endif