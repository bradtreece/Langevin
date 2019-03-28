#ifndef STATS_AND_HISTOGRAM_H
#define STATS_AND_HISTOGRAM_H
class histogram {
public:
    int bin_edges;
    int constructor_samples;
    double *x;
    double *f;
    bool f_is_probability;
    void init(double *observable, int samples, int bins = 0);
    void extend_histogram_left(int numbins);
    void extend_histogram_right(int numbins);
    void add_observations(double *observations, int num_of_observations);
    void consolidate_histogram(int reduction_amount, bool centered);
    void chop_histogram(int leave_first_num);
    double probability(double observation);
};
#endif