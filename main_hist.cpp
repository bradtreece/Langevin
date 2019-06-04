# include "brownian_particle.h"
# include "histogram.h"
# include <stdio.h>
# include <math.h>
# include <vector>
# include <time.h>

int main()
{

    srand(time(0));

    std::vector<int> indices;
    for (int i=100; i<1001; i+=100)
    {
        indices.push_back(i);
    }

    FILE *fp;
    fp = fopen("/home/btreece/Programs/BASIC_MPI/OUT.txt","w");

    
    brownian_particle b = brownian_particle(1001, 0.001, 1.0, 1.0, 1.0, 0.0, 0.0, 0.01);

    stat_list s_l = e_x | e_x2 | e_v | e_v2 | e_xv;
    brownian_histogram h = brownian_histogram(b, indices, s_l);
//    stat_list s_l = e_x2 | e_v2;
//    brownian_averages h = brownian_averages(b, s_l);

    h.add_samples(100000000, 1.0);
//    h.add_samples(1000, 1.0);
    h.print_data_to_file(fp);
    fclose(fp);

/*
    histogram h = histogram(0.0, 1.0, 10);
    printf("xmin\tb_w\tbins\tc_s\n");
    printf("%.3e\t%.3e\t%i\t%i\n\n", h.xmin, h.bin_width, h.bins, h.constructor_samples);

    h.increment_single_observation(1.05);
    printf("xmin\tb_w\tbins\tc_s\n");
    printf("%.3e\t%.3e\t%i\t%i\n\n", h.xmin, h.bin_width, h.bins, h.constructor_samples);
*/
    return 0;
}
