#include "stats_and_histogram.h"
#include <math.h>
#include <stdio.h>

enum stat_list {
    e_x  = (1 << 0),
    e_x2 = (1 << 1),
    e_v  = (1 << 2),
    e_v2 = (1 << 3),
    e_xv = (1 << 4)
};
inline stat_list operator|(stat_list a, stat_list b)
{return static_cast<stat_list>(static_cast<int>(a) | static_cast<int>(b));};

class brownian_particle {
public:
    int num_steps;
    double dt;
    double friction;
    double F;
    double m;
    double dT_out;
    double *x;
    double *v;
    double *f;

    brownian_particle(int num_steps, double dt, double friction,
    double kT, double m, double x0, double v0, double dT_output)
    : num_steps{num_steps}
    , dt{dt}
    , friction{friction}
    , F{pow(2.0*kT*friction/dT_output, 0.5)}
    , m{m}
    , dT_out{dT_output}
    {
        x = new double[num_steps];
        v = new double[num_steps];
        f = new double[num_steps];
        x[0] = x0;
        v[0] = v0;
        f[0] = 0.0;

        // Check if dT_out is a integer multiple of dt
        double res = remainder(dT_output, dt);
        if (not ( fabs(res) < 0.001*dt) )
        {
            printf("\ndT_out was not an integer multiple of dt.\n");
            dT_out = round(dT_output/dt)*dt;
            F = pow(2.0*kT*friction/dT_out, 0.5);
        }
    }

    brownian_particle(int num_steps, double dt, double friction,
    double kT, double m, double x0, double v0)
    : brownian_particle(num_steps, dt, friction, kT, m, x0, v0, dt) {}

    brownian_particle(int num_steps, double dt, double friction,
    double kT, double m)
    : brownian_particle(num_steps, dt, friction, kT, m, 0.0, 0.0, dt) {}

    void calculate_trajectory(double variance_of_r = NAN);
    void repeat_trajectory_to_init_xsquared_histogram(int number_of_repeats, histogram *h);
    void repeat_trajectory_to_update_xsquared_histogram(int number_of_repeats, histogram *h);
    double kT(){return 0.5*F*F*dT_out/friction;}
    ~brownian_particle(){};
};



class brownian_averages {
private:
    brownian_particle b;
public:
    int num_samples;
    double *xbar;
    double *x2bar;
    double *vbar;
    double *v2bar;
    double *xvbar;

    brownian_averages(brownian_particle b, stat_list s)
    : b{b}
    , num_samples(0)
    {
        xbar  = nullptr;
        x2bar = nullptr;
        vbar  = nullptr;
        v2bar = nullptr;
        xvbar = nullptr;

        if (e_x  & s)
        { xbar  = new double[b.num_steps]; }
        if (e_x2 & s)
        { x2bar = new double[b.num_steps]; }
        if (e_v  & s)
        { vbar  = new double[b.num_steps]; }
        if (e_v2 & s)
        { v2bar = new double[b.num_steps]; } 
        if (e_xv & s)
        { xvbar = new double[b.num_steps]; } 
    }

    int num_steps()
    { return b.num_steps; };

    void add_samples(int number_of_samples, double variance_of_r);

    void print_data_to_file(FILE *fp);

    ~brownian_averages(){};//{printf("BYE!\n");}
};