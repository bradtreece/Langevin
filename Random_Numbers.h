#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Generate A Random Number Sampled From U[0,1]
double rand_uniform_0_1();

// Generate A Random Number Sampled From N(mu, sigma) Given "p" Sampled From U[0,1]
double stdnormal_inv(double p, double mu=0.0, double sigma=1.0);

// Generate An Array Of Numbers Sampled From N(mu, sigma)
void Rand_Normal_Array(double_t *arr, int length, double *mu = nullptr, double *sigma = nullptr);
