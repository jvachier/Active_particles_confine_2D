#include <iostream>
#include <random>
#include <string>
#include <cmath>
#include <time.h>
#include <omp.h> //import library to use pragma
#include <tuple> //to output multiple components of a function
#define N_thread 4

using namespace std;

#include "initialization.h"

void initialization(
	double *x, double *y, int Particles,
	default_random_engine &generator, uniform_real_distribution<double> &distribution
)
{
#pragma omp parallel for simd num_threads(N_thread)
	for (int k = 0; k < Particles; k++)
	{
		x[k] = distribution(generator);
		y[k] = distribution(generator);
	}
}

void check_nooverlap(
	double *x, double *y, int Particles,
	double R, int L,
	default_random_engine &generator, uniform_real_distribution<double> &distribution
)
{
	int count = 0;
#pragma omp parallel for simd num_threads(N_thread)
	for (int k = 0; k < Particles; k++)
	{
		for (int j = 0; j < Particles; j++)
		{
			if (k != j)
			{
				R = sqrt((x[j]-x[k])*(x[j]-x[k]) + (y[j]-y[k])*(y[j]-y[k]));
				count = 0;
				while (R < 1.5 * L)
				{
					x[j] = distribution(generator);
					y[j] = distribution(generator);
					R = sqrt((x[j]-x[k])*(x[j]-x[k]) + (y[j]-y[k])*(y[j]-y[k]));
					count += 1;
					if (count > 3)
					{
						printf("Number of particle too high\n");
						exit(0);
					}
				}
			}
		}
	}

}