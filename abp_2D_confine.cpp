/*
 * Author: Jeremy Vachier
 * Purpose: ABP 2D confine in a square using an Euler-Mayurama algorithm
 * Language: C++
 * Date: 2023
 * Compilation line to use pragma: g++ name.cpp -fopenmp -o name.o (on mac run g++-12 ; 12 latest version obtain using brew list gcc)
 * Compilation line to use pragma, simd (vectorization) and tuple: g++ -O3 -std=c++17 name.cpp -fopenmp -o name.o
 */

#include <iostream>
#include <random>
#include <string>
#include <cmath>
#include <time.h>
#include <omp.h> //import library to use pragma
#include <tuple> //to output multiple components of a function

#define PI 3.141592653589793

using namespace std;

void update_position(
	double *x, double *y, double phi, double prefactor_e, int Particles, 
	double delta, double De, double Dt, double xi_e, double xi_px, 
	double xi_py, double vs, double prefactor_xi_px, double prefactor_xi_py,
	default_random_engine &generator, normal_distribution<double> &Gaussdistribution
)
{
#pragma omp parallel for simd num_threads(2)
	for (int k = 0; k < Particles; k++)
	{
		xi_e = Gaussdistribution(generator);
		xi_px = Gaussdistribution(generator);
		xi_py = Gaussdistribution(generator);

		phi = phi + prefactor_e * xi_e;
		x[k] = x[k] + vs * cos(phi) * delta + xi_px * prefactor_xi_px;
		y[k] = y[k] + vs * sin(phi) * delta + xi_py * prefactor_xi_py;
	}
}

void print_file(
	double *x, double *y,
	int Particles, int time,
	FILE *datacsv
)
{
	for(int k = 0; k < Particles; k++)
	{
		fprintf(datacsv,"Particles%d,%lf,%lf,%d,\n",k,x[k],y[k],time);
	}
}

void initialization(
	double *x, double *y, int Particles,
	double x_x, double y_y,
	default_random_engine &generator, uniform_real_distribution<double> &distribution
)
{
#pragma omp parallel for simd num_threads(2)
	for (int k = 0; k < Particles; k++)
	{
		x_x = distribution(generator);
		y_y = distribution(generator);
		x[k] = x_x;
		y[k] = y_y;
	}
}

tuple<double, double> moments(double *x, int Particles)
{
	double first = 0.0;
	double second = 0.0;
	for (int k = 0; k < Particles; k++)
	{
		second += x[k] * x[k];
		first += x[k];
	}
	return {first, second};
}

int main(int argc, char *argv[])
{
	// File
	FILE *datacsv;
	FILE *parameter;
	parameter = fopen("parameter.txt", "r");
	datacsv = fopen("./simulation_test.csv", "w");

	// check if the file parameter is exist
	if (parameter == NULL)
	{
		printf("no such file.");
		return 0;
	}

	// read the parameters from the file
	double epsilon, delta, Dt, De, vs;
	int Particles;
	fscanf(parameter, "%lf\t%lf\t%d\t%lf\t%lf\t%lf\n", &epsilon, &delta, &Particles, &Dt, &De, &vs);
	printf("%lf\t%lf\t%d\t%lf\t%lf\t%lf\n", epsilon, delta, Particles, Dt, De, vs);

	double *x = (double *)malloc(Particles * sizeof(double)); // x-position 
	double *y = (double *)malloc(Particles * sizeof(double)); // y-position

	// parameters
	const int N = 1E2; // number of iterations
	const int L = 1.0; // particle size

	// initialization of the random generator
	random_device rdev;
	default_random_engine generator(rdev()); // random seed -> rdev
	// default_random_engine generator(1); //same seed

	// Distributions Gaussian
	normal_distribution<double> Gaussdistribution(0.0, 1.0);
	// Distribution Uniform for initialization
	uniform_real_distribution<double> distribution(-10.0,10.0);

	double xi_px; // noise for x-position
	double xi_py; // noise for y-position
	double xi_e; // noise ortientation
	double x_x; // used to initialize
	double y_y; // used to initialize
	int i, j, k;

	double phi = 0.0;
	double prefactor_e = sqrt(2.0 * delta * De);
	double prefactor_xi_px = sqrt(2.0 * delta * Dt);
	double prefactor_xi_py = sqrt(2.0 * delta * Dt);


	clock_t tStart = clock(); // check time for one trajectory

	// initialization position and activity
	initialization(
		x, y, Particles,
		x_x, y_y,
		generator, distribution
	);
	// fprintf(datacsv,"Force,Deff,Drfit,\n");	// only for the first run as using 'a' for file
	// Time evoultion
	int time;
	for (time = 0; time < N; time++)
	{
		update_position(
			x, y, phi, prefactor_e, Particles, 
			delta, De, Dt, xi_e, xi_px, 
			xi_py, vs, prefactor_xi_px, prefactor_xi_py,
			generator, Gaussdistribution
		);
		if(time % 10 == 0 && time >= 0)
		{
			print_file(
				x, y,
				Particles, time,
				datacsv
			);
		}
	}
	// Compute the first and second moments
	//auto [r_mean, rr_mean] = moments(x, Taj);


	

	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC); // time for one trajectory

	free(x);
	free(y);

	fclose(datacsv);
	return 0;
}
