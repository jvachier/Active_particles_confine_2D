/*
 * Author: Jeremy Vachier
 * Purpose: ABP 2D confine in a square using an Euler-Mayurama algorithm
 * Language: C++
 * Date: 2023
 * Compilation line to use pragma: g++ name.cpp -fopenmp -o name (on mac run g++-12 ; 12 latest version obtain using brew list gcc)
 * Compilation line to use pragma, simd (vectorization) and tuple: g++ -O3 -std=c++17 name.cpp -fopenmp -o name
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

void update_position(double *x, double *orn, int Taj, double delta, double ratiol, double prefactor_va, double prefactor_xia, double inversetaua, double F, double xi_px, double xia, double prefactor_xi_px, default_random_engine &generator, normal_distribution<double> &Gaussdistribution)
{
#pragma omp parallel for simd num_threads(2)
	for (int k = 0; k < Taj; k++)
	{
		xia = Gaussdistribution(generator);
		xi_px = Gaussdistribution(generator);

		orn[k] = orn[k] - inversetaua * orn[k] * delta + prefactor_xia * xia;
		x[k] = x[k] + delta * (ratiol * sin(x[k] * ratiol) + F) + prefactor_va * orn[k] * delta + prefactor_xi_px * xi_px;
	}
}

void initialization(double *x, double *orn, int Taj)
{
#pragma omp parallel for simd num_threads(2)
	for (int k = 0; k < Taj; k++)
	{
		x[k] = 0.0;
		orn[k] = 0.0;
	}
}

tuple<double, double> moments(double *x, int Taj)
{
	double first = 0.0;
	double second = 0.0;
	for (int k = 0; k < Taj; k++)
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
	datacsv = fopen("./Da01_D01_01_09_2022.csv", "a");

	// check if the file parameter is exist
	if (parameter == NULL)
	{
		printf("no such file.");
		return 0;
	}

	// read the parameters from the file
	double F, delta, Da, taua;
	double ornn, xn;
	int Taj;
	fscanf(parameter, "%lf\t%lf\t%d\t%lf\t%lf\n", &F, &delta, &Taj, &Da, &taua);
	printf("%lf\t%lf\t%d\t%lf\t\t%lf\n", F, delta, Taj, Da, taua);

	double *x = (double *)malloc(Taj * sizeof(double));	  // x-position 1D
	double *orn = (double *)malloc(Taj * sizeof(double)); // activity

	// parameters
	const int N = 1E6; // number of iterations

	// initialization of the random generator
	random_device rdev;
	default_random_engine generator(rdev()); // random seed -> rdev
	// default_random_engine generator(1); //same seed

	// Distributions Gaussian
	normal_distribution<double> Gaussdistribution(0.0, 1.0);

	double D = 0.1;
	double xi_px; // noise for x-position
	double xi_py; // noise for y-position
	double xia;	  // noise activity
	double time = 0.0;
	double gamma = 1.0;
	double V0 = 1.0;
	double l0 = 2.0 * PI;
	double L = 2.0 * PI;
	double t0 = gamma * l0 * l0 / V0;
	double ratiol = l0 / L;
	int i, j, k;

	double inversetaua = t0 / taua;
	double prefactor_xia = sqrt(delta * inversetaua);
	double prefactor_va = sqrt(2.0 * Da * inversetaua);
	double prefactor_xi_px = sqrt(2.0 * D * delta);

	double Drift = 0.0;
	double Deff = 0.0;

	// safety check
	printf("%lf\t%lf\n", t0, prefactor_va / sqrt(t0));

	clock_t tStart = clock(); // check time for one trajectory

	// initialization position and activity
	initialization(x, orn, Taj);
	// fprintf(datacsv,"Force,Deff,Drfit,\n");	// only for the first run as using 'a' for file
	// Time evoultion
	for (j = 0; j < N; j++)
	{
		update_position(x, orn, Taj, delta, ratiol, prefactor_va, prefactor_xia, inversetaua, F, xi_px, xia, prefactor_xi_px, generator, Gaussdistribution);
	}
	// Compute the first and second moments
	auto [r_mean, rr_mean] = moments(x, Taj);

	rr_mean = rr_mean / Taj;
	r_mean = r_mean / Taj;
	Deff = (rr_mean - r_mean * r_mean) / (2.0 * (N - 1) * delta); // be careful when using 1D or 2D for the factor in the denominator
	Drift = r_mean / ((N - 1) * delta);
	fprintf(datacsv, "%lf,%lf,%lf\n", F, Deff, Drift);

	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC); // time for one trajectory

	free(x);
	free(orn);

	fclose(datacsv);
	return 0;
}
