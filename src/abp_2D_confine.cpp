/*
 * Author: Jeremy Vachier
 * Purpose: ABP 2D confine in a square using an Euler-Mayurama algorithm
 * Language: C++
 * Date: 2023
 */

#include <iostream>
#include <random>
#include <cstring>
#include <stdio.h>
#include <cmath>
#include <time.h>
#include <omp.h> //import library to use pragma
#include <tuple> //to output multiple components of a function

#include "headers/print_file.h"
#include "headers/reflective_boundary_conditions.h"
#include "headers/circular_reflective_boundary_conditions.h"
#include "headers/initialization.h"
#include "headers/update_position.h"
#include "headers/check_nooverlap.h"

#define PI 3.141592653589793
#define N_thread 6

using namespace std;

int main(int argc, char *argv[])
{
	// File
	FILE *datacsv;
	FILE *parameter;
	parameter = fopen("parameter.txt", "r");
	datacsv = fopen("./data/simulation.csv", "w");

	// check if the file parameter is exist
	if (parameter == NULL)
	{
		printf("no such file.");
		return 0;
	}

	omp_set_num_threads(N_thread);

	// read the parameters from the file
	double epsilon, delta, Dt, De, vs;
	double Wall;
	int Particles;
	int N; // number of iterations
	char name[100];
	char key1[] = "circular";
	char key2[] = "squared";

	bool flag = false;

	printf("Select confinement geometry, either squared or circular:");
	scanf("%s", name);

	while (flag == false)
	{
		if ((strcmp(name, key1) == 0) or (strcmp(name, key2) == 0))
		{
			flag = true;
		}
		else
		{
			printf("You have not selected the correct, please select again\n");
			printf("Select confinement geometry, either squared or circular:");
			scanf("%s", name);
			flag = false;
		}
	}
	fscanf(parameter, "%lf\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%d\n", &epsilon, &delta, &Particles, &Dt, &De, &vs, &Wall, &N);
	printf("%lf\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%d\n", epsilon, delta, Particles, Dt, De, vs, Wall, N);

	double *x = (double *)malloc(Particles * sizeof(double)); // x-position
	double *y = (double *)malloc(Particles * sizeof(double)); // y-position

	// parameters
	const int L = 1.0; // particle size

	// initialization of the random generator
	random_device rdev;
	default_random_engine generator(rdev()); // random seed -> rdev
	// default_random_engine generator(1); //same seed

	// Distributions Gaussian
	normal_distribution<double> Gaussdistribution(0.0, 1.0);
	// Distribution Uniform for initialization
	uniform_real_distribution<double> distribution(-Wall, Wall);
	// uniform_real_distribution<double> distribution_e(0.0,360.0*PI / 180.0); // directly in radian
	uniform_real_distribution<double> distribution_e(0.0, 360.0);

	double xi_px = 0.0; // noise for x-position
	double xi_py = 0.0; // noise for y-position
	double xi_e = 0.0;  // noise ortientation

	double phi = 0.0;
	double prefactor_e = sqrt(2.0 * delta * De);
	double prefactor_xi_px = sqrt(2.0 * delta * Dt);
	double prefactor_xi_py = sqrt(2.0 * delta * Dt);
	double prefactor_interaction = epsilon * 48.0;
	double r = 5.0 * L;

	double itime, ftime, exec_time;
    itime = omp_get_wtime(); 

	fprintf(datacsv, "Particles,x-position,y-position,time,%s\n", name);

	// initialization position and activity
	initialization(
		x, y, Particles,
		generator, distribution);

	check_nooverlap(
		x, y, Particles, L,
		generator, distribution);
	printf("Initialization done.\n");

	// Time evoultion
	int time;
	for (time = 0; time < N; time++)
	{
		update_position(
			x, y, phi, prefactor_e, Particles,
			delta, De, Dt, xi_e, xi_px,
			xi_py, vs, prefactor_xi_px, prefactor_xi_py,
			r, prefactor_interaction,
			generator, Gaussdistribution, distribution_e);
		if (strcmp(name, key1) == 0)
		{
			circular_reflective_boundary_conditions(
				x, y, Particles,
				Wall, L);
		}
		if (strcmp(name, key2) == 0)
		{
			reflective_boundary_conditions(
				x, y, Particles,
				Wall, L);
		}

		if (time % 100 == 0 && time >= 0)
		{
			print_file(
				x, y,
				Particles, time,
				datacsv);
		}
	}

	ftime = omp_get_wtime();
    exec_time = ftime - itime;
    printf("Time taken is %f", exec_time);

	free(x);
	free(y);

	fclose(datacsv);
	return 0;
}
