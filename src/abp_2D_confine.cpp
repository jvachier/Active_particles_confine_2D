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
#include <cstring>
#include <stdio.h>
#include <cmath>
#include <time.h>
#include <omp.h> //import library to use pragma
#include <tuple> //to output multiple components of a function

#define PI 3.141592653589793
#define N_thread 6

using namespace std;

void update_position(
	double *x, double *y, double phi, double prefactor_e, int Particles, 
	double delta, double De, double Dt, double xi_e, double xi_px, 
	double xi_py, double vs, double prefactor_xi_px, double prefactor_xi_py,
	double r, double R, double F, double prefactor_interaction,
	default_random_engine &generator, normal_distribution<double> &Gaussdistribution, uniform_real_distribution<double> &distribution_e
)
{
	double a = 0.0; // local variable - here check if no conflict elsewhere
#pragma omp parallel for simd num_threads(N_thread)
	for (int k = 0; k < Particles; k++)
	{
		xi_e = distribution_e(generator);
		xi_px = Gaussdistribution(generator);
		xi_py = Gaussdistribution(generator);

		phi = (prefactor_e * xi_e); // be careful with radian and degree
		F = 0.0;
		for (int j = 0; j < Particles; j++)
		{
			if (k!=j) // see how to improved the nested if conditions
			{
				R = sqrt((x[j]-x[k])*(x[j]-x[k]) + (y[j]-y[k])*(y[j]-y[k]));
				if (R < r)
				{
					a = prefactor_interaction / pow(R,14);
					if (a > 1.0)
					{
						a = 0.5; // this value needs to be checked 
					}
					F += a;
				}
			}
		}
		x[k] = x[k] + vs * cos(phi) * delta + F * x[k] * delta + xi_px * prefactor_xi_px;
		y[k] = y[k] + vs * sin(phi) * delta + F * y[k] * delta  + xi_py * prefactor_xi_py;
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
		fprintf(datacsv,"Particles%d,%lf,%lf,%d\n",k,x[k],y[k],time);
	}
}

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

void circular_reflective_boundary_conditions(
	double *x, double *y, int Particles,
	double Wall, int L
)
{
	double distance_squared = 0.0, Wall_squared = Wall * Wall;
#pragma omp parallel for simd num_threads(N_thread)	
	for (int k = 0; k < Particles; k++)
	{
		distance_squared = x[k]*x[k] + y[k]*y[k];
		if (distance_squared > Wall_squared)
		{
			x[k] = (sqrt(Wall_squared) / sqrt(distance_squared)) * x[k];
			y[k] = (sqrt(Wall_squared) / sqrt(distance_squared)) * y[k];
		}
	}
}

void reflective_boundary_conditions(
	double *x, double *y, int Particles,
	double Wall, int L
)
{
	double D_AW_x = 0.0;
	double D_AW_y = 0.0;
#pragma omp parallel for simd num_threads(N_thread)	
	for (int k = 0; k < Particles; k++)
	{
		D_AW_x = 0.0;
		D_AW_y = 0.0;
		if (abs(x[k]) > Wall)
		{
			D_AW_x = abs(x[k] + Wall); 

			if (D_AW_x > 4.0 * L )
			{
				if (x[k] > Wall)
				{
					x[k] = Wall - 2.0 * L;
				}
				else if (x[k] < -Wall)
				{
					x[k] = 2.0 * L - Wall;
				}
			}
			else
			{
				if (x[k] > Wall)
				{
					x[k] -= 2.0 * D_AW_x;
				}
				else if (x[k] < -Wall)
				{
					x[k] += 2.0 * D_AW_x;
				}
			}
			
		}
		if (abs(y[k]) > Wall)
		{
			D_AW_y = abs(y[k] + Wall);
			if (D_AW_y > 4.0 * L )
			{
				if (y[k] > Wall)
				{
					y[k] = Wall - 2.0 * L;
				}
				else if (y[k] < -Wall)
				{
					y[k] = 2.0 * L - Wall;
				}
			}
			else
			{
				if (y[k] > Wall)
				{
					y[k] -= 2.0 * D_AW_y;
				}
				else if (y[k] < -Wall)
				{
					y[k] += 2.0 * D_AW_y;
				}
			}
			
		}
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
	datacsv = fopen("./data/simulation_test.csv", "w");

	// check if the file parameter is exist
	if (parameter == NULL)
	{
		printf("no such file.");
		return 0;
	}

	// read the parameters from the file
	double epsilon, delta, Dt, De, vs;
	double F, R, Wall;
	int Particles;
	char name[100];
	char key1[] = "circular";
	char key2[] = "squared";

	bool flag = false;

	printf("Select confinement geometry, either squared or circular:");
	scanf("%s",&name);

	while(flag == false){
		if ((strcmp(name,key1) == 0) or (strcmp(name,key2) == 0)){
			flag = true;
		}
		else{
			printf("You have not selected the correct, please select again\n");
			printf("Select confinement geometry, either squared or circular:");
			scanf("%s",&name);
			flag = false;
		}
	}

	fscanf(parameter, "%lf\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\n", &epsilon, &delta, &Particles, &Dt, &De, &vs, &Wall);
	printf("%lf\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\n", epsilon, delta, Particles, Dt, De, vs, Wall);

	double *x = (double *)malloc(Particles * sizeof(double)); // x-position 
	double *y = (double *)malloc(Particles * sizeof(double)); // y-position

	// parameters
	const int N = 1E5; // number of iterations
	const int L = 1.0; // particle size

	// initialization of the random generator
	random_device rdev;
	default_random_engine generator(rdev()); // random seed -> rdev
	// default_random_engine generator(1); //same seed

	// Distributions Gaussian
	normal_distribution<double> Gaussdistribution(0.0, 1.0);
	// Distribution Uniform for initialization
	uniform_real_distribution<double> distribution(-Wall,Wall);
	//uniform_real_distribution<double> distribution_e(0.0,360.0*PI / 180.0); // directly in radian
	uniform_real_distribution<double> distribution_e(0.0,360.0); 

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
	double prefactor_interaction = epsilon * 48.0;
	double r = 5.0 * L;



	clock_t tStart = clock(); // check time for one trajectory

	fprintf(datacsv,"Particles,x-position,y-position,time\n");

	// initialization position and activity
	initialization(
		x, y, Particles,
		generator, distribution
	);

	check_nooverlap(
		x, y, Particles,
		R, L,
		generator, distribution
	);
	printf("Initialization done.\n");

	// Time evoultion
	int time;
	for (time = 0; time < N; time++)
	{
		update_position(
			x, y, phi, prefactor_e, Particles, 
			delta, De, Dt, xi_e, xi_px, 
			xi_py, vs, prefactor_xi_px, prefactor_xi_py,
			r, R, F, prefactor_interaction,
			generator, Gaussdistribution, distribution_e
		);
		if (strcmp(name,key1) == 0){
			circular_reflective_boundary_conditions(
				x, y, Particles,
				Wall, L
			);
		}
		if (strcmp(name,key2) == 0){
			reflective_boundary_conditions(
				x, y, Particles,
				Wall, L
			);
		}
		
		
		if(time % 100 == 0 && time >= 0)
		{
			print_file(
				x, y,
				Particles, time,
				datacsv
			);
		}
	}


	

	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC); // time for one trajectory

	free(x);
	free(y);

	fclose(datacsv);
	return 0;
}
