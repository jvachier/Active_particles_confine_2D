#include "headers/update_position.h"

using namespace std;

void update_position(
	double *x, double *y, double phi, double prefactor_e, int Particles,
	double delta, double De, double Dt, double xi_e, double xi_px,
	double xi_py, double vs, double prefactor_xi_px, double prefactor_xi_py,
	double r, double prefactor_interaction,
	default_random_engine &generator, normal_distribution<double> &Gaussdistribution, uniform_real_distribution<double> &distribution_e)
{
	double a = 0.0; // local variable - here check if no conflict elsewhere
	double F = 0.0, R = 0.0;
#pragma omp parallel for simd
	for (int k = 0; k < Particles; k++)
	{
		xi_e = distribution_e(generator);
		xi_px = Gaussdistribution(generator);
		xi_py = Gaussdistribution(generator);

		phi = (prefactor_e * xi_e); // be careful with radian and degree
		F = 0.0;
		for (int j = 0; j < Particles; j++)
		{
			if (k != j) // see how to improved the nested if conditions
			{
				R = sqrt((x[j] - x[k]) * (x[j] - x[k]) + (y[j] - y[k]) * (y[j] - y[k]));
				if (R < r)
				{
					a = prefactor_interaction / pow(R, 14);
					if (a > 1.0)
					{
						a = 1.0; // this value needs to be checked
					}
					F += a;
				}
			}
		}
		x[k] = x[k] + vs * cos(phi) * delta + F * x[k] * delta + xi_px * prefactor_xi_px;
		y[k] = y[k] + vs * sin(phi) * delta + F * y[k] * delta + xi_py * prefactor_xi_py;
	}
}