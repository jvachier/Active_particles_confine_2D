#include "print_file.h"

using namespace std;

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