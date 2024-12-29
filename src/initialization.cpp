#include "headers/initialization.h"

using namespace std;

void initialization(
  double *x, double *y, int Particles,
  default_random_engine &generator,
  uniform_real_distribution<double> &distribution) {
#pragma omp parallel for simd
  for (int k = 0; k < Particles; k++) {
    x[k] = distribution(generator);
    y[k] = distribution(generator);
  }
}
