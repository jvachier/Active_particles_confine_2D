#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <random>
#include <cstring>
#include <cmath>

void initialization(
  double *x, double *y, int Particles,
  std::default_random_engine &generator, 
  std::uniform_real_distribution<double> &distribution);
