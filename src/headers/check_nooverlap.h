#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <random>
#include <cstring>
#include <cmath>



void check_nooverlap(
  double *x, double *y, int Particles,
  int L,
  std::default_random_engine &generator, 
  std::uniform_real_distribution<double> &distribution);
