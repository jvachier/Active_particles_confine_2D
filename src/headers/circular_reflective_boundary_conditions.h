#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <random>
#include <cstring>
#include <cmath>

void circular_reflective_boundary_conditions(
  double *x, double *y, int Particles,
  double Wall, int L);
