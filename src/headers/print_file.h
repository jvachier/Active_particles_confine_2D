#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <random>
#include <cstring>
#include <cmath>

void print_file(
  double *x, double *y,
  int Particles, int time,
  FILE *datacsv);
