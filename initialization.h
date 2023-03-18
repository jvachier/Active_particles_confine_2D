// #ifndef INITIALIZATION_H
// #define INITIALIZATION_H

//using namespace std;

void initialization(
	double *, double *, int ,
	std::default_random_engine, std::uniform_real_distribution<double>
);

void check_nooverlap(
	double *, double *, int ,
	double , int ,
	std::default_random_engine, std::uniform_real_distribution<double>
);

// #endif