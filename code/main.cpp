#include <iostream>
#include <fstream>
#include <armadillo>
#include "solver.h"
using namespace std;
using namespace arma;


double f(double x){
	return -x;
}


int main(int argc, char* argv[])
{
	double x;

	int N = atoi(argv[1]);
	double Time = atof(argv[2]);

	double dx = 1./(N); // L = 1
	double dt = 0.5*dx*dx; // Stability condition for Forward Euler, Explicit
	double alpha = 0.5; // dt/(dx*dx)
	int timesteps = int(Time/dt) + 1;

	Solver solver; // Initializes Solver class

	// Intialize solution vectors
	vec v = zeros<vec>(N+1);
	vec v_new = zeros<vec> (N+1);


	// Initial and bounadry conditions
	v[0] = 0; v[N] = 0;
	for (int i = 1; i < N; i++){
		x = i*dx;
		v[i] = f(x); // Add f(x)
	}


	string outfile = "main.txt";
	solver.WriteToFile(outfile, 0, v, N, Time, dt, dx, f);
	for (int t = 1; t <= timesteps; t++){
		solver.Explicit(N, v, v_new, alpha); // Double loop
		solver.WriteToFile(outfile, t, v, N, Time, dt, dx, f);

		//solver.Implicit(N, u, alpha); // Tridiagonal, doesn't work properly yet, memory issues with small N, returns -nan for big N
		//solver.Crank_Nicolson(N, v, alpha); // Tridiagonal
		// ofile << v[timestep] + timestep*dx << endl; // u = v - f(x)

	}
}
