#include <iostream>
#include <fstream>
#include <armadillo>
#include <string>
#include "solver.h"
using namespace std;
using namespace arma;

// Disable warnings about unused code
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"


double f(double x, double L){
	return -x/L;
}


int main(int argc, char* argv[])
{
	double x,y;

	// int N = atoi(argv[1]);
	double dx = atof(argv[1]);
	double Time = atof(argv[2]);

	double L = 1;
	int N = L/dx - 1;
	double dt = 0.5*dx*dx; // Stability condition for Forward Euler, Explicit
	double alpha = 0.5; // dt/(dx*dx)
	int timesteps = int(Time/dt) + 1;

	Solver solver; // Initializes Solver class
	int dim = 1; // Choose dimension of the problem (1 or 2)


	if (dim == 1){
		// Intialize solution vectors
		vec v = zeros<vec>(N+2);
		vec v_new = zeros<vec>(N+2);

		// 1D Initial and bounadry conditions
		v[0] = 0; v[N+1] = 0;
		for (int i = 1; i < N+1; i++){
			x = i*dx;
			v[i] = f(x,L); // Add f(x)
		}

		// Open output file
		string method_name = "Explicit";
		string outfile = to_string(dim) + "D" + method_name + argv[1] + ".txt";
		solver.WriteToFile1D(outfile, 0, v, N, Time, dt, dx, f, L);

		// Main calculation loop
		for (int t = 1; t <= timesteps; t++){
			solver.Explicit(N, v, v_new, alpha); 					// Double loop
			// solver.Implicit(N, v, v_new, alpha); 				// Tridiagonal
			// solver.Crank_Nicolson(N, v, v_new, alpha); 	// Tridiagonal
			solver.WriteToFile1D(outfile, t, v, N, Time, dt, dx, f, L);
		}
	}


	if (dim == 2){
		// 2D matrix
		double pi = acos(-1.0);
		mat V = zeros<mat>(N+2,N+2);
		mat V_new = zeros<mat>(N+2,N+2);
		for (int i = 1; i < N+1; i++){
			for (int j = 1; j < N+1; j++){
				x = i*dx; y = j*dx;
				V[i,j] = sin(pi*x)*sin(pi*y);
			}
		}

		string method_name = "explicit";
			string outfile = to_string(dim) + "D" + method_name + to_string(dx) + ".txt";
		solver.WriteToFile2D(outfile, 0, V, N, Time, dt, dx);
		for (int t = 1; t <= timesteps; t++){
			solver.twoD_Explicit(N, V, V_new, alpha); // Triple loop
			solver.WriteToFile2D(outfile, t, V, N, Time, dt, dx);
		}
	}


} //end of main
