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


// Steady state function
double f(double x, double L){
	return -x/L;
}


int main(int argc, char* argv[])
{
	double x,y;


	// double dx = 0.01;
	// double dt = atof(argv[1]);
	double dx = atof(argv[1]);
	double Time = atof(argv[2]);

	Solver solver; // Initializes Solver class
	int dim = 1; // Choose dimension of the problem (1 or 2)

	double L = 1;
	int N = L/dx - 1;
	double dt = 0.5*dx*dx/dim; // Stability condition for explicit
	double alpha = dt/(dx*dx);
	int timesteps = int(Time/dt) + 1;


	if (dim == 1){ //One dimensional problem
		// Intialize 1D solution vectors
		vec v = zeros<vec>(N+2);
		vec v_new = zeros<vec>(N+2);

		// 1D Initial and bounadry conditions
		v[0] = 0; v[N+1] = 0;
		for (int i = 1; i < N+1; i++){
			x = i*dx;
			v[i] = f(x,L); // Add f(x)
		}

		// Open output file
		string method_name = "CN";
		string outfile = to_string(dim) + "D" + method_name + argv[1] + ".txt";
		solver.WriteToFile1D(outfile, 0, v, N, Time, dt, dx, f, L);

		// Main calculation loop
		for (int t = 1; t < timesteps; t++){
			solver.Explicit(N, v, v_new, alpha); 						// Double loop
			// solver.Implicit(N, v, v_new, alpha); 				// Tridiagonal
			// solver.Crank_Nicolson(N, v, v_new, alpha); 	// Tridiagonal
			solver.WriteToFile1D(outfile, t, v, N, Time, dt, dx, f, L);
		}
	}


	if (dim == 2){ //Two dimensional problem
		// Intialize 2D solution vectors
		mat u = zeros<mat>(N+2,N+2);
		mat u_new = zeros<mat>(N+2,N+2);


		// 2D Initial and bounadry conditions
		double pi = acos(-1.0);
		for (int i = 1; i < N+1; i++){
			x = i*dx;
			for (int j = 1; j < N+1; j++){
				y = j*dx;
				u(i,j) = sin(pi/L*x)*sin(pi/L*y);
			}
		}


		// Open output file
		string method_name = "Explicit";
		string outfile = to_string(dim) + "D" + method_name + argv[1] + ".txt";
		solver.WriteToFile2D(outfile, 0,  u,  N, Time, dt, dx);

		// Main calculation loop
		for (int t = 1; t < timesteps; t++){
			solver.twoD_Explicit(N, u, u_new, alpha); // Triple loop
			solver.WriteToFile2D(outfile, t,  u,  N, Time, dt, dx);
		}
	}


} //end of main
