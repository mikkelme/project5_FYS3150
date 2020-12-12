#include <iostream>
#include <fstream>
#include <armadillo>
#include <string>
#include "Gsolver.h"
using namespace std;
using namespace arma;

// Disable warnings about unused code
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"


double f(double z, double Ta, double Tb, double L){
	return z*(Ta - Tb)/L - Ta;
}

double q(double z, double o){
	if (z <= 20){
		o = 1.4; // mu W/m^3
	}
	if (z <= 40){
		o = 0.35; 
	}
	if (z <= 120){
		o = 0.05;
	}	
	return o;

}

int main(int argc, char* argv[])
{
	double x,y,z, o;

	// int N = atoi(argv[1]);
	double dx = atof(argv[1]);
	double Time = atof(argv[2]);

	double L = 120; // [km]
	int N = L/dx - 1;
	double dt = 0.5*dx*dx; // Stability condition for Forward Euler, Explicit
	double alpha = 0.5; // dt/(dx*dx)
	int timesteps = int(Time/dt) + 1;


	//Physical parameters
	double k = 2.5; 		// [W/(m*°C)
	double rho = 3.5e3; // [Kg/m^3]
	double cp =  1000; 	// [J/(kg*°C)]
	alpha = k/(rho*cp);

	double Q1 = 0;
	double Q2 = 0;
	double Q3 = 0;

	// Initial conditons
	double Ta = 8; 		// [°C]
	double Tb = 1300; // [°C]


	Solver solver; // Initializes Solver class
	// Intialize 1D solution vectors
	vec v = zeros<vec>(N+2);
	vec v_new = zeros<vec>(N+2);

	// 1D Initial and bounadry conditions
	v[0] = 0; v[N+1] = 0;
	for (int i = 1; i < N+1; i++){
		z = i*dx;
		v[i] = f(z, Ta, Tb ,L) + q(z, o); // Add f(x)
		cout << z <<" " << o << endl;
	}

	// Open output file
	string method_name = "PlainDist";
	string outfile = to_string	(1) + "D" + method_name + argv[1] + ".txt";
	solver.WriteToFile1D(outfile, 0, v, N, Time, dt, dx, f, L);

	// Main calculation loop
	for (int t = 1; t <= timesteps; t++){
		solver.Explicit(N, v, v_new, alpha); 					// Double loop
		// solver.Implicit(N, v, v_new, alpha); 				// Tridiagonal
		// solver.Crank_Nicolson(N, v, v_new, alpha); 	// Tridiagonal
		solver.WriteToFile1D(outfile, t, v, N, Time, dt, dx, f, L);
	}


} //end of main
