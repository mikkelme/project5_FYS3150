#include <iostream>
#include <fstream>
#include <armadillo>
#include <string>
#include "solver.h"
using namespace std;
using namespace arma;


double f(double x){
	return -x;
}


int main(int argc, char* argv[])
{
	double x,y;

	int N = atoi(argv[1]);
	double Time = atof(argv[2]);

	double dx = 1./(N+1); // L = 1
	double dt = 0.5*dx*dx; // Stability condition for Forward Euler, Explicit
	double alpha = 0.5; // dt/(dx*dx)
	int timesteps = int(Time/dt) + 1;

	Solver solver; // Initializes Solver class

	// Intialize solution vectors
	vec v = zeros<vec>(N+2);
	vec v_new = zeros<vec>(N+2);


	// 1D Initial and bounadry conditions
	v[0] = 0; v[N+1] = 0;
	for (int i = 1; i < N+1; i++){
		x = i*dx;
		v[i] = f(x); // Add f(x)
	}
	/*
	// 2D matrix, used same IC because idk what they are
	mat V = zeros<mat>(N+2,N+2);
	mat V_new = zeros<mat>(N+2,N+2);
	for (int i = 1; i < N; i++){
		x = i*dx; // dx = dy = h
		V[i,i] = f(x); 
	}*/

	// 2D matrix
	mat V = zeros<mat>(N+2,N+2);
	mat V_new = zeros<mat>(N+2,N+2);
	for (int i = 1; i < N; i++){
		for (int j = 1; j < N; j++){
			x = i*dx; y = j*dx;
			V[i,j] = sin(x)*sin(y);
		}
		//V[i,0] = 1;j =
		//V[i,0] = 1;
	}
	
	string outfile = "main.txt";
	solver.WriteToFile(outfile, 0, v, V, N, Time, dt, dx, f, 2);
	for (int t = 1; t <= timesteps; t++){
		//solver.Explicit(N, v, v_new, alpha); // Double loop
		//solver.Implicit(N, v, v_new, alpha); // Tridiagonal
		//solver.Crank_Nicolson(N, v, v_new, alpha); // Tridiagonal
		solver.twoD_Explicit(N, V, V_new, alpha); // Triple loop
		solver.WriteToFile(outfile, t, v, V, N, Time, dt, dx, f, 2);
	}
}
