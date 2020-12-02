#include <iostream>
#include <fstream>
#include "solver.h"

using namespace std;
//using namespace arma;

int main(int argc, char* argv[])
{

	int N = atoi(argv[1]);
	double dx = 1./(N+1); // L = 1
	double dt = 0.5*dx*dx; // Stability condition for Forward Euler, Explicit
	double alpha = 0.5; // dt/(dx*dx)

	Solver solver; // Initializes Solver class

	// Initial conditions
	double *v = new double [N];
	for (int i = 0; i < N; i++){
		v[i] = -i*dx; // Add f(x)
		//cout << v[i] << endl;
	}
	double *u = new double [N];
	u = v;

	ofstream ofile;
	ofile.open("main.txt");

	for (int timestep = 0; timestep < N; timestep++){
		solver.Explicit(N, v, alpha); // Double loop
		//solver.Implicit(N, u, alpha); // Tridiagonal, doesn't work properly yet, memory issues with small N, returns -nan for big N
		
		//solver.Crank_Nicolson(N, v, alpha); // Tridiagonal
		//solver.WriteToFile();
		ofile << v[timestep] + timestep*dx << endl; // u = v - f(x)
	}

	delete[] u,v;


}