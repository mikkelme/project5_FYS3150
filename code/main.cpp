#include <cmath>
#include <iostream>
#include "solver.h"

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{

	N = atoi(argv[1]);
	double dx = 1./(n+1);
	double dt = 0.5*dx*dx; // Stability condition for Forward Euler, Explicit



	solver Solver;

	/*
	Initial conditions: u(i) = dx*i
	*/

	solver.Init();

	for (int timestep = 1; timestep < N; timestep++){
		//solver.Explicit(); // Double loop
		solver.Implicit(); // Tridiagonal
		//solver.Crank_Nicolson(); // Tridiagonal
		solver.WriteToFile();
	}




}