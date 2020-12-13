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
	return z*(Ta - Tb) - Ta;
}

double Q(double z, double L, double alpha_const){
	double Q;
	if 			(z*L <= 20e3)	{ Q = 1.4;  }
	else if (z*L <= 40e3)	{ Q = 0.35; }
	else if (z*L <= 120e3){ Q = 0.05; }
	else {
		cout << "input z for function Q(z) were out of range 0-120" << endl;
		exit(1);
	}
	return Q*1e-6*(L*L)/alpha_const; // TO BE DONE: Check that scale are correct
}

int main(int argc, char* argv[])
{
	double z;
	double L = 120e3; // [km]

	//Physical parameters
	double k = 2.5; 		// [W/(m°C)
	double rho = 3.5e3; // [Kg/m^3]
	double cp =  1000; 	// [J/(kg*°C)]
	double alpha_const = k/(rho*cp);
	double year_to_sec = 31556926; //sec pr year

	// int N = atoi(argv[1]);
	double dz = atof(argv[1])*1e3/L; // [km]
	double Time = atof(argv[2])*year_to_sec*1e6*alpha_const/(L*L); // [My]
	int N = 1/dz - 1;
	double dt = 0.5*dz*dz; // Stability condition for Forward Euler, Explicit
	double alpha = 0.5; // dt/(dz*dz)
	int timesteps = int(Time/dt) + 1;



	// Initial conditons
	double Ta = 8; 		// [°C]
	double Tb = 1300; // [°C]


	Solver solver; // Initializes Solver class
	// Intialize 1D solution vectors
	vec v = zeros<vec>(N+2);
	vec v_new = zeros<vec>(N+2);
	vec Q_vec = zeros<vec>(N+2);

	// 1D Initial and bounadry conditions
	v[0] = 0; v[N+1] = 0;
	for (int i = 1; i < N+1; i++){
		z = i*dz;
		v[i] = f(z, Ta, Tb ,L); // Add f(x)
		Q_vec[i] = Q(z, L, alpha_const);
		// cout << Q_vec[i] << endl;
	}

	//Read initial state
	Solver.ReadState(infile, )

	// Open output file
	string method_name = "PlainDist";
	string outfile = to_string	(1) + "D" + method_name + argv[1] + ".txt";
	solver.WriteToFile(outfile, 0, v, N, Time, dt, dz, f, L, alpha_const, Ta, Tb);

	cout <<timesteps<< endl;
	// Main calculation loop
	for (int t = 1; t <= timesteps; t++){
		solver.Explicit(N, v, v_new, alpha); 					// Double loop
		solver.Add_QdT(v, Q_vec, dt, rho, cp);
		// solver.Implicit(N, v, v_new, alpha); 				// Tridiagonal
		// solver.Crank_Nicolson(N, v, v_new, alpha); 	// Tridiagonal
		solver.WriteToFile(outfile, t, v, N, Time, dt, dz, f, L, alpha_const, Ta, Tb);
	}
	solver.WriteLastState("SteadyState.txt", v, N, dz, f , L, alpha_const, Ta, Tb);





} //end of main
