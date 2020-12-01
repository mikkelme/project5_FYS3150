#include "solver.h"
#include "constants.h" // For geoscinece part later
#include <vector>
#include <iostream>

using namespace std;


void Solver::Init(int N, vector u, double dx, double dt){
	
	// Boundary conditions
	u(x,0) = 0;
	u(0,t) = 0;
	u(L,t) = 1;

	double alpha = dt/(dx*dx);

}


void Solver::Forward_Sub(int N, vector a, vector b, vector c, vector g){
	double k; // Reducing no. FLOPs
	for (int i = 1; i < N; i++){
		k = a(i)/b(i-1);
		b(i) = b(i) - k*c(i-1);
		g(i) = g(i) - k*g(i-1);
	}
}

void Solver::Backward_Sub(int N, vector b, vector c, vector g, vector v){
	//v(0) = v(N+1) = 0;
	v(N) = g(N)/b(N);
	for (int i = N-1; i > 0; i--){
		v(i) = (g(i) - c(i)*v(i+1))/b(i);
	}
}

void Solver::Explicit(){


	for (int i = 0; i < steps; i++){
		u[]
	}

	ut += (U(x, t + dt) - U(x, t))/dt

	uxx += (U(x + dx, t) - 2*U(x, t) +U(x - dx, t))/(dx*dx)
	alpha = dt/(dx*dx)
	uij += alpha*u[i-1,j] + (1-2*alpha)*u[i,j]+alpha*u[i+1,j]
}


void Solver::Implicit(){
	/* ut  += (U(x, t) - U(x, t - dt))/dt

	u[i,j-1] = -alpha*u[i-1,j] + (1+2*alpha)*u[i,j] - alpha*u[i+1,j]
	*/

	vector<N> a; fill(a.begin(),a.end(), -alpha); // Filling vector a with -alpha, a = [-alpha, ... ,-alpha]
	vector<N> b; fill(b.begin(),b.end(), 1+2*alpha);
	vector<N> c; fill(c.begin(),c.end(), -alpha);
	vector<N> g;
	vector<N> v;


	for (int t = 0; t> N; t++){
		Forward_Sub(N, a, b, c, g);
		Backward_Sub(N, b, c, g, v);
		u(t) = v;
	
	}
	

}

void::Crank_Nicolson(){

	u
}

void WriteToFile(string filename){
	ostream ofile;
	ofile.open(filename, ios::out | ios::app);

	ofile << setprecision(5) << v << endl;

	ofile.close();
}