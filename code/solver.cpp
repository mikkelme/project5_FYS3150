#include "solver.h"
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;


void Solver::Forward_Sub(int N, vec &a, vec &b, vec &c, vec &g){
	double k; // Reducing no. FLOPs
	for (int i = 1; i < N+1; i++){
		k = a[i]/b[i-1];
		b[i] = b[i] - k*c[i-1];
		g[i] = g[i] - k*g[i-1];
	}
}

void Solver::Backward_Sub(int N, vec &b, vec &c, vec &g, vec &v){
	v[N+1] = g[N+1]/b[N+1];
	for (int i = N+1; i > 0; i--){
		v[i] = (g[i] - c[i]*v[i+1])/b[i];
	}
}


void Solver::tridag(int N, vec &a, vec &b, vec &c, vec &g, vec &v){
	//Forward Sub
	double k; // Reducing no. FLOPs
	for (int i = 2; i < N+1; i++){
		k = a[i]/b[i-1];
		b[i] = b[i] - k*c[i-1];
		g[i] = g[i] - k*g[i-1];
	}

	//Backwards Sub
	v[0] = v[N+1] = 0;
	v[N] = g[N]/b[N];
	for (int i = N-1; i > 0; i--){
		v[i] = (g[i] - c[i]*v[i+1])/b[i];
	}
}



void Solver::Explicit(int N, vec &v, vec &v_new, double alpha){

		for (int i = 1; i < N+1; i++){
			v_new[i] = alpha*v[i-1] + (1-2*alpha)*v[i]+alpha*v[i+1];
		}
		v = v_new;
}


void Solver::Implicit(int N, vec &v, vec &v_new, double alpha){

	vec a = zeros<vec>(N+2);
	vec b = zeros<vec>(N+2);
	vec c = zeros<vec>(N+2);

	for (int i = 0; i < N+2; i++){
		a[i] = c[i] = -alpha;
		b[i] = 1+2*alpha;
	}

	tridag(N, a, b, c, v, v_new);

	v = v_new;

}


void Solver::Crank_Nicolson(int N, vec &v, vec &v_old, double alpha){

	vec a = zeros<vec>(N+2);
	vec b = zeros<vec>(N+2);
	vec c = zeros<vec>(N+2);

	a[0] = c[0] = -alpha; b[0] = 2+2*alpha;
	for (int i = 1; i < N+1; i++){
		a[i] = c[i] = -alpha;
		b[i] = 2+2*alpha;
	}
	a[N+1] = c[N+1] = -alpha; b[N+1] = 2+2*alpha;

	for (int i = 1; i < N+1; i++){
		v_old[i] = alpha*v[i-1] + (2-2*alpha)*v[i]+alpha*v[i+1];
	}

	tridag(N, a, b, c, v_old, v);
}


void Solver::twoD_Explicit(int N, mat &u, mat &u_new, double alpha){
	for (int i = 1; i < N+1; i++){
		for (int j = 1; j < N+1; j++){
			u_new(i,j) = u(i,j) + alpha*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j));
		}
	}
	u = u_new;
}


void Solver::WriteToFile1D(string outfile, int t, vec &v, int N, double Time, double dt, double dx, double func (double, double), double L){
	ofstream ofile;

	if (t == 0){
		ofile.open(outfile, ios::out);
		ofile << "N=" << N << endl;
		ofile << "Time=" << Time << endl;
		ofile << "dx=" << dx << endl;
		ofile << "dt=" << dt << endl;
	}
	else {
		ofile.open(outfile, ios::out | ios::app);
	}

	vec u = zeros<vec>(N+2);
	ofile << "t=" << t*dt << endl;
	for (int i = 0; i < N+2; i++){
		u[i] = v[i] - func(i*dx, L);
		ofile << setw(15) << setprecision(8) << u[i] << endl;
	}

	ofile.close();
}



void Solver::WriteToFile2D(string outfile, int t, mat &u, int N, double Time, double dt, double dx){
	ofstream ofile;

	if (t == 0){
		ofile.open(outfile, ios::out);
		ofile << "N=" << N << endl;
		ofile << "Time=" << Time << endl;
		ofile << "dx=" << dx << endl;
		ofile << "dt=" << dt << endl;
	}
	else {
		ofile.open(outfile, ios::out | ios::app);
	}

	ofile << "t=" << t*dt << endl;
	for (int i = 0; i < N+2; i++){
		for (int j = 0; j < N+2; j++){
			ofile << setw(15) << setprecision(6) << u(i,j) << endl;
		}
	}
	ofile.close();
}
