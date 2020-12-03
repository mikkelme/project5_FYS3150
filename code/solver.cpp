#include "solver.h"
#include <string>
#include <iostream>

using namespace std;



void Solver::Forward_Sub(int N, double *a, double *b, double *c, double *g){
	double k; // Reducing no. FLOPs
	for (int i = 1; i < N; i++){
		k = a[i]/b[i-1];
		b[i] = b[i] - k*c[i-1];
		g[i] = g[i] - k*g[i-1];
	}
}

void Solver::Backward_Sub(int N, double *b, double *c, double *g, double *v){
	//v[0] = v[N+1] = 0;
	v[N] = g[N]/b[N];
	for (int i = N-1; i > 0; i--){
		v[i] = (g[i] - c[i]*v[i+1])/b[i];
	}
}

void Solver::Explicit(int N, vec &v, vec &v_new, double alpha){

		for (int i = 1; i <= N; i++){
			v_new[i] = alpha*v[i-1] + (1-2*alpha)*v[i]+alpha*v[i+1];
		}
		v = v_new;

	/*
	double *v_new = new double[N];
	double *v_old; // Need of temporary memory so the memory doesnt overwrite


	for (int t = 0; t < N; t++){
		for (int i = 1; i < N; i++){
			v_new[i] = alpha*v[i-1] + (1-2*alpha)*v[i]+alpha*v[i+1];
		}

		v_old = v_new;
		v_new = v;
		v = v_old;
	}

	delete[] v_new;
	*/

}


void Solver::Implicit(int N, double *v, double alpha){

	double *a = new double [N];
	double *b = new double [N];
	double *c = new double [N];

	double *v_new = new double[N];
	double *v_old;

	for (int t = 0; t < N; t++){
		for (int i = 1; i < N; i++){
			a[i] = c[i] = -alpha;
			b[i] = 1+2*alpha;
		}
		Forward_Sub(N, a, b, c, v);
		Backward_Sub(N, b, c, v, v_new);

		// Keep boundary
		a[0] = c[0] = 0;
		b[0] = 1;
		a[N] = c[N] = 0;
		b[N] = 1;

		v_old = v_new;
		v_new = v;
		v = v_old;
	}


	// delete[] a, b, c, v_old;

}

void Solver::Crank_Nicolson(int N, double *v, double alpha){

	double *a = new double [N];
	double *b = new double [N];
	double *c = new double [N];

	double *v_new = new double[N];
	double *v_old = new double[N];

	for (int t = 0; t < N; t++){
		for (int i = 1; i < N-1; i++){
			v_old[i] = alpha*v[i-1] + (2-2*alpha)*v[i]+alpha*v[i+1];

		}


		for (int i = 0; i < N; i++){
			a[i] = c[i] = -alpha;
			b[i] = 2+2*alpha;
		}
		Forward_Sub(N, a, b, c, v_old);
		Backward_Sub(N, b, c, v_old, v);

		// Keep boundary
		a[0] = c[0] = 0;
		b[0] = 1;
		a[N] = c[N] = 0;
		b[N] = 1;

		v_new = v_old;
		v_old = v;
		v = v_new;
	}


	//delete[] a, b, c, v_new;

}

void Solver::WriteToFile(string outfile, int t, vec &v, int N, double Time, double dt, double dx, double func (double)){
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

	vec u;
	ofile << "t=" << t*dt << endl;
	for (int i = 0; i <= N+1; i++){
		u[i] = v[i] - func(i*dx);
		ofile << setw(15) << setprecision(8) << u[i] << endl;
	}
}







//
