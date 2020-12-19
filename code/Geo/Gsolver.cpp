#include "Gsolver.h"
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

void Solver::Explicit(int N, vec &v, vec &v_new, double alpha){
		for (int i = 1; i < N+1; i++){
			v_new[i] = alpha*v[i-1] + (1-2*alpha)*v[i]+alpha*v[i+1];
		}
		v = v_new;
}

void Solver::Add_Q(double t, vec &v, vec &Q_vec, double Q_rad (double, double, double, double), double dt, double dz, double rho, double cp, double alpha_const, double L){
	int N = size(v)[0]-2;
	double z;
	for (int i = 1; i < N+1; i++){
		z = i*dz;
		v[i] += (Q_vec[i] + Q_rad(z, t, L, alpha_const))*dt/(rho*cp);
	}
}


void Solver::Implicit(int N, vec &v, vec &v_new, double alpha){
	vec a = zeros<vec>(N+2);
	vec b = zeros<vec>(N+2);
	vec c = zeros<vec>(N+2);

	a[0] = c[0] = -alpha; b[0] = 1+2*alpha;
	for (int i = 1; i < N+1; i++){
		a[i] = c[i] = -alpha;
		b[i] = 1+2*alpha;
	}
	a[N+1] = c[N+1] = -alpha; b[N+1] = 1+2*alpha;

	Forward_Sub(N, a, b, c, v);
	Backward_Sub(N, b, c, v, v_new);

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


	Forward_Sub(N, a, b, c, v_old);
	Backward_Sub(N, b, c, v_old, v);
}


void Solver::twoD_Explicit(int N, mat &u, mat &u_new, double alpha){
	for (int i = 1; i < N+1; i++){
		for (int j = 1; j < N+1; j++){
			u_new(i,j) = u(i,j) + alpha*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j));
		}
	}
	u = u_new;
}


void Solver::WriteToFile(string outfile, int t, vec &v, int N, double Time, double dt, double dz, double func (double, double, double, double), double L, double alpha_const, double Ta, double Tb){
	ofstream ofile;


	if (t == 0){
		ofile.open(outfile, ios::out);
		ofile << "N=" << N << endl;
		ofile << "Time=" << Time << endl;
		ofile << "dz=" << dz*L << endl;
		ofile << "dt=" << dt*(L*L)/(31556926*1e6*alpha_const) << endl;
	}

	else {
		ofile.open(outfile, ios::out | ios::app);
	}
	vec u = zeros<vec>(N+2);
	ofile << "t=" << t*dt*(L*L)/(31556926*1e6*alpha_const) << endl;
	for (int i = 0; i < N+2; i++){
		u[i] = v[i] - func(i*dz,Ta,Tb, L);
		ofile << setw(15) << setprecision(8) << u[i] << endl;
	}

	ofile.close();
}

void Solver::WriteLastState(string outfile, vec &v, int N, double dz, double func (double, double, double, double), double L, double alpha_const, double Ta, double Tb){
	ofstream ofile;

	ofile.open(outfile, ios::out);
	ofile << N << endl;
	ofile << dz*L << endl;

	vec u = zeros<vec>(N+2);
	for (int i = 0; i < N+2; i++){
		u[i] = v[i] - func(i*dz,Ta,Tb, L);
		ofile << setw(15) << setprecision(8) << u[i] << endl;
	}
}


vec Solver::ReadState(string filename, double L, double &dz, int &N){
	string line;
	ifstream infile (filename);

	getline(infile, line);
	N = stof(line);

	getline(infile, line);
	dz = stof(line)/L;

	// cout << N << " " << dz << endl;
	vec v = zeros<vec>(N+2);
	for (int i = 0; i < N+2; i++){
		getline(infile, line);
		v[i] = stof(line);
		// cout << v[i] << endl;
	}
	return v;

}
