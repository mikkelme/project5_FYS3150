#ifndef SOLVER_H
#define SOLVER_H
#include <armadillo>

using namespace std;
using namespace arma;

//#include "constants.h" // For geoscinece part later

#include <string>

class Solver
{
public:
  void Forward_Sub(int N, vec &a, vec &b, vec &c, vec &g);
  void Backward_Sub(int N, vec &b, vec &c, vec &g, vec &v);
  void Explicit(int N, vec &v, vec &v_new, double alpha);
  void Implicit(int N, vec &v, vec &v_new, double alpha);
  void Crank_Nicolson(int N, vec &v, vec &v_old, double alpha);
  void twoD_Explicit(int N, mat &V, mat &V_new, double alpha);
  void WriteToFile1D(string outfile, int t, vec &v, int N, double Time, double dt, double dx, double func (double, double), double L);
  void WriteToFile2D(string outfile, int t, mat &u, int N, double Time, double dt, double dx);


private:
};

#endif //SOLVER_H
