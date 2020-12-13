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
  void Add_QdT(vec &v, vec &Q_vec, double dt, double rho, double cp);
  void twoD_Explicit(int N, mat &V, mat &V_new, double alpha);
  void WriteToFile(string outfile, int t, vec &v, int N, double Time, double dt, double dx, double func (double, double, double, double), double L, double alpha_const, double Ta, double Tb);
  void WriteLastState(string outfile, vec &v, int N, double dx, double func (double, double, double, double), double L, double alpha_const, double Ta, double Tb);


private:
};

#endif //SOLVER_H
