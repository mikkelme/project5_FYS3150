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
  void WriteToFile(string outfile, int t, vec &v, int N, double Time, double dt, double dx, double func (double));


private:
};

#endif //SOLVER_H
