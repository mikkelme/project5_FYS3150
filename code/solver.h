#ifndef SOLVER_H
#define SOLVER_H
#include <armadillo>

using namespace std;
using namespace arma;

//#include "constants.h" // For geoscinece part later



#include <vector>
#include <string>

class Solver
{
public:
  void Forward_Sub(int N, double *a, double *b, double *c, double *g);
  void Backward_Sub(int N, double *b, double *c, double *g, double *v);
  void Explicit(int N, vec &v, vec &v_new, double alpha);
  void Implicit(int N, double *v, double alpha);
  void Crank_Nicolson(int N, double *v, double alpha);
  void WriteToFile(string outfile, int t, vec &v, int N, double Time, double dt, double dx, double func (double));


private:
};

#endif //SOLVER_H
