#ifndef SOLVER_H
#define SOLVER_H

using namespace std;

//#include "constants.h" // For geoscinece part later



#include <vector>
#include <string>

class Solver
{
public:
  void Forward_Sub(int N, double *a, double *b, double *c, double *g);
  void Backward_Sub(int N, double *b, double *c, double *g, double *v);
  void Explicit(int N, double *v, double alpha);
  void Implicit(int N, double *v, double alpha);
  void Crank_Nicolson(int N, double *v, double alpha);
  //void WriteToFile(string filename);


private:
};

#endif //SOLVER_H
