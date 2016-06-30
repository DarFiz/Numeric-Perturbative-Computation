#pragma once
#include <vector>
#include "PsiNJ.h"
#include <armadillo>

namespace alg {
// using namespace  arma;

class PolynomialCurveFitting {
public:
  PolynomialCurveFitting(int jmax, long double rmin, long double h, int min_rzad_poly, int max_rzad_poly, Psi y);
  void Execute();
  Psi getWynik();
  
// private:
  void constructMatrixA();
  void constructMatrixB();
  long double SumPowX(int potega);
  long double SumPowXY(int potega);
  Psi x_;
  Psi y_;
  std::vector<long double> X_;

  arma::mat A_;
  arma::vec B_;
  
  int j_max_;
  long double h_;
  long double r_min_;
  int max_rzad_polynomial_;
  int min_rzad_polynomial_;
  Psi wynik_;
  std::vector<long double> v_wspocz_wielomianu_;
};
  
}