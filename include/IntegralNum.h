#pragma once
#include "PsiNJ.h"

namespace alg {
class IntegralNum {
public:
  IntegralNum(long double rmin, long double rmax, int jmax);
  Psi getPhinPlus1();
  Psi getMnPlus1();
  void setZeroOrder(Psi m, Psi phi);
  void setHorizonRadius(long double r_horizon);
  void executeIntegralInfinityZero(int rzad_perturbacji, long double eps);
private:
  long double r_min_;
  long double r_max_;
  int    j_max_;
  Psi m_n_plus_1_;
  Psi m_n_;
  Psi phi_n_plus_1_;
  Psi phi_n_;
  long double dr_;
  long double radius_horizon_;
};
}