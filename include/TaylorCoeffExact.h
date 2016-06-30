#pragma once
#include <vector>
#include <PsiNJ.h>

namespace alg {
  
class TaylorCoeffExact {
public:
  TaylorCoeffExact(long double rmin, long double rmax, long double horizon);

  Psi computeM0(int jmax);
  Psi computePhi0(int jmax);
  
  Psi computePhi1DN(int jmax, long double eps, int rzad_poch);
  Psi computeM1DN(int jmax, long double eps, int rzad_poch);
  
  Psi computeM2(int jmax, long double eps);
  Psi computePhi2(int jmax, long double eps);

  Psi computeM3(int jmax, long double eps);
  Psi computePhi3(int jmax, long double eps);
  
  long double r_min_;
  long double r_max_;
  long double r_horizon_;
};
}