#pragma once
#include <vector>
#include "PsiNJ.h"

namespace alg {
  
class Error {
public:
  
  static long double computeError(Psi psi1, Psi psi2);
  static long double computeErrorIvs2I(Psi psi1, Psi psi2);
  static long double computeErrorD(int rzad_poch, Psi psi1, Psi psi2);

};
}