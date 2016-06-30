#pragma once
#include "PsiNJ.h"

namespace alg {
  
class ConvFactor {
public:
  ConvFactor():psih_(0,0),psi2h_(0,0), psi4h_(0,0) {
            
          };
  
  void setPsih(Psi psih);
  void setPsi2h(Psi psi2h);
  void setPsi4h(Psi psi4h);
  
  long double computeConvFactor();
  long double computeConvFactorD(int rzad_poch);
  
private:
  Psi psih_;
  Psi psi2h_;
  Psi psi4h_;
};
}