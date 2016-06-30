#include <cmath>
#include "PsiNJ.h"
#include "ConvFactor.h"

namespace alg {
void alg::ConvFactor::setPsi4h(Psi psi4h){
  psi4h_ = psi4h;
}

void ConvFactor::setPsi2h(Psi psi2h){
  psi2h_ = psi2h;
}

void ConvFactor::setPsih(Psi psih){
  psih_ = psih;
}

long double ConvFactor::computeConvFactor(){
  
  long double err1 = 0;
  for(int i = 0; i <= psi4h_.getMaxj(); ++i){
    
    err1 += (psi4h_.getPsiElem(i) - psi2h_.getPsiElem(2*i)) *
            (psi4h_.getPsiElem(i) - psi2h_.getPsiElem(2*i));
  }
  long double norm_egz4_min_egz2 = sqrt(err1 / (psi4h_.getMaxj()+1));
  
  err1 = 0;
  for(int i = 0; i <= psi2h_.getMaxj(); ++i){
    err1 += (psi2h_.getPsiElem(i) - psih_.getPsiElem(2*i)) *
            (psi2h_.getPsiElem(i) - psih_.getPsiElem(2*i));
  }
  long double norm_egz2_min_egz1 = sqrt(err1 / (psi2h_.getMaxj()+1));
  
  
  return norm_egz4_min_egz2 / norm_egz2_min_egz1;
}

long double ConvFactor::computeConvFactorD(int rzad_poch){
  
  long double err1 = 0;
  for(int i = 0; i <= psi4h_.getMaxj(); ++i){
    
    err1 += (psi4h_.getDPsiElem(i, rzad_poch) - psi2h_.getDPsiElem(2*i, rzad_poch)) *
            (psi4h_.getDPsiElem(i, rzad_poch) - psi2h_.getDPsiElem(2*i, rzad_poch));
  }
  long double norm_egz4_min_egz2 = sqrt(err1 / (psi4h_.getMaxj()+1));
  
  err1 = 0;
  for(int i = 0; i <= psi2h_.getMaxj(); ++i){
    err1 += (psi2h_.getDPsiElem(i, rzad_poch) - psih_.getDPsiElem(2*i, rzad_poch)) *
            (psi2h_.getDPsiElem(i, rzad_poch) - psih_.getDPsiElem(2*i, rzad_poch));
  }
  long double norm_egz2_min_egz1 = sqrt(err1 / (psi2h_.getMaxj()+1));
  
  
  return norm_egz4_min_egz2 / norm_egz2_min_egz1;
}


}