#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include "PsiNJ.h"
#include "Error.h"

namespace alg {

long double Error::computeError(Psi psi1, Psi psi2) {
  int jmax = psi1.getMaxj();

  long double err = 0;
  for(int i = 0; i <= jmax; ++i){
    err += sqrt((psi1.getPsiElem(i) - psi2.getPsiElem(i)) *
                (psi1.getPsiElem(i) - psi2.getPsiElem(i)));
  }
  std::cout << err / (jmax+1)<<std::endl; 
  
  return err/ (jmax+1);

}

long double Error::computeErrorIvs2I(Psi psiI, Psi psi2I){
  
  int jI = psiI.getMaxj();

  long double err = 0;
  for(int i = 0; i <= jI; ++i){
    err += sqrt((psiI.getPsiElem(i) - psi2I.getPsiElem(2*i)) *
                (psiI.getPsiElem(i) - psi2I.getPsiElem(2*i)));
  }
  std::cout << err / (jI+1)<<std::endl; 
  
  return err/ (jI+1);
}

long double Error::computeErrorD(int rzad_poch, Psi psi1, Psi psi_pochodnej){

  int jmax = psi1.getMaxj();
  long double err = 0;
  for(int i = 0; i <= jmax; ++i){
    err += sqrt((psi1.getDPsiElem(i, rzad_poch) - psi_pochodnej.getPsiElem(i)) *
                (psi1.getDPsiElem(i, rzad_poch) - psi_pochodnej.getPsiElem(i)));
  }
  std::cout << err / (jmax+1) << std::endl; 

  return err / (jmax+1);
}

}