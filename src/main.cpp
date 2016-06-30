#include <iostream>
#include "PsiNJ.h"
#include "TabFunc.h"
#include "IntegralNum.h"
#include "PolynomialCurveFitting.h"
#include "ITaylorCoeff.h"
#include "TaylorCoeffClassic.h"
#include "TaylorCoeffExact.h"
#include "TaylorCoeffFit.h"
#include "FactoryMethodTaylorCoeff.h"
#include <cmath>
#include <limits>       // std::numeric_limits


int main(int argc, char **argv) {
  
  long  double r_horizon = 2;
  long  double rmin = r_horizon;
  long double rmax = 200.5;
  long double eps = 1./10000.;
  
  alg::FactoryMethodTaylorCoeff factory_taylor_coeff;
  std::unique_ptr<alg::ITaylorCoeff> up_itaylor_coeff = factory_taylor_coeff.factoryMethod("classic", rmin, rmax, r_horizon);
//   std::unique_ptr<alg::ITaylorCoeff> up_itaylor_coeff = factory_taylor_coeff.factoryMethod("fit", rmin, rmax, r_horizon);
  
  
  up_itaylor_coeff->plotM1vsM1ExactoPs(0.02, 2.5, eps);
  up_itaylor_coeff->plotM1vsM1ExactoPs(0.02, 10., eps);
  
  up_itaylor_coeff->plotPhi1vsPhi1ExactoPs(0.02, 2.5, eps);
  up_itaylor_coeff->plotPhi1vsPhi1ExactoPs(0.02, 10., eps);

  up_itaylor_coeff->plotErrorM1DNToPs(0.0001, 0.002, 0.0001, eps, 0);
  up_itaylor_coeff->plotErrorPhi1DNToPs(0.0001, 0.002, 0.0001, eps, 0);
  
  up_itaylor_coeff->plotM2vsM2ExactoPs(0.02, 2.5, eps);
  up_itaylor_coeff->plotM2vsM2ExactoPs(0.02, 10., eps);
  
  up_itaylor_coeff->plotPhi2vsPhi2ExactoPs(0.003, 2.5, eps);
  up_itaylor_coeff->plotPhi2vsPhi2ExactoPs(0.003, 10., eps);
  
  
  up_itaylor_coeff->plotErrorM2ToPs(0.001, 0.007, 0.0001, eps);
  up_itaylor_coeff->plotErrorPhi2ToPs(0.0005, 0.007, 0.0001, eps);
  
  alg::TaylorCoeffFit * p_taylor_fit ;
  
  if ((p_taylor_fit = dynamic_cast<alg::TaylorCoeffFit*>(up_itaylor_coeff.get()))) {
  
  p_taylor_fit->plotM1AnCoeff(0.01, 0.1, 0.01, eps, -108, 5);
  p_taylor_fit->plotM1AnCoeff(0.01,0.1, 0.01, eps, 196, 6);
  p_taylor_fit->plotM2AnCoeff(0.01, 0.1, 0.01, eps, -98496./11., 11);
  p_taylor_fit->plotM2AnCoeff(0.01, 0.1, 0.01, eps, 14808, 12);
    
  p_taylor_fit->plotPhi1AnCoeff(0.01,0.1, 0.01, eps, -108, 6);
  p_taylor_fit->plotPhi2AnCoeff(0.01, 0.1, 0.01, eps, 165888/11, 11);
  p_taylor_fit->plotPhi2AnCoeff(0.01, 0.1, 0.01, eps, -29808, 12);
  }
  
  return 0;
}
