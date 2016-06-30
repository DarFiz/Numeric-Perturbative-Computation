#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "PsiNJ.h"
#include "TaylorCoeffExact.h"
#include "PolynomialCurveFitting.h"
#include "IntegralNum.h"
#include "TabFunc.h"
#include "gnuplot_only_declaration.h"

namespace alg {

TaylorCoeffExact::TaylorCoeffExact(long double rmin, long double 
                                           rmax, long double horizon)
    :r_min_(rmin), r_max_(rmax), r_horizon_(horizon){

}


Psi TaylorCoeffExact::computeM0(int jmax) {

  long double h = (r_max_-r_min_)/jmax;

  Psi m_zero_order1h(jmax, h);
  for(int j = 0; j <= jmax; ++j){
  m_zero_order1h.setPsiElem(j,r_horizon_/2.);
  }
  return m_zero_order1h;
}

Psi TaylorCoeffExact::computePhi0(int jmax){

  long double h = (r_max_-r_min_) / jmax;

  Psi phi_zero_order1h(jmax, h);
  for(int j = 0; j <= jmax; ++j){
    phi_zero_order1h.setPsiElem(j,0);
  }
  return phi_zero_order1h;
}


Psi TaylorCoeffExact::computeM1DN(int jmax, long double eps, int rzad_poch){

  long double h = (r_max_-r_min_) / jmax;
  
  Psi m_first_oder_exact(jmax, h);
  for(int j = 0; j <= jmax; ++j){
    if (rzad_poch == 0){m_first_oder_exact.setPsiElem(j,M1(eps, r_horizon_/2.,r_min_ + j*h));}
    if (rzad_poch == 1){m_first_oder_exact.setPsiElem(j,M1D1(eps, r_horizon_/2.,r_min_ + j*h));}
    if (rzad_poch == 2){m_first_oder_exact.setPsiElem(j,M1D2(eps, r_horizon_/2.,r_min_ + j*h));}
    if (rzad_poch == 3){m_first_oder_exact.setPsiElem(j,M1D3(eps, r_horizon_/2.,r_min_ + j*h));}
    if (rzad_poch == 4){m_first_oder_exact.setPsiElem(j,M1D4(eps, r_horizon_/2.,r_min_ + j*h));}
  }
  return m_first_oder_exact;
}

Psi TaylorCoeffExact::computePhi1DN(int jmax, long double eps, int rzad_poch){

  long double h = (r_max_-r_min_) / jmax;
  
  Psi m_first_oder_exact(jmax, h);
  for(int j = 0; j <= jmax; ++j){
    if (rzad_poch == 0){m_first_oder_exact.setPsiElem(j,Phi1(eps, r_horizon_/2.,r_min_ + j*h));}
    if (rzad_poch == 1){m_first_oder_exact.setPsiElem(j,Phi1D1(eps, r_horizon_/2.,r_min_ + j*h));}
    if (rzad_poch == 2){m_first_oder_exact.setPsiElem(j,Phi1D2(eps, r_horizon_/2.,r_min_ + j*h));}
    if (rzad_poch == 3){m_first_oder_exact.setPsiElem(j,Phi1D3(eps, r_horizon_/2.,r_min_ + j*h));}
    if (rzad_poch == 4){m_first_oder_exact.setPsiElem(j,Phi1D4(eps, r_horizon_/2.,r_min_ + j*h));}
  }
  return m_first_oder_exact;
}


Psi TaylorCoeffExact::computeM2(int jmax, long double eps){
  
  long double h = (r_max_-r_min_)/jmax;
  
  Psi m_first_oder_exact(jmax, h);
  for(int j = 0; j <= jmax; ++j){
    m_first_oder_exact.setPsiElem(j,M2(r_horizon_/2.,r_min_ + j*h));
  }
  return m_first_oder_exact;
}

Psi TaylorCoeffExact::computeM3(int jmax, long double eps){
  
  long double h = (r_max_-r_min_)/jmax;
  
  Psi m_first_oder_exact(jmax, h);
  for(int j = 0; j <= jmax; ++j){
    m_first_oder_exact.setPsiElem(j,M3(r_horizon_/2.,r_min_ + j*h));
  }
  return m_first_oder_exact;
}


Psi TaylorCoeffExact::computePhi3(int jmax, long double eps){
  
  long double h = (r_max_-r_min_)/jmax;
  
  Psi m_first_oder_exact(jmax, h);
  for(int j = 0; j <= jmax; ++j){
    m_first_oder_exact.setPsiElem(j,Phi3(r_horizon_/2.,r_min_ + j*h));
  }
  return m_first_oder_exact;
}


Psi TaylorCoeffExact::computePhi2(int jmax, long double eps) {
  
  long double h = (r_max_-r_min_)/jmax;
  
  Psi phi_first_oder_exact(jmax, h);
  for(int j = 0; j <= jmax; ++j){
    phi_first_oder_exact.setPsiElem(j,Phi2( r_horizon_/2.,r_min_ + j*h));
  }
  
  return phi_first_oder_exact;
}

}