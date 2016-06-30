#include <iostream>
#include <cmath>
#include "IntegralNum.h"
#include "TabFunc.h"

namespace alg {
IntegralNum::IntegralNum(long double rmin, long double rmax, int jmax)
    : m_n_plus_1_(jmax,(rmax - rmin) / jmax), m_n_(jmax,(rmax - rmin) / jmax), 
                  phi_n_plus_1_(jmax,(rmax - rmin) / jmax), 
                  phi_n_(jmax,(rmax - rmin) / jmax) {
  r_min_ = rmin;
  r_max_ = rmax;
  dr_ = (rmax - rmin) / jmax;
  j_max_ = jmax;
};
Psi IntegralNum::getMnPlus1() {
  return m_n_plus_1_;
}

Psi IntegralNum::getPhinPlus1() {
  return phi_n_plus_1_;
}


void IntegralNum::setZeroOrder(Psi m, Psi phi) {
  m_n_ = m;
  m_n_plus_1_ = m;
  phi_n_ = phi;
  phi_n_plus_1_ = phi;
}

void IntegralNum::setHorizonRadius(long double r_horizon){
  radius_horizon_ = r_horizon;
}

void IntegralNum::executeIntegralInfinityZero(int rzad_perturbacji, long double eps)
{
for(int pert = 1; pert <= rzad_perturbacji; ++pert){
  m_n_ =  m_n_plus_1_;
  phi_n_ = phi_n_plus_1_;
    
  for(int j = 0; j <= j_max_; ++j){
     
    if(j == 0){m_n_plus_1_.setPsiElem(j, 0);}
    else {
    long double right = m_n_plus_1_.getPsiElem(j-1) - eps/2.*dr_* (Ttt(m_n_.getPsiElem(j-1),
                                            m_n_.getDPsiElem(j-1,1),
                                            m_n_.getDPsiElem(j-1,2),m_n_.getDPsiElem(j-1,3),
                                            m_n_.getDPsiElem(j-1,4), phi_n_.getPsiElem(j-1),
                                            phi_n_.getDPsiElem(j-1,1),phi_n_.getDPsiElem(j-1,2),
                                            phi_n_.getDPsiElem(j-1,3),phi_n_.getDPsiElem(j-1,4),
                                            r_min_+dr_*(j-1)) + 
                                            Ttt(m_n_.getPsiElem(j),m_n_.getDPsiElem(j,1),
                                            m_n_.getDPsiElem(j,2),m_n_.getDPsiElem(j,3),
                                            m_n_.getDPsiElem(j,4), phi_n_.getPsiElem(j),
                                            phi_n_.getDPsiElem(j,1),phi_n_.getDPsiElem(j,2),
                                            phi_n_.getDPsiElem(j,3),phi_n_.getDPsiElem(j,4),
                                            r_min_+dr_*j));
    m_n_plus_1_.setPsiElem(j, right);
    }
  }
  for(int j = 0; j <= j_max_; ++j){
    long double right = m_n_plus_1_.getPsiElem(j) - m_n_plus_1_.getPsiElem(j_max_);
    m_n_plus_1_.setPsiElem(j, radius_horizon_/2. + right);
  }
}

for(int j = 0; j <= j_max_; ++j){
  if(j == 0){phi_n_plus_1_.setPsiElem(j, 0);}
  else {
  
    long double right = phi_n_plus_1_.getPsiElem(j-1) +
                      eps/2.* dr_* (TrrMTttPrzez(m_n_.getPsiElem(j-1),m_n_.getDPsiElem(j-1,1),
                                            m_n_.getDPsiElem(j-1,2),m_n_.getDPsiElem(j-1,3),
                                            m_n_.getDPsiElem(j-1,4), phi_n_.getPsiElem(j-1),
                                            phi_n_.getDPsiElem(j-1,1),phi_n_.getDPsiElem(j-1,2),
                                            phi_n_.getDPsiElem(j-1,3),phi_n_.getDPsiElem(j-1,4),
                                            r_min_+dr_*(j-1))+
                                            TrrMTttPrzez(m_n_.getPsiElem(j),m_n_.getDPsiElem(j,1),
                                            m_n_.getDPsiElem(j,2),m_n_.getDPsiElem(j,3),
                                            m_n_.getDPsiElem(j,4), phi_n_.getPsiElem(j),
                                            phi_n_.getDPsiElem(j,1),phi_n_.getDPsiElem(j,2),
                                            phi_n_.getDPsiElem(j,3),phi_n_.getDPsiElem(j,4),
                                            r_min_+dr_*j));
    phi_n_plus_1_.setPsiElem(j, right);
  }
}
  for(int j = 0; j <= j_max_; ++j){
    long double right = phi_n_plus_1_.getPsiElem(j) - phi_n_plus_1_.getPsiElem(j_max_);
    phi_n_plus_1_.setPsiElem(j, right);
  }
}

}