#include "PsiNJ.h"

Psi::Psi(int maxj, long double h)
    :h_(h), j_max_(maxj) {
  v_psi_.resize(maxj + 1);
}

Psi Psi::operator+(const Psi& psi1) {
  Psi psi_egz(j_max_, h_);
  for(int j = 0; j<=j_max_; ++j ) {
    psi_egz.v_psi_[j] = v_psi_[j] + psi1.v_psi_[j];
  }
  
  return psi_egz;
}

Psi Psi::operator-(const Psi& psi1) {
  Psi psi_egz(j_max_, h_);
  for(int j = 0; j<=j_max_; ++j ) {
    psi_egz.v_psi_[j] = v_psi_[j] - psi1.v_psi_[j];
  }
  
  return psi_egz;
}

Psi Psi::operator*(const long double float_number) {
  Psi psi_egz(j_max_, h_);
  for(int j = 0; j<=j_max_; ++j ) {
    psi_egz.v_psi_[j] = float_number*v_psi_[j];
  }
  
  return psi_egz;
}


long double Psi::getPsiElem(int j) {
  return v_psi_[j];
}
int Psi::getMaxj() {
  return j_max_;
}

long double Psi::geth() {
  return h_;
}


long double Psi::getDPsiElem(int j, int rzadpoch) {
  
  switch(rzadpoch) {
    case 1:
        if(j == 0) {return (-3/2.*getPsiElem(0)+2.*getPsiElem(1)-1/2.*getPsiElem(2))/(h_);}
        else if(j == j_max_){return (1/2.*getPsiElem(j_max_-2)-2.*getPsiElem(j_max_-1)+
                                     3/2.*getPsiElem(j_max_))/h_;}
        else {return (-1/2.*getPsiElem(j-1)+1/2.*getPsiElem(j+1))/h_;}
      break;
    case 2:
         if(j == 0) {return (2.*getPsiElem(j)-5.*getPsiElem(j+1)+4.*getPsiElem(j+2)-1.*getPsiElem(j+3))/(h_*h_);}
         else if(j == j_max_ ){return (-1.*getPsiElem(j-3)+4.*getPsiElem(j-2)-5.*getPsiElem(j-1)+2.*getPsiElem(j))/(h_*h_);}
         else {return (getPsiElem(j-1)-2.*getPsiElem(j)+getPsiElem(j+1))/(h_*h_);}      
      break;
    case 3:
       if(j == 0 || j == 1) {return (-5/2.*getPsiElem(j)+9.*getPsiElem(j+1)-12.*getPsiElem(j+2)+7.*getPsiElem(j+3)-3./2.*getPsiElem(j+4))/(h_*h_*h_);}
      else if (j == j_max_ || j == (j_max_-1)) {return (3./2.*getPsiElem(j-4)-7.*getPsiElem(j-3)+12.*getPsiElem(j-2)-9.*getPsiElem(j-1)+5./2.*getPsiElem(j))/(h_*h_*h_);}
      else {return (-1./2.*getPsiElem(j-2)+getPsiElem(j-1)-getPsiElem(j+1)+1./2.*getPsiElem(j+2))/(h_*h_*h_);}
      
      break;
    case 4:
      if(j == 0 || j == 1) {return (3.*getPsiElem(j)-14.*getPsiElem(j+1)+26.*getPsiElem(j+2)-24.*getPsiElem(j+3)+11.*getPsiElem(j+4)-2.*getPsiElem(j+5))/(h_*h_*h_*h_);}
      else if (j == j_max_ || j == (j_max_-1)){return (-2.*getPsiElem(j-5)+11.*getPsiElem(j-4)-24.*getPsiElem(j-3)+26.*getPsiElem(j-2)-14.*getPsiElem(j-1)+3.*getPsiElem(j))/(h_*h_*h_*h_);}
      else {return (getPsiElem(j-2)-4.*getPsiElem(j-1)+6.*getPsiElem(j)-4.*getPsiElem(j+1)+getPsiElem(j+2))/(h_*h_*h_*h_);}
      break;
    default:
      return 10000000;
      break;
  }

}

void Psi::setPsiElem(int j, long double wartoscPsi) {
  v_psi_[j] = wartoscPsi;
}