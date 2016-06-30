#include <iostream>
#include <cmath>
#include <armadillo>

#include "PsiNJ.h"
#include "PolynomialCurveFitting.h"

namespace alg {

PolynomialCurveFitting::PolynomialCurveFitting(int jmax, long double rmin, long double h, int min_rzad_poly, int max_rzad_poly, Psi y) 
    : x_(jmax, h), y_(y), A_(max_rzad_poly+1,max_rzad_poly+1), h_(h), r_min_(rmin), wynik_(jmax, h), v_wspocz_wielomianu_(max_rzad_poly+1) {
      
  X_.resize(max_rzad_poly+1);
  max_rzad_polynomial_ = max_rzad_poly;
  min_rzad_polynomial_ = min_rzad_poly;
  j_max_ = jmax;
  
  Psi x(jmax,h);
  for(int j = 0; j <= jmax; ++j){
    x.setPsiElem(j,rmin + h*j);
  }
  
  x_ = x;
}

long double PolynomialCurveFitting::SumPowX(int potega) {

  long double tymcz_double = 0;
  
  for(int i = 0; i<=j_max_; ++i) {
    tymcz_double = tymcz_double + (1./pow(x_.getPsiElem(i), potega));  
  }
  
  return tymcz_double;
}

long double PolynomialCurveFitting::SumPowXY(int potega){
  long double tymcz_double = 0;
  
  for(int i = 0; i<=j_max_; ++i) {
    tymcz_double = tymcz_double + y_.getPsiElem(i)*(1./pow(x_.getPsiElem(i), potega));
  }
  
  return tymcz_double;

}

void PolynomialCurveFitting::constructMatrixA() {

  arma::mat A(max_rzad_polynomial_+1, max_rzad_polynomial_+1);
  A.zeros();
   
  for (int j = min_rzad_polynomial_; j <= max_rzad_polynomial_; ++j) {
    for (int k = min_rzad_polynomial_; k <= max_rzad_polynomial_; ++k) {
      A(j,k) = SumPowX(k+j);
    }
  }
   
   A_ = A;
}

void PolynomialCurveFitting::constructMatrixB() {
  arma::vec B(max_rzad_polynomial_+1);
  B.zeros();
  
  for (int j = min_rzad_polynomial_; j <= max_rzad_polynomial_; ++j) {
    B(j) = SumPowXY(j); 
  }
   
  B_ = B;
}

Psi PolynomialCurveFitting::getWynik(){
  return wynik_;
}


void PolynomialCurveFitting::Execute() {
  constructMatrixA();
  constructMatrixB();
  
  arma::vec x = solve( A_, B_ );
  
  Psi wynik(j_max_, h_);
  
  for(int j = 0; j <= j_max_; ++j){
    long double liczba = 0;
  
    for(int k = min_rzad_polynomial_; k <= max_rzad_polynomial_; ++k){
      liczba = liczba + x(k)/pow(r_min_ + h_*j, k);
    }
    
    wynik.setPsiElem(j,liczba);
  }
      
  wynik_ = wynik;
      
  for(int j = min_rzad_polynomial_; j <= max_rzad_polynomial_; ++j){
    v_wspocz_wielomianu_[j] = x(j);
  }

}

}