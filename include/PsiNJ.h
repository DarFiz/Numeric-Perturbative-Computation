#pragma once
#include <vector>
#include "IPsiNJ.h"

class Psi : public IPsi {
public:
  Psi(int maxj, long double h);
  void setPsiElem(int j, long double wartoscPsi);
  long double getPsiElem(int j);
  long double getDPsiElem(int j, int rzadpoch);
  int getMaxj();
  long double geth();
  Psi operator+(const Psi &psi1);
  Psi operator-(const Psi &psi1);
  Psi operator* (const long double float_number);
private:
  std::vector<long double> v_psi_;
  long double h_;
  int j_max_;
};