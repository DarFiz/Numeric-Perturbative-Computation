#pragma once
#include <vector>
#include <PsiNJ.h>

namespace alg {
  
class ITaylorCoeff {
public:

  virtual ~ITaylorCoeff(){};

  virtual void computePierwszyRzad(int jmax, long double eps) = 0;
  virtual void computeDrugiRzad(int jmax, long double eps) = 0;
  
  virtual long double convFactorM2(int jmax, long double eps) = 0;
  
  virtual void plotM1vsM1ExactoPs(long double h, long double r_gorny_zakres, 
                                  long double eps) = 0;
  virtual void plotPhi1vsPhi1ExactoPs(long double h, long double r_gorny_zakres, 
                                  long double eps) = 0;
  virtual void plotM2vsM2ExactoPs(long double h, long double r_gorny_zakres, 
                          long double eps) = 0;
  virtual void plotPhi2vsPhi2ExactoPs(long double h, long double r_gorny_zakres, 
                          long double eps) = 0;
  virtual void plotErrorM1ExactVsExactDNToPs(long double hdol, long double hgora, long double hstep, long double eps, int rzad_poch) = 0;
  
  virtual void plotErrorM1DNToPs(long double hdol, long double hgora, long double hstep, long double eps, int rzad_poch) = 0;
  virtual void plotErrorPhi1DNToPs(long double hdol, long double hgora, long double hstep, long double eps, int rzad_poch) = 0;
  
  virtual void plotRelativeErrorM1ToPs(long double hdol, long double hgora, long double hstep, long double eps) = 0;
  
  virtual void plotErrorM2ToPs(long double hdol, long double hgora, long double hstep, long double eps) = 0;
  virtual void plotErrorPhi2ToPs(long double hdol, long double hgora, long double hstep, long double eps) = 0;
  
  virtual void plotRelativeErrorM2ToPs(long double hdol, long double hgora, long double hstep, long double eps) = 0;
  
  virtual void plotConvFactorM2ToPs(long double hdol, long double hgora, long double hstep, long double eps) = 0;
  
  virtual void plotTwoPsiToPs(long double rmin, long double rmax, long double r_gorne, Psi psi1, Psi psi2, std::string nazwa_pliku_ps) = 0;

};

}
