#pragma once
#include <vector>
#include <PsiNJ.h>
#include "ITaylorCoeff.h"

namespace alg {
  
class TaylorCoeffFit: public ITaylorCoeff {
public:
  TaylorCoeffFit (long double rmin, long double rmax, long double horizon);

  Psi getM1asPsi(int jmax, long double eps);
  Psi getPhi1asPsi(int jmax, long double eps);
  Psi getM2asPsi(int jmax, long double eps);
  Psi getPhi2asPsi(int jmax, long double eps);

  std::vector<long double> getWspWielomianuM(int rzadRozw);
  std::vector<long double> getWspWielomianuPhi(int rzadRozw);

  void computePierwszyRzad(int jmax, long double eps);
  void computeDrugiRzad(int jmax, long double eps);
  
  void plotM1vsM1ExactoPs(long double h, long double r_gorny_zakres, 
                          long double eps);
  void plotPhi1vsPhi1ExactoPs(long double h, long double r_gorny_zakres, 
                              long double eps);
  void plotM2vsM2ExactoPs(long double h, long double r_gorny_zakres, 
                          long double eps);
  void plotPhi2vsPhi2ExactoPs(long double h, long double r_gorny_zakres, 
                          long double eps);
  
  void plotErrorM1ExactVsExactDNToPs(long double hdol, long double hgora, 
                                     long double hstep, long double eps, int rzad_poch);
  
  void plotErrorM1DNToPs(long double hdol, long double hgora, long double hstep, 
                         long double eps, int rzad_poch);
  void plotErrorPhi1DNToPs(long double hdol, long double hgora, long double hstep, 
                           long double eps, int rzad_poch);
  void plotRelativeErrorM1ToPs(long double hdol, long double hgora, long double hstep, 
                               long double eps);
  void plotErrorM2ToPs(long double hdol, long double hgora, long double hstep, 
                       long double eps);
  void plotErrorPhi2ToPs(long double hdol, long double hgora, long double hstep, 
                         long double eps);
  void plotRelativeErrorM2ToPs(long double hdol, long double hgora, long double hstep, 
                               long double eps);
  
  
  long double convFactorM1DN(int jmax, long double eps, int rzad_poch);
  long double convFactorPhi1DN(int jmax, long double eps,int rzad_poch);
  long double convFactorM2(int jmax, long double eps);
  
  
  
  
  
  void plotConvFactorM1DNToPs(long double hdol, long double hgora, long double hstep, 
                                   long double eps, int rzad_poch);
  void plotConvFactorPhi1DNToPs(long double hdol, long double hgora,long double hstep, 
                                     long double eps, int rzad_poch);
  
  void plotConvFactorM2ToPs(long double hdol, long double hgora, long double hstep, long double eps);
  
  void plotTwoPsiToPs(long double rmin, long double rmax, long double r_gorne, Psi psi_przybliz, Psi psi_dokl, std::string nazwa_pliku_ps);
  
  void plotM1AnCoeff(long double hdol, long double hgora, long double hstep, long double eps,long double wsp_dokladny, int numer_wsp);
  void plotPhi1AnCoeff(long double hdol, long double hgora, long double hstep, long double eps,long double wsp_dokladny, int numer_wsp);
  
  void plotM1AnCoeffError(long double hdol, long double hgora, long double hstep, long double eps,long double wsp_dokladny, int numer_wsp);
  void plotPhi1AnCoeffError(long double hdol, long double hgora, long double hstep, long double eps,long double wsp_dokladny, int numer_wsp);

  void plotM2AnCoeff(long double hdol, long double hgora, long double hstep, long double eps,long double wsp_dokladny, int numer_wsp);
  void plotPhi2AnCoeff(long double hdol, long double hgora, long double hstep, long double eps,long double wsp_dokladny, int numer_wsp);
  
  
  std::vector<std::vector<long double>> wsp_wielomianu_rozwM_;
  std::vector<std::vector<long double>> wsp_wielomianu_rozwPhi_;
  
  
  long double r_min_;
  long double r_max_;

  long double r_horizon_;
  
};

}