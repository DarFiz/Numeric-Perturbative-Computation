#pragma once
#include <vector>
#include <PsiNJ.h>
#include "ITaylorCoeff.h"
#include "Error.h"

namespace alg {
  
class TaylorCoeffClassic: public ITaylorCoeff {
public:
  TaylorCoeffClassic(long double rmin, long double rmax, long double horizon);

  Psi getMNasPsi(int n);
  Psi getPhiNasPsi(int n);
  
  void computePierwszyRzad(int jmax, long double eps);
  void computeDrugiRzad(int jmax, long double eps);
  void computeTrzeciRzad(int jmax, long double eps);
  
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
  long double convFactorPhi1DN(int jmax, long double eps, int rzad_poch);
  long double convFactorM2(int jmax, long double eps);

  void plotConvFactorM1DNToPs(long double hdol, long double hgora, long double hstep, 
                              long double eps, int rzad_poch);
  void plotConvFactorPhi1DNToPs(long double hdol, long double hgora, long double hstep, 
                                long double eps, int rzad_poch);
  void plotConvFactorM2ToPs(long double hdol, long double hgora, long double hstep, 
                            long double eps);
  
  void plotTwoPsiToPs(long double rmin, long double rmax, long double r_gorne, Psi psi1, 
                      Psi psi2, std::string nazwa_pliku_ps);

  std::vector<Psi> wsp_wielomianu_rozwMasPsi_;
  std::vector<Psi> wsp_wielomianu_rozwPhiasPsi_;
  
  long double r_min_;
  long double r_max_;

  long double r_horizon_;
};
}