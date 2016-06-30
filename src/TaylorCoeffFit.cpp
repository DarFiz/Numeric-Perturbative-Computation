#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "PsiNJ.h"
#include "TaylorCoeffFit.h"
#include "TaylorCoeffExact.h"
#include "PolynomialCurveFitting.h"
#include "IntegralNum.h"
#include "TabFunc.h"
#include "Error.h"
#include "ConvFactor.h"
#include "gnuplot_only_declaration.h"


namespace alg {

TaylorCoeffFit::TaylorCoeffFit(long double rmin, long double 
                                           rmax, long double horizon)
    :r_min_(rmin), r_max_(rmax), r_horizon_(horizon){

}


std::vector<long double> TaylorCoeffFit::getWspWielomianuM(int rzadRozw) {
  return wsp_wielomianu_rozwM_[rzadRozw];
}

std::vector<long double> TaylorCoeffFit::getWspWielomianuPhi(int rzadRozw) {
  return wsp_wielomianu_rozwPhi_[rzadRozw];
}


void TaylorCoeffFit::computePierwszyRzad(int jmax, long double eps){
  
  TaylorCoeffExact exact_taylor_coeff(r_min_, r_max_, r_horizon_);
  Psi m_zero_order1h = exact_taylor_coeff.computeM0(jmax);
  Psi phi_zero_order1h = exact_taylor_coeff.computePhi0(jmax); 
  
  alg::IntegralNum dob_egz1h(r_min_, r_max_, jmax);
  dob_egz1h.setZeroOrder(m_zero_order1h, phi_zero_order1h);
  dob_egz1h.setHorizonRadius(r_horizon_);
  dob_egz1h.executeIntegralInfinityZero(1, eps);
  Psi PhiFromZeroOrder = dob_egz1h.getPhinPlus1();
  Psi MFromZeroOrder = dob_egz1h.getMnPlus1();
  
  long double h = (r_max_-r_min_)/jmax;
  
  alg::PolynomialCurveFitting  curve_fitM1(jmax, r_min_, h, 5, 6, (m_zero_order1h*(-1.) +
                                                                   MFromZeroOrder)*pow(1/eps,1));
  curve_fitM1.Execute();
  
  std::vector<long double> wsp_wiel_M0(1);
  std::vector<long double> wsp_wiel_M1(7);
  
  
  wsp_wiel_M1[5] = curve_fitM1.v_wspocz_wielomianu_[5];
  wsp_wiel_M1[6] = curve_fitM1.v_wspocz_wielomianu_[6];
  
  if (wsp_wielomianu_rozwM_.size() == 0) {
    wsp_wielomianu_rozwM_.push_back(wsp_wiel_M0);
    wsp_wielomianu_rozwM_.push_back(wsp_wiel_M1);
  } 
  else {
    wsp_wielomianu_rozwM_[0] = wsp_wiel_M0;
    wsp_wielomianu_rozwM_[1] = wsp_wiel_M1;
  }
  
  
  alg::PolynomialCurveFitting  curve_fitPhi1(jmax, r_min_, h, 6, 6, (phi_zero_order1h*(-1.) +
                                                                     PhiFromZeroOrder)*pow(1/eps,1));
  curve_fitPhi1.Execute();
  std::vector<long double> wsp_wiel_Phi0(1);
  std::vector<long double> wsp_wiel_Phi1(7);
  wsp_wiel_Phi1[6] = curve_fitPhi1.v_wspocz_wielomianu_[6];
  
  if (wsp_wielomianu_rozwPhi_.size() == 0) {
    wsp_wielomianu_rozwPhi_.push_back(wsp_wiel_Phi0);
    wsp_wielomianu_rozwPhi_.push_back(wsp_wiel_Phi1);
  } 
  else {
    wsp_wielomianu_rozwPhi_[0] = wsp_wiel_Phi0;
    wsp_wielomianu_rozwPhi_[1] = wsp_wiel_Phi1;
  }
}




Psi TaylorCoeffFit::getM1asPsi(int jmax, long double eps){
  
  long double h = (r_max_-r_min_)/jmax;
  
  Psi wynik(jmax, h);
  for(int j = 0; j <= jmax; ++j){
    long double liczba = 0;
    liczba = liczba + wsp_wielomianu_rozwM_[1][5]/pow(r_min_ + h*j, 5);
    liczba = liczba + wsp_wielomianu_rozwM_[1][6]/pow(r_min_ + h*j, 6);
    wynik.setPsiElem(j,liczba);
  }
  
  return wynik;
}





Psi TaylorCoeffFit::getPhi1asPsi(int jmax, long double eps){
  
  long double h = (r_max_-r_min_)/jmax;
  
  Psi wynik(jmax, h);
  for(int j = 0; j <= jmax; ++j){
    long double liczba = 0;
    liczba = liczba + wsp_wielomianu_rozwPhi_[1][6]/pow(r_min_ + h*j, 6);
    wynik.setPsiElem(j,liczba);
  }
  
  return wynik;
}


void TaylorCoeffFit::computeDrugiRzad(int jmax, long double eps) {

  TaylorCoeffExact exact_taylor_coeff(r_min_, r_max_, r_horizon_);
  Psi zerowy_rzad_M = exact_taylor_coeff.computeM0(jmax);
  Psi zerowy_rzad_Phi = exact_taylor_coeff.computePhi0(jmax);
  
  Psi pierwszy_rzad_M_eps1 = getM1asPsi(jmax, eps);
  Psi pierwszy_rzad_Phi_eps1 = getPhi1asPsi(jmax, eps);


  alg::IntegralNum dob_M2_eps1(r_min_, r_max_, jmax);
  dob_M2_eps1.setZeroOrder(zerowy_rzad_M + pierwszy_rzad_M_eps1*pow(eps,1), zerowy_rzad_Phi + pierwszy_rzad_Phi_eps1*pow(eps,1));
  dob_M2_eps1.setHorizonRadius(r_horizon_);
  dob_M2_eps1.executeIntegralInfinityZero(1, eps);


  Psi pierwszy_rzad_M_eps2 = getM1asPsi(jmax, 2.*eps);
  Psi pierwszy_rzad_Phi_eps2 = getPhi1asPsi(jmax, 2.*eps);
  
  alg::IntegralNum dob_M2_eps2(r_min_, r_max_, jmax);
  dob_M2_eps2.setZeroOrder(zerowy_rzad_M + pierwszy_rzad_M_eps2*pow(2.*eps,1), zerowy_rzad_Phi + pierwszy_rzad_Phi_eps2*pow(2.*eps,1));
  dob_M2_eps2.setHorizonRadius(r_horizon_);
  dob_M2_eps2.executeIntegralInfinityZero(1, 2.*eps);
  
  
  Psi pierwszy_rzad_M_eps3 = getM1asPsi(jmax, 3.*eps);
  Psi pierwszy_rzad_Phi_eps3 = getPhi1asPsi(jmax, 3.*eps);
  
  alg::IntegralNum dob_M2_eps3(r_min_, r_max_, jmax);
  dob_M2_eps3.setZeroOrder(zerowy_rzad_M + pierwszy_rzad_M_eps3*pow(3.*eps,1), zerowy_rzad_Phi + pierwszy_rzad_Phi_eps3*pow(3.*eps,1));
  dob_M2_eps3.setHorizonRadius(r_horizon_);
  dob_M2_eps3.executeIntegralInfinityZero(1, 3.*eps);
  
  
  
  long double h = (r_max_-r_min_)/jmax;

  alg::PolynomialCurveFitting  curve_fitM2(jmax, r_min_, h, 11, 
                                           12, (zerowy_rzad_M*(2.) + 
                                                dob_M2_eps1.getMnPlus1()*(-5.) +
                                                dob_M2_eps2.getMnPlus1()*(4.) +
                                                dob_M2_eps3.getMnPlus1()*(-1.))*pow(1/eps,2)*(1./2.));
  curve_fitM2.Execute();
  
  std::vector<long double> wsp_wiel_M2(13);
  wsp_wiel_M2[11] = curve_fitM2.v_wspocz_wielomianu_[11];
  wsp_wiel_M2[12] = curve_fitM2.v_wspocz_wielomianu_[12];
  
  if(wsp_wielomianu_rozwM_.size() == 2){
    wsp_wielomianu_rozwM_.push_back(wsp_wiel_M2);
    
  } 
  else if (wsp_wielomianu_rozwM_.size() == 3) {
    wsp_wielomianu_rozwM_[2] = wsp_wiel_M2;
  }
  
  
  alg::PolynomialCurveFitting  curve_fitPhi2(jmax, r_min_, h, 11, 
                                           12, (zerowy_rzad_Phi*(2.) + 
                                                dob_M2_eps1.getPhinPlus1()*(-5.) +
                                                dob_M2_eps2.getPhinPlus1()*(4.) +
                                                dob_M2_eps3.getPhinPlus1()*(-1.))*pow(1/eps,2)*(1./2.));
  curve_fitPhi2.Execute();
  
  std::vector<long double> wsp_wiel_Phi2(13);
  wsp_wiel_Phi2[11] = curve_fitPhi2.v_wspocz_wielomianu_[11];
  wsp_wiel_Phi2[12] = curve_fitPhi2.v_wspocz_wielomianu_[12];
  
  if(wsp_wielomianu_rozwPhi_.size() == 2){
    wsp_wielomianu_rozwPhi_.push_back(wsp_wiel_Phi2);
    
  } 
  else if (wsp_wielomianu_rozwPhi_.size() == 3) {
    wsp_wielomianu_rozwPhi_[2] = wsp_wiel_Phi2;
  }
  
  
}






Psi TaylorCoeffFit::getM2asPsi(int jmax, long double eps){
  
  long double h = (r_max_-r_min_)/jmax;
  
  Psi wynik(jmax, h);
  for(int j = 0; j <= jmax; ++j){
    long double liczba = 0;
    liczba = liczba + wsp_wielomianu_rozwM_[2][11]/pow(r_min_ + h*j, 11);
    liczba = liczba + wsp_wielomianu_rozwM_[2][12]/pow(r_min_ + h*j, 12);
    wynik.setPsiElem(j,liczba);
  }
  
  return wynik;
}

Psi TaylorCoeffFit::getPhi2asPsi(int jmax, long double eps){
  
  long double h = (r_max_-r_min_)/jmax;
  
  Psi wynik(jmax, h);
  for(int j = 0; j <= jmax; ++j){
    long double liczba = 0;
    liczba = liczba + wsp_wielomianu_rozwPhi_[2][11]/pow(r_min_ + h*j, 11);
    liczba = liczba + wsp_wielomianu_rozwPhi_[2][12]/pow(r_min_ + h*j, 12);
    
    wynik.setPsiElem(j,liczba);
  }
  
  return wynik;
}

long double TaylorCoeffFit::convFactorM1DN(int jmax, long double eps, int rzad_poch)
{
  computePierwszyRzad(jmax * 1, eps);
  Psi m_first_order4h = getM1asPsi(jmax, eps);

  computePierwszyRzad(jmax * 2, eps);
  Psi m_first_order2h = getM1asPsi(jmax * 2, eps);
  
  computePierwszyRzad(jmax * 4, eps);
  Psi m_first_order1h = getM1asPsi(jmax * 4, eps);

  ConvFactor conv_factor;
  conv_factor.setPsih(m_first_order1h);
  conv_factor.setPsi2h(m_first_order2h);
  conv_factor.setPsi4h(m_first_order4h);
  
  if (rzad_poch == 0){
    return conv_factor.computeConvFactor();}
  else {
    return conv_factor.computeConvFactorD(rzad_poch);}
}

long double TaylorCoeffFit::convFactorM2(int jmax, long double eps) {
  computePierwszyRzad(jmax * 1., eps);
  computeDrugiRzad(jmax * 1., eps);
  Psi m_second_order4h = getM2asPsi(jmax, eps);
  
  computePierwszyRzad(jmax * 2., eps);
  computeDrugiRzad(jmax * 2., eps);
  Psi m_second_order2h = getM2asPsi(jmax*2, eps);
  
  computePierwszyRzad(jmax * 4., eps);
  computeDrugiRzad(jmax * 4., eps);
  Psi m_second_order1h = getM2asPsi(jmax*4, eps);

  ConvFactor conv_factor;
  conv_factor.setPsih(m_second_order1h);
  conv_factor.setPsi2h(m_second_order2h);
  conv_factor.setPsi4h(m_second_order4h);
  
  return conv_factor.computeConvFactor();
}


long double TaylorCoeffFit::convFactorPhi1DN(int jmax, long double eps, int  rzad_poch) {

  computePierwszyRzad(jmax / 4., eps);
  Psi phi_first_order4h = getPhi1asPsi(jmax / 4., eps);
  
  computePierwszyRzad(jmax / 2., eps);
  Psi phi_first_order2h = getPhi1asPsi(jmax / 2., eps);
  
  computePierwszyRzad(jmax / 1., eps);
  Psi phi_first_order1h = getPhi1asPsi(jmax / 1., eps);

  ConvFactor conv_factor;
  conv_factor.setPsih(phi_first_order1h);
  conv_factor.setPsi2h(phi_first_order2h);
  conv_factor.setPsi4h(phi_first_order4h);
  
  if (rzad_poch == 0){
    return conv_factor.computeConvFactor();}
  else {
    return conv_factor.computeConvFactorD(rzad_poch);}
}

void TaylorCoeffFit::plotM1vsM1ExactoPs(long double h, long double r_gorny_zakres, long double eps) {
  
  std::cout << "plotM1vsM1ExactoPs" << std::endl;
  
  int jmax = (r_max_ - r_min_) / h;
  
  TaylorCoeffExact exact(r_min_, r_max_, r_horizon_);
  computePierwszyRzad(jmax, eps);
  
  plotTwoPsiToPs(r_min_, r_max_, r_gorny_zakres, getM1asPsi(jmax, eps), exact.computeM1DN(jmax, eps, 0), 
                 "M1 fit h="+std::to_string(h)+" zak="+std::to_string(r_gorny_zakres));
}

void TaylorCoeffFit::plotPhi1vsPhi1ExactoPs(long double h, long double r_gorny_zakres, long double eps) {
  
  std::cout << "plotPhi1vsPhi1ExactoPs" << std::endl;
  
  int jmax = (r_max_ - r_min_) / h;
  
  TaylorCoeffExact exact(r_min_, r_max_, r_horizon_);
  computePierwszyRzad(jmax, eps);
  
  plotTwoPsiToPs(r_min_, r_max_, r_gorny_zakres, getPhi1asPsi(jmax, eps), exact.computePhi1DN(jmax, eps, 0), 
                 "Phi1 fit h="+std::to_string(h)+" zak="+std::to_string(r_gorny_zakres));
}


void TaylorCoeffFit::plotM2vsM2ExactoPs(long double h, long double r_gorny_zakres, long double eps) {
  
  std::cout << "plotM2vsM2ExactoPs" << std::endl;
  
  int jmax = (r_max_ - r_min_) / h;
  
  TaylorCoeffExact exact(r_min_, r_max_, r_horizon_);
  computePierwszyRzad(jmax, eps);
  computeDrugiRzad(jmax, eps);
  
  plotTwoPsiToPs(r_min_, r_max_, r_gorny_zakres, getM2asPsi(jmax, eps), exact.computeM2(jmax, eps), 
                 "M2 fit h="+std::to_string(h)+" zak="+std::to_string(r_gorny_zakres));
}

void TaylorCoeffFit::plotPhi2vsPhi2ExactoPs(long double h, long double r_gorny_zakres, long double eps) {
  
  std::cout << "plotPhi2vsPhi2ExactoPs" << std::endl;
  
int jmax = (r_max_ - r_min_) / h;
  
  TaylorCoeffExact exact(r_min_, r_max_, r_horizon_);
  computePierwszyRzad(jmax, eps);
  computeDrugiRzad(jmax, eps);
  
  plotTwoPsiToPs(r_min_, r_max_, r_gorny_zakres, getPhi2asPsi(jmax, eps), exact.computePhi2(jmax, eps), 
                 "Phi2 fit h="+std::to_string(h)+" zak="+std::to_string(r_gorny_zakres));
}


void TaylorCoeffFit::plotErrorM1ExactVsExactDNToPs(long double hdol, long double hgora, long double hstep, long double eps, int rzad_poch) {
        Gnuplot g1("lines");
     
        std::vector<long double> error_y, h_x;
        TaylorCoeffExact exact(r_min_, r_max_, r_horizon_);
        for (long double h = hdol; h < hgora; h += hstep) {
          int i = (r_max_ - r_min_) / h;
          h_x.push_back(h);
          if (rzad_poch == 1) {error_y.push_back(Error::computeErrorD( 1,exact.computeM1DN(i, eps, 0), exact.computeM1DN(i, eps, 1)));}
          if (rzad_poch == 2) {error_y.push_back(Error::computeErrorD( 2,exact.computeM1DN(i, eps, 0), exact.computeM1DN(i, eps, 2)));}
          if (rzad_poch == 3) {error_y.push_back(Error::computeErrorD( 3,exact.computeM1DN(i, eps, 0), exact.computeM1DN(i, eps, 3)));}
          if (rzad_poch == 4) {error_y.push_back(Error::computeErrorD( 4,exact.computeM1DN(i, eps, 0), exact.computeM1DN(i, eps, 4)));}
          
          std::cout <<"h="+std::to_string(h) <<std::endl;
        }

        g1.reset_all();
        g1.savetops("PlotErrorM1ExactD"+ std::to_string(rzad_poch));
        g1.reset_plot();
        g1.set_grid();
        g1.set_style("linespoints").plot_xy(h_x, error_y,"Error M1ExactD"+std::to_string(rzad_poch)+" od h");
        g1.showonscreen();
}

void TaylorCoeffFit::plotErrorM1DNToPs(long double hdol, long double hgora, long double hstep, long double eps, int rzad_poch){
  
  std::cout << "plotErrorM1DNToPs" << std::endl;
  
Gnuplot g1("lines");
  std::vector<long double> error_y, h_x;
  TaylorCoeffExact exact(r_min_, r_max_, r_horizon_);
  
  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    h_x.push_back(h);
    computePierwszyRzad(i, eps);
    if (rzad_poch == 0) {error_y.push_back(Error::computeError(    getM1asPsi(i, eps), exact.computeM1DN(i, eps, 0)));}
    if (rzad_poch == 1) {error_y.push_back(Error::computeErrorD( 1,getM1asPsi(i, eps), exact.computeM1DN(i, eps, 1)));}
    if (rzad_poch == 2) {error_y.push_back(Error::computeErrorD( 2,getM1asPsi(i, eps), exact.computeM1DN(i, eps, 2)));}
    if (rzad_poch == 3) {error_y.push_back(Error::computeErrorD( 3,getM1asPsi(i, eps), exact.computeM1DN(i, eps, 3)));}
    if (rzad_poch == 4) {error_y.push_back(Error::computeErrorD( 4,getM1asPsi(i, eps), exact.computeM1DN(i, eps, 4)));}
    
    std::cout <<"conv factor" <<std::endl;
  }

  g1.reset_all();
  g1.savetops("Err M1D"+ std::to_string(rzad_poch)+" fit dh=" + std::to_string(hstep));
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("linespoints").plot_xy(h_x, error_y,"Error M1D"+std::to_string(rzad_poch) + " fit dh=" + std::to_string(hstep));
  g1.showonscreen();
}

void TaylorCoeffFit::plotErrorPhi1DNToPs(long double hdol, long double hgora, long double hstep, long double eps, int rzad_poch){
  
  std::cout << "plotErrorPhi1DNToPs" << std::endl;
  
  Gnuplot g1("lines");
  std::vector<long double> error_y, h_x;
  TaylorCoeffExact exact_taylor_coeff(r_min_, r_max_, r_horizon_);

  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    h_x.push_back(h);
    computePierwszyRzad(i, eps);
    if (rzad_poch == 0) {error_y.push_back(Error::computeError(    getPhi1asPsi(i, eps), exact_taylor_coeff.computePhi1DN(i, eps, 0)));}
    if (rzad_poch == 1) {error_y.push_back(Error::computeErrorD( 1,getPhi1asPsi(i, eps), exact_taylor_coeff.computePhi1DN(i, eps, 1)));}
    if (rzad_poch == 2) {error_y.push_back(Error::computeErrorD( 2,getPhi1asPsi(i, eps), exact_taylor_coeff.computePhi1DN(i, eps, 2)));}
    if (rzad_poch == 3) {error_y.push_back(Error::computeErrorD( 3,getPhi1asPsi(i, eps), exact_taylor_coeff.computePhi1DN(i, eps, 3)));}
    if (rzad_poch == 4) {error_y.push_back(Error::computeErrorD( 4,getPhi1asPsi(i, eps), exact_taylor_coeff.computePhi1DN(i, eps, 4)));}
    
    std::cout <<"conv factor" <<std::endl;
  }

  g1.reset_all();
  g1.savetops("Err Phi1D"+ std::to_string(rzad_poch)+ " fit dh=" + std::to_string(hstep));
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("linespoints").plot_xy(h_x, error_y,"Error Phi1D"+std::to_string(rzad_poch)+ " fit dh=" + std::to_string(hstep));
  g1.showonscreen();
}

void TaylorCoeffFit::plotRelativeErrorM1ToPs(long double hdol, long double hgora, long double hstep, long double eps){
  Gnuplot g1("lines");
  std::vector<long double> error_y, h_x;

  for (long double h = hdol; h < hgora; h += hstep) {
    
    std::cout <<"h="+std::to_string(h) <<std::endl;
    
    h_x.push_back(h);
    
    int i = (r_max_ - r_min_) / h;
    computePierwszyRzad(i, eps);
    Psi MnI(i, h);
    MnI = getM1asPsi(i, eps);
    
    computePierwszyRzad(2*i, eps);
    Psi Mn2I(2*i,h/2);
    Mn2I = getM1asPsi(i, eps);
    
    error_y.push_back(Error::computeErrorIvs2I( MnI, Mn2I));
    
  }

  g1.reset_all();
  g1.savetops("PlotRelativeErrorM1 curve fit");
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("linespoints").plot_xy(h_x, error_y,"Relative Error M1 od h");
  g1.showonscreen();
}

void TaylorCoeffFit::plotErrorM2ToPs(long double hdol, long double hgora, long double hstep, long double eps){
  
  std::cout << "plotErrorM2ToPs" << std::endl;
  
  Gnuplot g1("lines");
  std::vector<long double> error_y, h_x;

  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    
    std::cout <<"h="<<std::to_string(h) <<std::endl;
    TaylorCoeffExact exact_taylor_coeff(r_min_, r_max_, r_horizon_);
    
    h_x.push_back(h);
    computePierwszyRzad(i, eps);
    computeDrugiRzad(i, eps);
    error_y.push_back(Error::computeError(getM2asPsi(i, eps), exact_taylor_coeff.computeM2(i, eps)));
    
  }

  g1.reset_all();
  g1.savetops("Err M2 fit dh=" + std::to_string(hstep));
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("linespoints").plot_xy(h_x, error_y,"Error M2 fit dh=" + std::to_string(hstep));
  g1.showonscreen();
}

void TaylorCoeffFit::plotErrorPhi2ToPs(long double hdol, long double hgora, long double hstep, long double eps){
  
  std::cout << "plotErrorPhi2ToPs" << std::endl;
  
  Gnuplot g1("lines");
  std::vector<long double> error_y, h_x;

  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    
    std::cout <<"h="<<std::to_string(h) <<std::endl;
    TaylorCoeffExact exact_taylor_coeff(r_min_, r_max_, r_horizon_);
    
    h_x.push_back(h);
    computePierwszyRzad(i, eps);
    computeDrugiRzad(i, eps);
    error_y.push_back(Error::computeError(getPhi2asPsi(i, eps), exact_taylor_coeff.computePhi2(i, eps)));
    
  }

  g1.reset_all();
  g1.savetops("Err Phi2 fit dh=" + std::to_string(hstep));
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("linespoints").plot_xy(h_x, error_y,"Error Phi2 fit dh=" + std::to_string(hstep));
  g1.showonscreen();

}

void TaylorCoeffFit::plotRelativeErrorM2ToPs(long double hdol, long double hgora, long double hstep, long double eps){
  Gnuplot g1("lines");
  std::vector<long double> error_y, h_x;

  for (long double h = hdol; h < hgora; h += hstep) {
    
    std::cout <<"h="+std::to_string(h) <<std::endl;
    
    h_x.push_back(h);
    
    int i = (r_max_ - r_min_) / h;
    computePierwszyRzad(i, eps);
    computeDrugiRzad(i, eps);
    Psi MnI = getM2asPsi(i, eps);
    
    computePierwszyRzad(2*i, eps);
    computeDrugiRzad(2*i, eps);
    Psi Mn2I = getM2asPsi(i, eps);
    
    error_y.push_back(Error::computeErrorIvs2I( MnI, Mn2I));
    
  }

  g1.reset_all();
  g1.savetops("PlotRelativeErrorM2 Curve Fit");
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("linespoints").plot_xy(h_x, error_y,"Relative Error M2 od h");
  g1.showonscreen();
}

void TaylorCoeffFit::plotConvFactorM1DNToPs(long double hdol, long double hgora, 
                                                 long double hstep, long double eps, 
                                                 int rzad_poch){
        Gnuplot g1("lines");
     
        std::vector<long double> error_y, h_x;

        for (long double h = hdol; h < hgora; h += hstep) {
          int i = (r_max_ - r_min_) / h;
          h_x.push_back(h);
          long double convFactor = convFactorM1DN(i, eps, rzad_poch);
          error_y.push_back(convFactor);
          
          std::cout <<"conv factor" <<std::endl;
        }

        g1.reset_all();
        g1.savetops("PlotConvergenceFactorM1CurveFit");
        g1.reset_plot();
        g1.set_grid();
        g1.set_style("linespoints").plot_xy(h_x, error_y,"convFactorM1DN od h");
        g1.showonscreen();
}

void TaylorCoeffFit::plotConvFactorPhi1DNToPs(long double hdol, long double hgora, 
                                                   long double hstep, long double eps, 
                                                   int rzad_poch){
        Gnuplot g1("lines");
     
        std::vector<long double> error_y, h_x;

        for (long double h = hdol; h < hgora; h += hstep) {
          int i = (r_max_ - r_min_) / h;
          h_x.push_back(h);
          long double convFactor = convFactorPhi1DN(i, eps, rzad_poch);
          error_y.push_back(convFactor);
          
          std::cout <<"conv factor" <<std::endl;
        }

        g1.reset_all();
        g1.savetops("PlotConvergenceFactorPhi1CurveFit");
        g1.reset_plot();
        g1.set_grid();
        g1.set_style("linespoints").plot_xy(h_x, error_y,"convFactorPhi1DN od h");
        g1.showonscreen();
}

void TaylorCoeffFit::plotConvFactorM2ToPs(long double hdol, long double hgora, long double hstep, long double eps){
Gnuplot g1("lines");
     
        std::vector<long double> error_y, h_x;

        for (long double h = hdol; h <= hgora; h += hstep) {
          int i = (r_max_ - r_min_) / h;
          h_x.push_back(h);
          std::cout <<"h="<<std::to_string(h) <<std::endl;
          long double convFactor = convFactorM2(i, eps);
          error_y.push_back(convFactor);
        }

        g1.reset_all();
        g1.savetops("PlotConvergenceFactorM2CurveFit");
        g1.reset_plot();
        g1.set_grid();
        g1.set_style("linespoints").plot_xy(h_x, error_y,"convFactorM2 od h");
        g1.showonscreen();
}

void TaylorCoeffFit::plotTwoPsiToPs(long double rmin, long double rmax, long double r_gorne, Psi psi_przybliz, Psi psi_dokl, std::string nazwa_pliku_ps){
  Gnuplot g1("lines");

  std::vector<double> y_dokl, y_przybliz, x;
  int jmax1 = psi_przybliz.getMaxj();
  long double h=(rmax-rmin)/jmax1;

  int jgora = ((r_gorne - rmin) / (rmax - rmin))*jmax1;

  for (int i = 0; i < jgora; i++)
  {
      x.push_back(rmin+h*i);
      y_przybliz.push_back(psi_przybliz.getPsiElem(i));
      y_dokl.push_back(psi_dokl.getPsiElem(i));
  }

  g1.reset_all();
  g1.savetops(nazwa_pliku_ps);
  g1.reset_plot();
  std::cout << std::endl << std::endl << "*** user-defined lists of points (x,y)" << std::endl;
  g1.set_grid();
  g1.set_style("points").plot_xy(x, y_przybliz,"Przyblizony " + nazwa_pliku_ps);
  g1.set_style("lines").plot_xy(x, y_dokl,"Dokladny " + nazwa_pliku_ps);
  g1.showonscreen();
}

void TaylorCoeffFit::plotM1AnCoeff(long double hdol, long double hgora, long double hstep, long double eps,long double wsp_dokladny, int numer_wsp){
  std::cout << "plotM1AnCoeff" <<std::endl;
  Gnuplot g1("lines");
  std::vector<long double> y_przybliz,y_dokladny, h_x;

  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    h_x.push_back(1/h);
    computePierwszyRzad(i, eps);
    
    std::cout << "h=" << h << std::endl;
    std::cout << getWspWielomianuM(1)[numer_wsp] <<std::endl;
    
    y_przybliz.push_back(getWspWielomianuM(1)[numer_wsp]);
    y_dokladny.push_back(wsp_dokladny);
  }

  g1.reset_all();
  g1.savetops("M1 fit Wsp "+std::to_string(numer_wsp) + " dh=" + std::to_string(hstep));

  g1.reset_plot();
  g1.set_grid();
//   g1.set_yrange(-108.02, -107.998);
  g1.set_style("points").plot_xy(h_x, y_przybliz,"M1 fit Przybliz Wsp a "+std::to_string(numer_wsp)+" od 1/h dh=" + std::to_string(hstep));
  g1.set_style("lines").plot_xy(h_x, y_dokladny,"M1 fit Dokladny Wsp a "+std::to_string(numer_wsp)+" od 1/h dh=" + std::to_string(hstep));
  g1.showonscreen();
}

void TaylorCoeffFit::plotPhi1AnCoeff(long double hdol, long double hgora, long double hstep, long double eps, long double wsp_dokladny, int numer_wsp){
  
  std::cout << "plotPhi1AnCoeff" <<std::endl;
  
  Gnuplot g1("lines");
  std::vector<long double> y_przybliz,y_dokladny, h_x;

  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    h_x.push_back(1/h);
    computePierwszyRzad(i, eps);
    
    std::cout << "h=" << h << std::endl;
    std::cout << getWspWielomianuPhi(1)[numer_wsp] <<std::endl;
    
    y_przybliz.push_back(getWspWielomianuPhi(1)[numer_wsp]);
    y_dokladny.push_back(wsp_dokladny);
  }

  g1.reset_all();
  g1.savetops("Phi1 fit Wsp b"+std::to_string(numer_wsp)+ " dh=" + std::to_string(hstep));
  g1.reset_plot();
  g1.set_grid();
//   g1.set_yrange(-108.02, -107.998);
  g1.set_style("points").plot_xy(h_x, y_przybliz,"Phi1 fit Przybliz Wsp b "+std::to_string(numer_wsp) + " od 1/h dh=" + std::to_string(hstep));
  g1.set_style("lines").plot_xy(h_x, y_dokladny,"Phi1 fit Dokladny Wsp b "+std::to_string(numer_wsp) + " od 1/h dh=" + std::to_string(hstep));
  g1.showonscreen();
}

void TaylorCoeffFit::plotM1AnCoeffError(long double hdol, long double hgora, long double hstep, long double eps, long double wsp_dokladny, int numer_wsp){

Gnuplot g1("lines");
  std::vector<long double> y_przybliz, h_x;

  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    h_x.push_back(1/h);
    computePierwszyRzad(i, eps);
    y_przybliz.push_back(sqrt((getWspWielomianuM(1)[numer_wsp]-wsp_dokladny)*
                              (getWspWielomianuM(1)[numer_wsp]-wsp_dokladny)));
  }

  g1.reset_all();
  g1.savetops("Err M1 fit a"+std::to_string(numer_wsp));
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("points").plot_xy(h_x, y_przybliz,"Error M1 fit a "+std::to_string(numer_wsp)+" od 1/h");
  g1.showonscreen();
}

void TaylorCoeffFit::plotPhi1AnCoeffError(long double hdol, long double hgora, long double hstep, long double eps, long double wsp_dokladny, int numer_wsp){
  Gnuplot g1("lines");
  std::vector<long double> y_przybliz, h_x;

  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    h_x.push_back(1/h);
    computePierwszyRzad(i, eps);
    y_przybliz.push_back(sqrt((getWspWielomianuPhi(1)[numer_wsp]-wsp_dokladny)*
                              (getWspWielomianuPhi(1)[numer_wsp]-wsp_dokladny)));
  }

  g1.reset_all();
  g1.savetops("Err Phi1 fit a"+std::to_string(numer_wsp));
  g1.reset_plot();
  g1.set_grid();
//   g1.set_yrange(-108.02, -107.998);
  g1.set_style("points").plot_xy(h_x, y_przybliz,"Error Phi1 fit a "+std::to_string(numer_wsp)+" od 1/h");
  g1.showonscreen();
}

void TaylorCoeffFit::plotM2AnCoeff(long double hdol, long double hgora, long double hstep, long double eps, long double wsp_dokladny, int numer_wsp){
  
  std::cout << "plotM2AnCoeff" << std::endl;
  
  Gnuplot g1("lines");
  std::vector<long double> y_przybliz,y_dokladny, h_x;

  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    h_x.push_back(1/h);
    computePierwszyRzad(i, eps);
    computeDrugiRzad(i, eps);
    
    std::cout << "h=" << h << std::endl;
    std::cout << getWspWielomianuM(2)[numer_wsp] <<std::endl;
    
    y_przybliz.push_back(getWspWielomianuM(2)[numer_wsp]);
    y_dokladny.push_back(wsp_dokladny);
  }

  g1.reset_all();
  g1.savetops("M2 fit Wsp a"+std::to_string(numer_wsp)+ " dh=" + std::to_string(hstep));
  g1.reset_plot();
  g1.set_grid();
//   g1.set_yrange(-108.02, -107.998);
  g1.set_style("points").plot_xy(h_x, y_przybliz,"M2 fit Przybliz Wsp a "+std::to_string(numer_wsp)+ " od 1/h dh=" + std::to_string(hstep));
  g1.set_style("lines").plot_xy(h_x, y_dokladny,"M2 fit Dokladny Wsp a "+std::to_string(numer_wsp)+ " od 1/h dh=" + std::to_string(hstep));
  g1.showonscreen();
}

void TaylorCoeffFit::plotPhi2AnCoeff(long double hdol, long double hgora, long double hstep, long double eps, long double wsp_dokladny, int numer_wsp){
  
  std::cout << "plotPhi2AnCoeff" << std::endl;
  
  Gnuplot g1("lines");
  std::vector<long double> y_przybliz,y_dokladny, h_x;

  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    h_x.push_back(1/h);
    computePierwszyRzad(i, eps);
    computeDrugiRzad(i, eps);
    
    std::cout << "h=" << h << std::endl;
    std::cout << getWspWielomianuPhi(2)[numer_wsp] <<std::endl;
    
    y_przybliz.push_back(getWspWielomianuPhi(2)[numer_wsp]);
    y_dokladny.push_back(wsp_dokladny);
  }

  g1.reset_all();
  g1.savetops("Phi2 fit Wsp b"+std::to_string(numer_wsp) + " dh=" + std::to_string(hstep));
  g1.reset_plot();
  g1.set_grid();
//   g1.set_yrange(-108.02, -107.998);
  g1.set_style("points").plot_xy(h_x, y_przybliz,"Phi2 fit Przybliz Wsp b "+std::to_string(numer_wsp)+ " od 1/h dh=" + std::to_string(hstep));
  g1.set_style("lines").plot_xy(h_x, y_dokladny,"Phi2 fit Dokladny Wsp b "+std::to_string(numer_wsp)+ " od 1/h dh=" + std::to_string(hstep));
  g1.showonscreen();
}


}