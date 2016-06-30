#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "PsiNJ.h"
#include "TaylorCoeffClassic.h"
#include "TaylorCoeffExact.h"
#include "PolynomialCurveFitting.h"
#include "IntegralNum.h"
#include "TabFunc.h"
#include "Error.h"
#include "ConvFactor.h"
#include "gnuplot_i.hpp"

namespace alg {

TaylorCoeffClassic::TaylorCoeffClassic(long double rmin, long double 
                                           rmax, long double horizon)
    :r_min_(rmin), r_max_(rmax), r_horizon_(horizon){

}




void TaylorCoeffClassic::computePierwszyRzad(int jmax, long double eps){

  TaylorCoeffExact exact_taylor_coeff(r_min_, r_max_, r_horizon_);
  Psi m_zero_order1h = exact_taylor_coeff.computeM0(jmax);
  Psi phi_zero_order1h = exact_taylor_coeff.computePhi0(jmax); 

  alg::IntegralNum dob_egz1h(r_min_, r_max_, jmax);
  dob_egz1h.setZeroOrder(m_zero_order1h, phi_zero_order1h);
  dob_egz1h.setHorizonRadius(r_horizon_);
  dob_egz1h.executeIntegralInfinityZero(1, eps);
  Psi PhiFromZeroOrder = dob_egz1h.getPhinPlus1();
  Psi MFromZeroOrder = dob_egz1h.getMnPlus1();
  
  
  Psi m_first_order1h = (m_zero_order1h*(-1.) + MFromZeroOrder)*pow(1/eps,1);
 
  
  if (wsp_wielomianu_rozwMasPsi_.size() == 0) {
    wsp_wielomianu_rozwMasPsi_.push_back(m_zero_order1h);
    wsp_wielomianu_rozwMasPsi_.push_back(m_first_order1h);
  } 
  else {
    wsp_wielomianu_rozwMasPsi_[0] = m_zero_order1h;
    wsp_wielomianu_rozwMasPsi_[1] = m_first_order1h;
  }
  
  
  Psi phi_first_order1h = (phi_zero_order1h*(-1.) + PhiFromZeroOrder)*pow(1/eps,1);
  


  
  if (wsp_wielomianu_rozwPhiasPsi_.size() == 0) {
    wsp_wielomianu_rozwPhiasPsi_.push_back(phi_zero_order1h);
    wsp_wielomianu_rozwPhiasPsi_.push_back(phi_first_order1h);
  } 
  else {
    wsp_wielomianu_rozwPhiasPsi_[0] = phi_zero_order1h;
    wsp_wielomianu_rozwPhiasPsi_[1] = phi_first_order1h;
  }
}



Psi TaylorCoeffClassic::getMNasPsi(int n){
  return wsp_wielomianu_rozwMasPsi_[n];
}

Psi TaylorCoeffClassic::getPhiNasPsi(int n){

   return wsp_wielomianu_rozwPhiasPsi_[n];
}


void TaylorCoeffClassic::computeDrugiRzad(int jmax, long double eps){
  
  TaylorCoeffExact exact_taylor_coeff(r_min_, r_max_, r_horizon_);
  Psi zerowy_rzad_M = exact_taylor_coeff.computeM0(jmax);
  Psi zerowy_rzad_Phi = exact_taylor_coeff.computePhi0(jmax);
  
  Psi pierwszy_rzad_M_eps1 = getMNasPsi(1);
  Psi pierwszy_rzad_Phi_eps1 = getPhiNasPsi(1);

  alg::IntegralNum dob_M2_eps1(r_min_, r_max_, jmax);
  dob_M2_eps1.setZeroOrder(zerowy_rzad_M + pierwszy_rzad_M_eps1*pow(eps,1), zerowy_rzad_Phi + pierwszy_rzad_Phi_eps1*pow(eps,1));
  dob_M2_eps1.setHorizonRadius(r_horizon_);
  dob_M2_eps1.executeIntegralInfinityZero(1, eps);

  alg::IntegralNum dob_M2_eps2(r_min_, r_max_, jmax);
  dob_M2_eps2.setZeroOrder(zerowy_rzad_M + pierwszy_rzad_M_eps1*pow(2.*eps,1), zerowy_rzad_Phi + pierwszy_rzad_Phi_eps1*pow(2.*eps,1));
  dob_M2_eps2.setHorizonRadius(r_horizon_);
  dob_M2_eps2.executeIntegralInfinityZero(1, 2.*eps);
  
  
  alg::IntegralNum dob_M2_eps3(r_min_, r_max_, jmax);
  dob_M2_eps3.setZeroOrder(zerowy_rzad_M + pierwszy_rzad_M_eps1*pow(3.*eps,1), zerowy_rzad_Phi + pierwszy_rzad_Phi_eps1*pow(3.*eps,1));
  dob_M2_eps3.setHorizonRadius(r_horizon_);
  dob_M2_eps3.executeIntegralInfinityZero(1, 3.*eps);
  
  
  Psi m_second_orderh1 = (zerowy_rzad_M*(2.) + 
                          dob_M2_eps1.getMnPlus1()*(-5.) +
                          dob_M2_eps2.getMnPlus1()*(4.) +
                          dob_M2_eps3.getMnPlus1()*(-1.))*pow(1/eps,2)*(1./2.);
  
  
  if(wsp_wielomianu_rozwMasPsi_.size() == 2){
    wsp_wielomianu_rozwMasPsi_.push_back(m_second_orderh1);
  }
  else if (wsp_wielomianu_rozwMasPsi_.size() > 2) {
    wsp_wielomianu_rozwMasPsi_[2] = m_second_orderh1;
  }
  
  Psi phi_second_orderh1 = (zerowy_rzad_Phi*(2.) + 
                            dob_M2_eps1.getPhinPlus1()*(-5.) +
                            dob_M2_eps2.getPhinPlus1()*(4.) +
                            dob_M2_eps3.getPhinPlus1()*(-1.))*pow(1/eps,2)*(1./2.);
  
  
  if(wsp_wielomianu_rozwPhiasPsi_.size() == 2){
    wsp_wielomianu_rozwPhiasPsi_.push_back(phi_second_orderh1);
  }
  else if (wsp_wielomianu_rozwPhiasPsi_.size() > 2) {
    wsp_wielomianu_rozwPhiasPsi_[2] = phi_second_orderh1;
  }
}

void TaylorCoeffClassic::computeTrzeciRzad(int jmax, long double eps) {
  
  TaylorCoeffExact exact_taylor_coeff(r_min_, r_max_, r_horizon_);
  Psi zerowy_rzad_M = exact_taylor_coeff.computeM0(jmax);
  Psi zerowy_rzad_Phi = exact_taylor_coeff.computePhi0(jmax);
  
  Psi pierwszy_rzad_M_eps1 = getMNasPsi(1);
  Psi pierwszy_rzad_Phi_eps1 = getPhiNasPsi(1);
  
  Psi drugi_rzad_M_eps1 = getMNasPsi(2);
  Psi drugi_rzad_Phi_eps1 = getPhiNasPsi(2);
  
  
  alg::IntegralNum dob_M3_eps1(r_min_, r_max_, jmax);
  dob_M3_eps1.setZeroOrder(zerowy_rzad_M + pierwszy_rzad_M_eps1*pow(eps,1) + 
                           drugi_rzad_M_eps1*pow(eps,2), zerowy_rzad_Phi + 
                           pierwszy_rzad_Phi_eps1*pow(eps,1) + 
                           drugi_rzad_Phi_eps1*pow(eps,2));
  
  dob_M3_eps1.setHorizonRadius(r_horizon_);
  dob_M3_eps1.executeIntegralInfinityZero(1, eps);
  

  alg::IntegralNum dob_M3_eps2(r_min_, r_max_, jmax);
  dob_M3_eps2.setZeroOrder(zerowy_rzad_M + pierwszy_rzad_M_eps1*pow(2.*eps,1) + 
                           drugi_rzad_M_eps1*pow(2.*eps,2), zerowy_rzad_Phi + 
                           pierwszy_rzad_Phi_eps1*pow(2.*eps,1) + 
                           drugi_rzad_Phi_eps1*pow(2.*eps,2));
  
  dob_M3_eps2.setHorizonRadius(r_horizon_);
  dob_M3_eps2.executeIntegralInfinityZero(1, 2.*eps);
  
  
  alg::IntegralNum dob_M3_eps3(r_min_, r_max_, jmax);
  dob_M3_eps3.setZeroOrder(zerowy_rzad_M + pierwszy_rzad_M_eps1*pow(3.*eps,1) + 
                           drugi_rzad_M_eps1*pow(3.*eps,2), zerowy_rzad_Phi + 
                           pierwszy_rzad_Phi_eps1*pow(3.*eps,1) + 
                           drugi_rzad_Phi_eps1*pow(3.*eps,2));
  dob_M3_eps3.setHorizonRadius(r_horizon_);
  dob_M3_eps3.executeIntegralInfinityZero(1, 3.*eps);
 
  
  alg::IntegralNum dob_M3_eps4(r_min_, r_max_, jmax);
  dob_M3_eps4.setZeroOrder(zerowy_rzad_M + pierwszy_rzad_M_eps1*pow(4.*eps,1) + 
                           drugi_rzad_M_eps1*pow(4.*eps,2), zerowy_rzad_Phi + 
                           pierwszy_rzad_Phi_eps1*pow(4.*eps,1) + 
                           drugi_rzad_Phi_eps1*pow(4.*eps,2));
  dob_M3_eps4.setHorizonRadius(r_horizon_);
  dob_M3_eps4.executeIntegralInfinityZero(1, 4.*eps);
  
  Psi m_third_orderh1 = (zerowy_rzad_M*(-5./2.) + 
                          dob_M3_eps1.getMnPlus1()*(9.) +
                          dob_M3_eps2.getMnPlus1()*(-12.) +
                          dob_M3_eps3.getMnPlus1()*(7.) +
                          dob_M3_eps4.getMnPlus1()*(-3./2.) )*pow(1./eps,3)*(1./6.);
                          
  
  if(wsp_wielomianu_rozwMasPsi_.size() == 3){
    wsp_wielomianu_rozwMasPsi_.push_back(m_third_orderh1);
    
  } 
  else if (wsp_wielomianu_rozwMasPsi_.size() > 3) {
    wsp_wielomianu_rozwMasPsi_[3] = m_third_orderh1;
  }
  
  Psi phi_third_orderh1 = (zerowy_rzad_Phi*(-5./2.) + 
                          dob_M3_eps1.getPhinPlus1()*(9.) +
                          dob_M3_eps2.getPhinPlus1()*(-12.) +
                          dob_M3_eps3.getPhinPlus1()*(7.) + 
                          dob_M3_eps4.getPhinPlus1()*(-3./2.) )*pow(1./eps,3)*(1./6.);
  
                          
  if(wsp_wielomianu_rozwPhiasPsi_.size() == 3){
    wsp_wielomianu_rozwPhiasPsi_.push_back(phi_third_orderh1);
  } 
  else if (wsp_wielomianu_rozwPhiasPsi_.size() > 3) {
    wsp_wielomianu_rozwPhiasPsi_[3] = phi_third_orderh1;
  }
}




long double TaylorCoeffClassic::convFactorM1DN(int jmax, long double eps, int rzad_poch) {

  computePierwszyRzad(jmax * 1., eps);
  Psi m_first_order4h = getMNasPsi(1);
  
  computePierwszyRzad(jmax * 2., eps);
  Psi m_first_order2h = getMNasPsi(1);
  
  computePierwszyRzad(jmax * 4., eps);
  Psi m_first_order1h = getMNasPsi(1);

  ConvFactor conv_factor;
  conv_factor.setPsih(m_first_order1h);
  conv_factor.setPsi2h(m_first_order2h);
  conv_factor.setPsi4h(m_first_order4h);
  
  if (rzad_poch == 0){
    return conv_factor.computeConvFactor();}
  else {
    return conv_factor.computeConvFactorD(rzad_poch);}
}


long double TaylorCoeffClassic::convFactorPhi1DN(int jmax, long double eps, int rzad_poch){
  computePierwszyRzad(jmax * 1., eps);
  Psi phi_first_order4h = getPhiNasPsi(1);
  
  computePierwszyRzad(jmax * 2., eps);
  Psi phi_first_order2h = getPhiNasPsi(1);
  
  computePierwszyRzad(jmax * 4., eps);
  Psi phi_first_order1h = getPhiNasPsi(1);

  ConvFactor conv_factor;
  conv_factor.setPsih(phi_first_order1h);
  conv_factor.setPsi2h(phi_first_order2h);
  conv_factor.setPsi4h(phi_first_order4h);
  
  if (rzad_poch == 0){
    return conv_factor.computeConvFactor();}
  else {
    return conv_factor.computeConvFactorD(rzad_poch);}

}


long double TaylorCoeffClassic::convFactorM2(int jmax, long double eps) {
  computePierwszyRzad(jmax * 1., eps);
  computeDrugiRzad(jmax * 1., eps);
  Psi m_second_order4h = getMNasPsi(2);
  
  computePierwszyRzad(jmax * 2., eps);
  computeDrugiRzad(jmax * 2., eps);
  Psi m_second_order2h = getMNasPsi(2);
  
  computePierwszyRzad(jmax * 4., eps);
  computeDrugiRzad(jmax * 4., eps);
  Psi m_second_order1h = getMNasPsi(2);


  ConvFactor conv_factor;
  conv_factor.setPsih(m_second_order1h);
  conv_factor.setPsi2h(m_second_order2h);
  conv_factor.setPsi4h(m_second_order4h);
  
  return conv_factor.computeConvFactor();
}


void TaylorCoeffClassic::plotConvFactorM2ToPs(long double hdol, long double hgora, long double hstep, long double eps){
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
   g1.savetops("PlotConvergenceFactorM2");
   g1.reset_plot();
   g1.set_grid();
   g1.set_style("linespoints").plot_xy(h_x, error_y,"convFactorM2 od h");
   g1.showonscreen();
}

void TaylorCoeffClassic::plotConvFactorM1DNToPs(long double hdol, long double hgora, long double hstep, long double eps, int rzad_poch){

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
  std::string nazwa_pliku_ps = "PlotConvergenceFactorDM1poch"+std::to_string(rzad_poch);
  g1.savetops(nazwa_pliku_ps);
  g1.reset_plot();
  g1.set_grid();
  std::string nazwa_wykresu = "converFactorDM1poch" + std::to_string(rzad_poch);
  g1.set_style("linespoints").plot_xy(h_x, error_y,nazwa_wykresu);
  g1.showonscreen();
}

void TaylorCoeffClassic::plotConvFactorPhi1DNToPs(long double hdol, long double hgora, long double hstep, long double eps, int rzad_poch){
        
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
  std::string nazwa_pliku_ps = "PlotConvergenceFactorDPhi1poch"+std::to_string(rzad_poch);
  g1.savetops(nazwa_pliku_ps);
  g1.reset_plot();
  g1.set_grid();
  std::string nazwa_wykresu = "converFactorDPhi1poch" + std::to_string(rzad_poch);
  g1.set_style("linespoints").plot_xy(h_x, error_y,nazwa_wykresu);
  g1.showonscreen();
}


void TaylorCoeffClassic::plotErrorM2ToPs(long double hdol, long double hgora, long double hstep, long double eps) {
  
  std::cout << "plotErrorM2ToPs" << std::endl;
  
  Gnuplot g1("lines");
  std::vector<long double> error_y, h_x;

  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    
    std::cout <<"h="<<std::to_string(h) <<std::endl;
    
    h_x.push_back(h);
    
    TaylorCoeffExact exact_taylor_coeff(r_min_, r_max_, r_horizon_);
    computePierwszyRzad(i, eps);
    computeDrugiRzad(i, eps);
    error_y.push_back(Error::computeError(getMNasPsi(2), exact_taylor_coeff.computeM2(i, eps)));
    
  }

  g1.reset_all();
  g1.savetops("Err M2 classic dh=" + std::to_string(hstep));
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("linespoints").plot_xy(h_x, error_y,"Error M2 classic dh=" + std::to_string(hstep));
  g1.showonscreen();
}

void TaylorCoeffClassic::plotErrorPhi2ToPs(long double hdol, long double hgora, long double hstep, long double eps) {
  
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
    error_y.push_back(Error::computeError(getPhiNasPsi(2), exact_taylor_coeff.computePhi2(i, eps)));
    
  }

  g1.reset_all();
  g1.savetops("Err Phi2 classic dh=" + std::to_string(hstep));
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("linespoints").plot_xy(h_x, error_y,"Error Phi2 classic dh=" + std::to_string(hstep));
  g1.showonscreen();

}

void TaylorCoeffClassic::plotErrorM1DNToPs(long double hdol, long double hgora, long double hstep, long double eps, int  rzad_poch) {
  
  std::cout << "plotErrorM1DNToPs" << std::endl;
  
  Gnuplot g1("lines");
  std::vector<long double> error_y, h_x;
  
  TaylorCoeffExact exact(r_min_, r_max_, r_horizon_);
  
  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    
    std::cout << "h=" << h << std::endl;
    
    h_x.push_back(h);
    computePierwszyRzad(i, eps);
    if (rzad_poch == 0) {error_y.push_back(Error::computeError( getMNasPsi(1), exact.computeM1DN(i, eps, 0)));}
    if (rzad_poch == 1) {error_y.push_back(Error::computeErrorD( 1,getMNasPsi(1), exact.computeM1DN(i, eps, 1)));}
    if (rzad_poch == 2) {error_y.push_back(Error::computeErrorD( 2,getMNasPsi(1), exact.computeM1DN(i, eps, 2)));}
    if (rzad_poch == 3) {error_y.push_back(Error::computeErrorD( 3,getMNasPsi(1), exact.computeM1DN(i, eps, 3)));}
    if (rzad_poch == 4) {error_y.push_back(Error::computeErrorD( 4,getMNasPsi(1), exact.computeM1DN(i, eps, 4)));}
    
  }

  g1.reset_all();
  g1.savetops("Err M1D"+ std::to_string(rzad_poch)+" classic dh=" + std::to_string(hstep));
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("linespoints").plot_xy(h_x, error_y,"Error M1D"+std::to_string(rzad_poch)+" classic dh=" + std::to_string(hstep));
  g1.showonscreen();
}



void TaylorCoeffClassic::plotRelativeErrorM1ToPs(long double hdol, long double hgora, long double hstep, long double eps) {
  Gnuplot g1("lines");
  std::vector<long double> error_y, h_x;

  for (long double h = hdol; h < hgora; h += hstep) {
    
    std::cout <<"h="+std::to_string(h) <<std::endl;
    
    h_x.push_back(h);
    
    int i = (r_max_ - r_min_) / h;
    computePierwszyRzad(i, eps);
    Psi MnI(i, h);
    MnI = getMNasPsi(1);
    
    computePierwszyRzad(2*i, eps);
    Psi Mn2I(2*i,h/2);
    Mn2I = getMNasPsi(1);
    
    error_y.push_back(Error::computeErrorIvs2I( MnI, Mn2I));
  }

  g1.reset_all();
  g1.savetops("Rel Error M1 classic dh=" + std::to_string(hstep));
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("linespoints").plot_xy(h_x, error_y,"Relative Error M1 classic dh=" + std::to_string(hstep));
  g1.showonscreen();
}

void TaylorCoeffClassic::plotRelativeErrorM2ToPs(long double hdol, long double hgora, long double hstep, long double eps) {
  Gnuplot g1("lines");
  std::vector<long double> error_y, h_x;

  for (long double h = hdol; h < hgora; h += hstep) {
    
    std::cout <<"h="+std::to_string(h) <<std::endl;
    
    h_x.push_back(h);
    
    int i = (r_max_ - r_min_) / h;
    computePierwszyRzad(i, eps);
    computeDrugiRzad(i, eps);
    Psi MnI = getMNasPsi(2);
    
    computePierwszyRzad(2*i, eps);
    computeDrugiRzad(2*i, eps);
    Psi Mn2I = getMNasPsi(2);
    
    error_y.push_back(Error::computeErrorIvs2I( MnI, Mn2I));
    
  }

  g1.reset_all();
  g1.savetops("PlotRelativeErrorM2");
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("linespoints").plot_xy(h_x, error_y,"Relative Error M2 od h");
  g1.showonscreen();
}


void TaylorCoeffClassic::plotErrorM1ExactVsExactDNToPs(long double hdol, long double hgora, long double hstep, long double eps, int rzad_poch) {
  
  Gnuplot g1("lines");

  std::vector<long double> error_y, h_x;
  TaylorCoeffExact exact_taylor_coeff(r_min_, r_max_, r_horizon_);
  
  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    h_x.push_back(h);

    if (rzad_poch == 1) {error_y.push_back(Error::computeErrorD( 1,exact_taylor_coeff.computeM1DN(i, eps, 0), exact_taylor_coeff.computeM1DN(i, eps, 1)));}
    if (rzad_poch == 2) {error_y.push_back(Error::computeErrorD( 2,exact_taylor_coeff.computeM1DN(i, eps, 0), exact_taylor_coeff.computeM1DN(i, eps, 2)));}
    if (rzad_poch == 3) {error_y.push_back(Error::computeErrorD( 3,exact_taylor_coeff.computeM1DN(i, eps, 0), exact_taylor_coeff.computeM1DN(i, eps, 3)));}
    if (rzad_poch == 4) {error_y.push_back(Error::computeErrorD( 4,exact_taylor_coeff.computeM1DN(i, eps, 0), exact_taylor_coeff.computeM1DN(i, eps, 4)));}
    
    std::cout <<"h="+std::to_string(h) <<std::endl;
  }

  g1.reset_all();
  g1.savetops("PlotErrorM1ExactD"+ std::to_string(rzad_poch));
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("linespoints").plot_xy(h_x, error_y,"Error M1ExactD"+std::to_string(rzad_poch)+" od h");
  g1.showonscreen();
}

void TaylorCoeffClassic::plotM1vsM1ExactoPs(long double h, long double r_gorny_zakres, long double eps) {
  
  std::cout << "plotM1vsM1ExactoPs" << std::endl;
  
  int jmax = (r_max_ - r_min_) / h;
  
  TaylorCoeffExact exact(r_min_, r_max_, r_horizon_);
  computePierwszyRzad(jmax, eps);
  
  plotTwoPsiToPs(r_min_, r_max_, r_gorny_zakres, getMNasPsi(1), exact.computeM1DN(jmax, eps, 0), 
                 "M1 classic h="+std::to_string(h)+" zak="+std::to_string(r_gorny_zakres));
}

void TaylorCoeffClassic::plotPhi1vsPhi1ExactoPs(long double h, long double r_gorny_zakres, long double eps) {
  
  std::cout << "plotPhi1vsPhi1ExactoPs" << std::endl;
  
  int jmax = (r_max_ - r_min_) / h;
  
  TaylorCoeffExact exact(r_min_, r_max_, r_horizon_);
  computePierwszyRzad(jmax, eps);
  
  plotTwoPsiToPs(r_min_, r_max_, r_gorny_zakres, getPhiNasPsi(1), exact.computePhi1DN(jmax, eps, 0), 
                 "Phi1 classic h="+std::to_string(h)+" zak="+std::to_string(r_gorny_zakres));
}

void TaylorCoeffClassic::plotM2vsM2ExactoPs(long double h, long double r_gorny_zakres, long double eps) {
  
  std::cout << "plotM2vsM2ExactoPs" << std::endl;
  
  int jmax = (r_max_ - r_min_) / h;
  
  TaylorCoeffExact exact(r_min_, r_max_, r_horizon_);
  computePierwszyRzad(jmax, eps);
  computeDrugiRzad(jmax, eps);
  
  plotTwoPsiToPs(r_min_, r_max_, r_gorny_zakres, getMNasPsi(2), exact.computeM2(jmax, eps), 
                 "M2 classic h="+std::to_string(h)+" zak="+std::to_string(r_gorny_zakres));
}

void TaylorCoeffClassic::plotPhi2vsPhi2ExactoPs(long double h, long double r_gorny_zakres, long double eps) {
  
  std::cout << "plotPhi2vsPhi2ExactoPs" << std::endl;
  
  int jmax = (r_max_ - r_min_) / h;
  
  TaylorCoeffExact exact(r_min_, r_max_, r_horizon_);
  computePierwszyRzad(jmax, eps);
  computeDrugiRzad(jmax, eps);
  
  plotTwoPsiToPs(r_min_, r_max_, r_gorny_zakres, getPhiNasPsi(2), exact.computePhi2(jmax, eps), 
                 "Phi2 classic h="+std::to_string(h)+" zak="+std::to_string(r_gorny_zakres));
}


void TaylorCoeffClassic::plotErrorPhi1DNToPs(long double hdol, long double hgora, long double hstep, long double eps, int rzad_poch) {
  
  std::cout << "plotErrorPhi1DNToPs" << std::endl;
  
  Gnuplot g1("lines");
  std::vector<long double> error_y, h_x;
  TaylorCoeffExact exact_taylor_coeff(r_min_, r_max_, r_horizon_);

  for (long double h = hdol; h < hgora; h += hstep) {
    int i = (r_max_ - r_min_) / h;
    
    std::cout << "h=" << h << std::endl;
    
    h_x.push_back(h);
    computePierwszyRzad(i, eps);
    if (rzad_poch == 0) {error_y.push_back(Error::computeError(    getPhiNasPsi(1), exact_taylor_coeff.computePhi1DN(i, eps, 0)));}
    if (rzad_poch == 1) {error_y.push_back(Error::computeErrorD(1 ,getPhiNasPsi(1), exact_taylor_coeff.computePhi1DN(i, eps, 1)));}
    if (rzad_poch == 2) {error_y.push_back(Error::computeErrorD(2 ,getPhiNasPsi(1), exact_taylor_coeff.computePhi1DN(i, eps, 2)));}
    if (rzad_poch == 3) {error_y.push_back(Error::computeErrorD(3 ,getPhiNasPsi(1), exact_taylor_coeff.computePhi1DN(i, eps, 3)));}
    if (rzad_poch == 4) {error_y.push_back(Error::computeErrorD(4 ,getPhiNasPsi(1), exact_taylor_coeff.computePhi1DN(i, eps, 4)));}
  }

  g1.reset_all();
  g1.savetops("Err Phi1D"+ std::to_string(rzad_poch) + " classic dh=" + std::to_string(hstep));
  g1.reset_plot();
  g1.set_grid();
  g1.set_style("linespoints").plot_xy(h_x, error_y,"Error Phi1D"+std::to_string(rzad_poch)+ " classic dh=" + std::to_string(hstep));
  g1.showonscreen();
}


void TaylorCoeffClassic::plotTwoPsiToPs(long double rmin, long double rmax, long double r_gorne, Psi psi_przybliz, Psi psi_dokl, std::string nazwa_pliku_ps) {
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
  g1.set_grid();
  g1.set_style("points").plot_xy(x, y_przybliz,"Przyblizony " + nazwa_pliku_ps);
  g1.set_style("lines").plot_xy(x, y_dokl,"Dokladny " + nazwa_pliku_ps);
  g1.showonscreen();
}

}