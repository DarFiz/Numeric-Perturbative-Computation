#pragma once

long double M1(long double eps,long double DM0, long double r);
long double M1D1(long double eps,long double DM0, long double r);
long double M1D2(long double eps,long double DM0, long double r);
long double M1D3(long double eps,long double DM0, long double r);
long double M1D4(long double eps,long double DM0, long double r);
long double Phi1(long double eps, long double DM0, long double r);
long double Phi1D1(long double eps, long double DM0, long double r);
long double Phi1D2(long double eps, long double DM0, long double r);
long double Phi1D3(long double eps, long double DM0, long double r);
long double Phi1D4(long double eps, long double DM0, long double r);

long double M2(long double DM0, long double r);
long double Phi2(long double DM0, long double r);
long double Phi3(long double DM0, long double r);

long double M3(long double DM0, long double r);

long double Ttt(long double DM0, long double DM1, long double DM2, long double DM3, long double DM4,
           long double Dphi0, long double Dphi1, long double Dphi2, long double Dphi3, long double Dphi4,
           long double r);
long double TrrMTttPrzez(long double DM0, long double DM1, long double DM2, long double DM3, long double DM4, 
           long double Dphi0, long double Dphi1, long double Dphi2, long double Dphi3, long double Dphi4,
           long double r);