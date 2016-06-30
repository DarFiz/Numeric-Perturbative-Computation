#include  <cmath>
#include "TabFunc.h"

long double M1(long double eps, long double DM0, long double r){
  return (196.*pow(DM0, 3)/pow(r,6)-108.*pow(DM0,2)/pow(r,5));
}

long double M1D1(long double eps, long double DM0, long double r){
  return (540.*pow(DM0, 2)/pow(r,6)-1176.*pow(DM0,3)/pow(r,7));
}

long double M1D2(long double eps, long double DM0, long double r){
  return (-3240.*pow(DM0, 2)/pow(r,7)+8232.*pow(DM0,3)/pow(r,8));
}

long double M1D3(long double eps, long double DM0, long double r){
  return (22680.*pow(DM0, 2)/pow(r,8)-65856.*pow(DM0,3)/pow(r,9));
}


long double M1D4(long double eps, long double DM0, long double r){
  return (-181440.*pow(DM0, 2)/pow(r,9)+592704.*pow(DM0,3)/pow(r,10));
}

long double Phi1(long double eps, long double DM0, long double r){
return (-108.*pow(DM0, 2)/pow(r, 6));
}

long double Phi1D1(long double eps, long double DM0, long double r){
return (648.*pow(DM0, 2)/pow(r, 7));
}

long double Phi1D2(long double eps, long double DM0, long double r){
return (-4536.*pow(DM0, 2)/pow(r, 8));
}

long double Phi1D3(long double eps, long double DM0, long double r){
return (36288.*pow(DM0, 2)/pow(r, 9));
}

long double Phi1D4(long double eps, long double DM0, long double r){
return (-326592.*pow(DM0, 2)/pow(r, 10));
}

long double M2(long double DM0, long double r){
  return (14808.*pow(DM0, 5)/pow(r,12)-98496./11.*pow(DM0,4)/pow(r,11));
}

long double Phi2(long double DM0, long double r){
return (165888./11.*pow(DM0, 3)/pow(r, 11)-29808.*pow(DM0, 4)/pow(r, 12));
}


long double M3(long double DM0, long double r){
  return (672841728.*pow(DM0, 4)/pow(r,15)-52295621760./11.*pow(DM0,5)/pow(r,16)+2063604936960./187.*pow(DM0, 6)/pow(r,17)-8433504096.*pow(DM0,7)/pow(r,18));
}

long double Phi3(long double DM0, long double r){
  return (626030208.*pow(DM0,4)/pow(r,16)-603943547904./187.*pow(DM0, 5)/pow(r,17)+4083264000.*pow(DM0,6)/pow(r,18));
}


long double Ttt(long double DM0, long double DM1, long double DM2, long double DM3, long double DM4, 
           long double Dphi0, long double Dphi1, long double Dphi2, long double Dphi3, long double Dphi4,
           long double r){
  return 
(48/r*DM0*Dphi1*Dphi3-48*DM1*Dphi1*Dphi3-12*pow(Dphi1,4)-2316/pow(r,4)*pow(DM0,2
)*pow(Dphi1,2)*DM1-2100/pow(r,4)*pow(DM0,3)*Dphi1*Dphi2-132*DM3*Dphi2*DM1-96*DM3
*Dphi3*DM0+4344/pow(r,5)*pow(DM0,2)*DM1*Dphi1-72/pow(r,2)*pow(DM0,2)*Dphi1*DM4+
48/r*pow(Dphi1,2)*pow(DM0,2)*DM4+336/pow(r,3)*DM1*DM0*DM3+96/r*pow(Dphi1,2)*pow(
DM0,3)*Dphi4-144/pow(r,2)*pow(DM0,3)*Dphi1*Dphi4-1104/pow(r,4)*DM1*DM0*DM2-3048/
pow(r,4)*DM1*pow(DM0,2)*Dphi2-2136/pow(r,4)*pow(DM1,2)*Dphi1*DM0-1716/pow(r,4)*
pow(DM0,2)*Dphi1*DM2+240/r*pow(Dphi1,2)*DM0*pow(DM2,2)+192/pow(r,2)*pow(DM0,3)*
pow(Dphi1,3)*Dphi2-108/pow(r,3)*DM1*pow(DM0,2)*pow(Dphi1,3)+60/r*pow(DM1,2)*DM2*
Dphi2+1560/pow(r,3)*pow(DM1,2)*DM0*Dphi2+720/r*pow(DM1,2)*pow(Dphi2,2)*DM0-1296/
pow(r,2)*DM1*pow(DM0,2)*pow(Dphi2,2)+1092/pow(r,3)*pow(DM0,2)*DM2*Dphi2-816/pow(
r,2)*pow(DM0,3)*Dphi1*pow(Dphi2,2)-360/pow(r,2)*DM0*Dphi1*pow(DM2,2)+960/pow(r,3
)*DM1*pow(DM0,2)*Dphi3+24/r*DM2*DM4*DM0-48/pow(r,2)*DM1*DM4*DM0+96/r*DM3*Dphi3*
pow(DM0,2)+96/r*Dphi2*pow(DM0,3)*Dphi4-96/pow(r,2)*DM1*Dphi4*pow(DM0,2)+48/r*DM2
*Dphi4*pow(DM0,2)+48/r*Dphi2*pow(DM0,2)*DM4+12/r*DM3*DM2*DM1-156/pow(r,2)*DM3*
DM2*DM0-288/pow(r,2)*Dphi3*pow(DM0,2)*DM2-480/pow(r,2)*Dphi3*pow(DM0,3)*Dphi2-
264/pow(r,2)*DM3*Dphi2*pow(DM0,2)+432/r*DM2*pow(Dphi2,2)*pow(DM0,2)+168/r*pow(
DM2,2)*Dphi2*DM0+144/r*pow(Dphi1,3)*pow(DM0,2)*DM3+48/r*pow(Dphi1,4)*pow(DM0,2)*
DM2+192/r*pow(Dphi1,3)*pow(DM0,3)*Dphi3-96/r*pow(Dphi1,4)*pow(DM0,3)*Dphi2-336/
pow(r,2)*pow(DM0,2)*pow(Dphi1,2)*DM3-216/pow(r,2)*pow(DM0,2)*pow(Dphi1,3)*DM2+
480/pow(r,2)*pow(DM0,2)*pow(Dphi1,4)*DM1-480/pow(r,2)*pow(DM0,3)*pow(Dphi1,2)*
Dphi3-144/r*DM1*pow(Dphi1,5)*pow(DM0,2)+468/pow(r,3)*pow(DM0,2)*Dphi1*DM3+984/
pow(r,3)*pow(DM0,2)*pow(Dphi1,2)*DM2+768/pow(r,3)*pow(DM0,3)*Dphi1*Dphi3+672/pow
(r,3)*pow(DM0,3)*pow(Dphi1,2)*Dphi2+108/pow(r,2)*pow(DM1,2)*DM0*pow(Dphi1,3)+36/
r*pow(DM1,2)*Dphi1*DM3-36/pow(r,2)*pow(DM1,2)*Dphi1*DM2+24/r*pow(DM1,2)*pow(
Dphi1,2)*DM2+804/pow(r,3)*pow(DM1,2)*pow(Dphi1,2)*DM0+180/r*pow(DM1,3)*Dphi1*
Dphi2-72/r*pow(DM1,2)*pow(Dphi1,4)*DM0-384/pow(r,2)*pow(DM1,2)*DM0*Dphi3+480/r*
pow(Dphi1,2)*pow(DM0,3)*pow(Dphi2,2)+1464/r*DM0*DM1*DM2*Dphi2*Dphi1+24/r*Dphi1*
Dphi2-24/pow(r,2)*Dphi1*DM2+72/pow(r,3)*Dphi1*DM1-72/pow(r,4)*Dphi1*DM0-24*Dphi3
*DM2+48/r*Dphi3*DM1-48/pow(r,2)*Dphi3*DM0-48*r*Dphi3*Dphi2-48*r*Dphi3*pow(Dphi1,
2)-96*r*pow(Dphi2,2)*Dphi1-96*r*Dphi2*pow(Dphi1,3)+36*r*pow(Dphi1,3)*DM3+36*pow(
Dphi1,3)*DM2-48/pow(r,2)*pow(Dphi1,3)*DM0-24*r*pow(Dphi1,6)*DM0+12*r*pow(Dphi1,2
)*DM4+24/r*pow(Dphi1,2)*DM2-132/pow(r,2)*pow(Dphi1,2)*DM1+204/pow(r,3)*pow(Dphi1
,2)*DM0+48*pow(Dphi1,6)*pow(DM0,2)-60*pow(r,2)*pow(Dphi2,2)*pow(Dphi1,2)-12*pow(
r,2)*Dphi2*Dphi4+12*pow(r,2)*Dphi2*pow(Dphi1,4)-24*pow(r,2)*pow(Dphi1,3)*Dphi3-
12*pow(r,2)*pow(Dphi1,2)*Dphi4-28/pow(r,4)*pow(Dphi1,3)*pow(DM0,3)-576/pow(r,4)*
pow(DM0,3)*Dphi3-312/pow(r,4)*pow(DM0,2)*DM3+48/pow(r,3)*pow(DM0,2)*DM4+96/pow(r
,3)*pow(DM0,3)*Dphi4+96/r*pow(Dphi3,2)*pow(DM0,3)+24/r*pow(DM3,2)*DM0+672/pow(r,
3)*pow(DM0,3)*pow(Dphi2,2)+240/pow(r,3)*DM0*pow(DM2,2)-24/pow(r,2)*pow(DM1,2)*
DM3+48/pow(r,3)*pow(DM1,2)*DM2+160/r*pow(Dphi2,3)*pow(DM0,3)-32/r*pow(Dphi1,6)*
pow(DM0,3)-120/pow(r,2)*pow(DM1,3)*Dphi2-312/pow(r,3)*pow(DM0,3)*pow(Dphi1,4)+
144/pow(r,2)*pow(DM0,3)*pow(Dphi1,5)+60/pow(r,2)*pow(DM1,3)*pow(Dphi1,2)-36/r*
pow(DM1,3)*pow(Dphi1,3)-2400/pow(r,6)*pow(DM0,2)*DM1-2280/pow(r,6)*pow(DM0,3)*
Dphi1+1056/pow(r,5)*pow(DM0,2)*DM2+1320/pow(r,5)*DM0*pow(DM1,2)+1608/pow(r,5)*
pow(DM0,3)*Dphi2+1548/pow(r,5)*pow(DM0,3)*pow(Dphi1,2)+72/pow(r,3)*pow(DM1,3)*
Dphi1-24*Dphi2*pow(Dphi1,2)+72*r*pow(Dphi3,2)*DM0-144*pow(Dphi3,2)*pow(DM0,2)+24
*r*DM3*Dphi3-360*pow(Dphi2,2)*pow(DM1,2)+180/r*pow(DM2,2)*Dphi1-84*pow(DM2,2)*
Dphi2-120*pow(DM2,2)*pow(Dphi1,2)-12*DM2*DM4+72/r*DM2*DM3+504/pow(r,3)*DM2*DM1-
504/pow(r,4)*DM2*DM0+12*r*DM2*Dphi4+108*r*DM2*pow(Dphi2,2)+12*r*DM2*pow(Dphi1,4)
+960/pow(r,3)*pow(DM1,2)*Dphi1-600/pow(r,2)*pow(DM1,2)*Dphi2-492/pow(r,2)*pow(
DM1,2)*pow(Dphi1,2)+192/r*pow(DM1,2)*Dphi3+24/r*DM1*DM4-144/pow(r,2)*DM1*DM3+
1080/pow(r,5)*DM1*DM0-24*DM1*Dphi4+36*DM1*pow(Dphi2,2)+84*DM1*pow(Dphi1,4)+1248/
pow(r,5)*pow(DM0,2)*Dphi1-888/pow(r,4)*pow(DM0,2)*Dphi2-1164/pow(r,4)*pow(DM0,2)
*pow(Dphi1,2)+384/pow(r,3)*pow(DM0,2)*Dphi3+96/pow(r,3)*pow(DM0,2)*pow(Dphi1,3)-
96/pow(r,2)*pow(DM0,2)*Dphi4-360/pow(r,2)*pow(DM0,2)*pow(Dphi2,2)+228/pow(r,2)*
pow(DM0,2)*pow(Dphi1,4)-24/pow(r,2)*DM0*DM4+144/pow(r,3)*DM0*DM3+24/r*DM0*Dphi4+
36/r*DM0*pow(Dphi2,2)-12/r*DM0*pow(Dphi1,4)+36*pow(Dphi1,4)*pow(DM1,2)-36*r*pow(
Dphi1,5)*DM1-144/r*pow(Dphi1,5)*pow(DM0,2)+36*pow(Dphi1,5)*DM0+120*r*pow(Dphi2,3
)*DM0+12*r*Dphi2*DM4+24/r*Dphi2*DM2-72/pow(r,2)*Dphi2*DM1+72/pow(r,3)*Dphi2*DM0-
240*pow(Dphi2,3)*pow(DM0,2)+216/r*Dphi1*DM0*DM3*DM2+336/r*pow(DM0,2)*Dphi2*Dphi1
*DM3+264/r*DM3*Dphi2*DM0*DM1+192/r*Dphi3*DM0*DM2*DM1+864/r*Dphi3*pow(DM0,2)*
Dphi2*DM1+144/r*DM1*Dphi1*pow(DM0,2)*Dphi4+72/r*DM1*Dphi1*DM0*DM4-1344/pow(r,2)*
DM1*pow(DM0,2)*Dphi1*Dphi3+336/r*DM1*pow(Dphi1,2)*DM0*DM3+864/r*DM1*pow(Dphi1,2)
*pow(DM0,2)*Dphi3+384/r*Dphi1*pow(DM0,2)*Dphi3*DM2+576/r*Dphi1*pow(DM0,3)*Dphi3*
Dphi2+1656/pow(r,3)*DM1*Dphi1*DM2*DM0-1884/pow(r,2)*pow(DM1,2)*DM0*Dphi2*Dphi1+
576/r*pow(DM1,2)*Dphi1*Dphi3*DM0+1152/r*pow(DM1,2)*pow(Dphi1,2)*Dphi2*DM0+3612/
pow(r,3)*DM1*pow(DM0,2)*Dphi2*Dphi1-1632/pow(r,2)*pow(DM0,2)*pow(Dphi1,2)*Dphi2*
DM1+1584/r*DM1*Dphi1*pow(Dphi2,2)*pow(DM0,2)-1512/pow(r,2)*pow(DM0,2)*Dphi1*DM2*
Dphi2+960/r*pow(Dphi1,2)*pow(DM0,2)*DM2*Dphi2-912/pow(r,2)*DM1*DM0*pow(Dphi1,2)*
DM2+360/r*DM1*pow(Dphi1,3)*DM0*DM2+576/r*pow(Dphi1,3)*pow(DM0,2)*Dphi2*DM1-504/
pow(r,2)*DM1*DM0*Dphi1*DM3-1056/pow(r,2)*DM1*DM0*DM2*Dphi2-732*DM2*Dphi2*Dphi1*
DM1+780/r*DM2*Dphi2*Dphi1*DM0-384*DM2*Dphi3*Dphi1*DM0-960*DM2*Dphi2*pow(Dphi1,2)
*DM0-1728/pow(r,2)*DM1*Dphi2*Dphi1*DM0+768/r*DM1*Dphi3*Dphi1*DM0+480/r*DM1*Dphi2
*pow(Dphi1,2)*DM0-144*Dphi1*DM1*Dphi4*DM0-1584*Dphi1*DM1*pow(Dphi2,2)*DM0-864*
pow(Dphi1,2)*DM1*Dphi3*DM0-576*pow(Dphi1,3)*DM1*Dphi2*DM0+432*r*Dphi2*Dphi3*
Dphi1*DM0-336*Dphi2*DM0*Dphi1*DM3-864*Dphi2*DM0*Dphi3*DM1-864*Dphi2*pow(DM0,2)*
Dphi3*Dphi1-72*r*Dphi2*pow(Dphi1,4)*DM0+360*r*pow(Dphi2,2)*pow(Dphi1,2)*DM0+288/
r*Dphi2*pow(DM0,2)*Dphi3-144*Dphi2*pow(DM0,2)*Dphi4+144*Dphi2*pow(DM0,2)*pow(
Dphi1,4)-48*Dphi2*DM0*DM4+132/r*Dphi2*DM0*DM3-720*pow(Dphi2,2)*pow(DM0,2)*pow(
Dphi1,2)+72*r*pow(Dphi1,2)*Dphi4*DM0+144*r*pow(Dphi1,3)*Dphi3*DM0-144*pow(Dphi1,
3)*DM0*DM3-144*pow(Dphi1,2)*pow(DM0,2)*Dphi4-48*pow(Dphi1,2)*DM0*DM4-288*pow(
Dphi1,3)*pow(DM0,2)*Dphi3-72*pow(r,2)*Dphi2*Dphi3*Dphi1-108*DM2*Dphi1*DM3-792/
pow(r,2)*DM2*Dphi1*DM1+888/pow(r,3)*DM2*Dphi1*DM0+468/r*DM2*Dphi2*DM1-564/pow(r,
2)*DM2*Dphi2*DM0+432/r*DM2*pow(Dphi1,2)*DM1-528/pow(r,2)*DM2*pow(Dphi1,2)*DM0-96
*DM2*Dphi3*DM1+192/r*DM2*Dphi3*DM0-180*DM2*pow(Dphi1,3)*DM1+36/r*DM2*pow(Dphi1,3
)*DM0-48*DM2*Dphi4*DM0+96*r*DM2*Dphi3*Dphi1-432*DM2*pow(Dphi2,2)*DM0+240*r*DM2*
Dphi2*pow(Dphi1,2)-48*DM2*pow(Dphi1,4)*DM0+216/r*DM1*Dphi1*DM3-2208/pow(r,4)*DM1
*Dphi1*DM0+1488/pow(r,3)*DM1*Dphi2*DM0+1512/pow(r,3)*DM1*pow(Dphi1,2)*DM0-576/
pow(r,2)*DM1*Dphi3*DM0+96/r*DM1*Dphi4*DM0+576/r*DM1*pow(Dphi2,2)*DM0+168*DM1*
Dphi2*pow(Dphi1,2)-408/r*DM1*pow(Dphi1,4)*DM0+672/r*pow(DM1,2)*Dphi2*Dphi1-216/
pow(r,2)*DM0*Dphi1*DM3-24/r*DM0*Dphi2*pow(Dphi1,2)+1344/pow(r,3)*pow(DM0,2)*
Dphi2*Dphi1-480/pow(r,2)*pow(DM0,2)*Dphi3*Dphi1-192/pow(r,2)*pow(DM0,2)*Dphi2*
pow(Dphi1,2)-168*pow(Dphi1,2)*DM1*DM3-288*Dphi1*pow(DM1,2)*Dphi3+216*r*pow(Dphi1
,2)*DM1*Dphi3+144*r*pow(Dphi1,3)*DM1*Dphi2+144*pow(Dphi1,5)*DM1*DM0-36*Dphi1*DM1
*DM4+36*r*Dphi1*DM1*Dphi4+396*r*Dphi1*DM1*pow(Dphi2,2)-576*pow(Dphi1,2)*pow(DM1,
2)*Dphi2+168/r*pow(Dphi1,2)*DM0*DM3+144/r*Dphi1*pow(DM0,2)*Dphi4+72*pow(Dphi1,2)
*DM0*Dphi3+432/r*Dphi1*pow(DM0,2)*pow(Dphi2,2)+432*pow(Dphi1,3)*DM0*Dphi2+36/r*
Dphi1*DM0*DM4-36*Dphi1*DM0*Dphi4+180*Dphi1*DM0*pow(Dphi2,2)+288/r*pow(Dphi1,2)*
pow(DM0,2)*Dphi3-576/r*pow(Dphi1,3)*pow(DM0,2)*Dphi2+84*r*Dphi2*Dphi1*DM3-12*
Dphi2*Dphi1*DM2+96/r*Dphi2*Dphi1*DM1-240/pow(r,2)*Dphi2*Dphi1*DM0+216*r*Dphi2*
Dphi3*DM1+72*Dphi2*Dphi3*DM0+72*r*Dphi2*Dphi4*DM0-120/pow(r,2)*pow(DM2,2)-540/
pow(r,4)*pow(DM1,2)-540/pow(r,6)*pow(DM0,2)-20*pow(r,2)*pow(Dphi2,3)+4*pow(r,2)*
pow(Dphi1,6)+8/r*pow(Dphi1,3)-12*pow(DM3,2)-12*pow(r,2)*pow(Dphi3,2)-12*pow(
Dphi2,2)-12/pow(r,2)*pow(Dphi1,2)+1176/pow(r,7)*pow(DM0,3)-4/r*pow(DM2,3)-80/pow
(r,4)*pow(DM1,3));
}

long double TrrMTttPrzez(long double DM0, long double DM1, long double DM2, long double DM3, long double DM4, 
           long double Dphi0, long double Dphi1, long double Dphi2, long double Dphi3, long double Dphi4,
           long double r){
  return 
(24*DM1*pow(Dphi1,5)-24*r*pow(Dphi1,4)*Dphi2-24/r*pow(DM1,2)*pow(Dphi1,4)+24/r*
DM0*pow(Dphi1,5)+12*DM2*pow(Dphi1,4)+12*r*pow(Dphi1,3)*Dphi3+48*r*pow(Dphi1,2)*
pow(Dphi2,2)-48/r*DM1*pow(Dphi1,4)+96/r*pow(DM0,2)*pow(Dphi2,3)-96*DM0*pow(Dphi2
,3)-24*DM3*pow(Dphi1,3)+72*pow(Dphi1,3)*Dphi2+12*r*pow(Dphi1,2)*Dphi4+12/pow(r,2
)*pow(DM1,2)*pow(Dphi1,3)+300/r*pow(DM1,2)*pow(Dphi2,2)-72/pow(r,3)*pow(DM0,2)*
pow(Dphi1,4)+48/r*pow(DM0,2)*pow(Dphi3,2)-48*DM0*pow(Dphi3,2)+96/r*pow(DM2,2)*
pow(Dphi1,2)-48/r*DM2*pow(Dphi1,3)-120*DM2*pow(Dphi2,2)+48*pow(Dphi1,2)*Dphi3-12
*pow(Dphi1,2)*DM4+84*Dphi1*pow(Dphi2,2)+12*r*Dphi2*Dphi4+24/pow(r,2)*DM1*pow(
Dphi1,3)-24*DM3*Dphi3+96/r*pow(DM2,2)*Dphi2-12*DM2*Dphi4+24/r*pow(Dphi1,2)*Dphi2
+48*Dphi2*Dphi3-12*Dphi2*DM4+300/pow(r,3)*pow(DM1,2)*pow(Dphi1,2)-168/pow(r,2)*
pow(DM1,2)*Dphi3+24/r*DM1*Dphi4+108/pow(r,4)*pow(DM0,2)*pow(Dphi1,3)+348/pow(r,3
)*pow(DM0,2)*pow(Dphi2,2)-72/pow(r,3)*DM0*pow(Dphi1,3)-48/pow(r,2)*DM0*pow(Dphi2
,2)-156/pow(r,2)*pow(DM2,2)*Dphi1-24/pow(r,2)*DM2*pow(Dphi1,2)+24/r*DM2*Dphi3+12
/r*DM2*DM4+624/pow(r,3)*pow(DM1,2)*Dphi2+72/pow(r,3)*DM1*pow(Dphi1,2)-48/pow(r,2
)*DM1*Dphi3-24/pow(r,2)*DM1*DM4+48/pow(r,3)*pow(DM0,2)*Dphi4-24/pow(r,2)*DM0*
Dphi4-72/pow(r,2)*DM3*DM2-24/pow(r,2)*DM2*Dphi2-24/pow(r,2)*Dphi1*Dphi2-816/pow(
r,4)*pow(DM1,2)*Dphi1+144/pow(r,3)*DM1*DM3+72/pow(r,3)*DM1*Dphi2+492/pow(r,5)*
pow(DM0,2)*pow(Dphi1,2)-264/pow(r,4)*pow(DM0,2)*Dphi3-120/pow(r,4)*DM0*pow(Dphi1
,2)+48/pow(r,3)*DM0*Dphi3+24/pow(r,3)*DM0*DM4+24/pow(r,3)*DM2*Dphi1-528/pow(r,4)
*DM1*DM2-72/pow(r,4)*DM1*Dphi1+768/pow(r,5)*pow(DM0,2)*Dphi2-144/pow(r,4)*DM0*
DM3-72/pow(r,4)*DM0*Dphi2-960/pow(r,6)*pow(DM0,2)*Dphi1+528/pow(r,5)*DM0*DM2+72/
pow(r,5)*DM0*Dphi1-1200/pow(r,6)*DM1*DM0+24*r*pow(Dphi2,3)+12/r*pow(Dphi1,4)+12*
r*pow(Dphi3,2)+12/r*pow(DM3,2)+12/pow(r,2)*pow(Dphi1,3)+12/r*pow(Dphi2,2)+120/
pow(r,3)*pow(DM2,2)+12/pow(r,3)*pow(Dphi1,2)+600/pow(r,5)*pow(DM1,2)+648/pow(r,7
)*pow(DM0,2)-48/r*DM1*DM0*pow(Dphi1,5)-96/r*pow(DM0,2)*pow(Dphi1,4)*Dphi2+96*DM0
*pow(Dphi1,4)*Dphi2-36*DM1*pow(Dphi1,3)*Dphi2+48/r*pow(DM0,2)*pow(Dphi1,3)*Dphi3
+192/r*pow(DM0,2)*pow(Dphi1,2)*pow(Dphi2,2)-24/r*DM0*DM2*pow(Dphi1,4)-48*DM0*pow
(Dphi1,3)*Dphi3-192*DM0*pow(Dphi1,2)*pow(Dphi2,2)+348/r*pow(DM1,2)*pow(Dphi1,2)*
Dphi2+144/pow(r,2)*DM1*DM0*pow(Dphi1,4)+84/r*DM1*DM2*pow(Dphi1,3)-168*DM1*pow(
Dphi1,2)*Dphi3-348*DM1*Dphi1*pow(Dphi2,2)+216/pow(r,2)*pow(DM0,2)*pow(Dphi1,3)*
Dphi2+48/r*pow(DM0,2)*pow(Dphi1,2)*Dphi4+48/r*DM0*DM3*pow(Dphi1,3)-252/r*DM0*pow
(Dphi1,3)*Dphi2-48*DM0*pow(Dphi1,2)*Dphi4-204*DM2*pow(Dphi1,2)*Dphi2+60*r*Dphi1*
Dphi2*Dphi3+252/r*pow(DM1,2)*Dphi1*Dphi3+120/r*DM1*DM3*pow(Dphi1,2)-144/r*DM1*
pow(Dphi1,2)*Dphi2-36*DM1*Dphi1*Dphi4-204*DM1*Dphi2*Dphi3-144/pow(r,2)*pow(DM0,2
)*pow(Dphi1,2)*Dphi3-360/pow(r,2)*pow(DM0,2)*Dphi1*pow(Dphi2,2)+48/r*pow(DM0,2)*
Dphi2*Dphi4+12/pow(r,2)*DM0*DM2*pow(Dphi1,3)+240/r*DM0*DM2*pow(Dphi2,2)-24/r*DM0
*pow(Dphi1,2)*Dphi3+24/r*DM0*pow(Dphi1,2)*DM4+12/r*DM0*Dphi1*pow(Dphi2,2)-48*DM0
*Dphi2*Dphi4-72*DM3*Dphi1*Dphi2-84*DM2*Dphi1*Dphi3-648/pow(r,2)*pow(DM1,2)*Dphi1
*Dphi2-72/pow(r,3)*DM1*DM0*pow(Dphi1,3)-600/pow(r,2)*DM1*DM0*pow(Dphi2,2)+120/r*
DM1*DM3*Dphi2-288/pow(r,2)*DM1*DM2*pow(Dphi1,2)+84/r*DM1*DM2*Dphi3+24/r*DM1*
Dphi1*Dphi3+36/r*DM1*Dphi1*DM4+156/pow(r,3)*pow(DM0,2)*pow(Dphi1,2)*Dphi2-72/pow
(r,2)*pow(DM0,2)*Dphi1*Dphi4-216/pow(r,2)*pow(DM0,2)*Dphi2*Dphi3-120/pow(r,2)*
DM0*DM3*pow(Dphi1,2)+48/r*DM0*DM3*Dphi3+24/r*DM0*DM2*Dphi4+48/pow(r,2)*DM0*pow(
Dphi1,2)*Dphi2+36/r*DM0*Dphi1*Dphi4+12/r*DM0*Dphi2*Dphi3+24/r*DM0*Dphi2*DM4+96/r
*DM3*DM2*Dphi1-48/pow(r,2)*DM1*DM0*Dphi4-192/pow(r,2)*DM1*DM3*Dphi1-504/pow(r,2)
*DM1*DM2*Dphi2-72/pow(r,2)*DM1*Dphi1*Dphi2+300/pow(r,3)*pow(DM0,2)*Dphi1*Dphi3-
120/pow(r,2)*DM0*DM3*Dphi2+336/pow(r,3)*DM0*DM2*pow(Dphi1,2)-132/pow(r,2)*DM0*
DM2*Dphi3-24/pow(r,2)*DM0*Dphi1*Dphi3-36/pow(r,2)*DM0*Dphi1*DM4-744/pow(r,4)*DM1
*DM0*pow(Dphi1,2)+432/pow(r,3)*DM1*DM0*Dphi3+696/pow(r,3)*DM1*DM2*Dphi1-888/pow(
r,4)*pow(DM0,2)*Dphi1*Dphi2+192/pow(r,3)*DM0*DM3*Dphi1+552/pow(r,3)*DM0*DM2*
Dphi2+168/pow(r,3)*DM0*Dphi1*Dphi2-1392/pow(r,4)*DM1*DM0*Dphi2-744/pow(r,4)*DM0*
DM2*Dphi1+1776/pow(r,5)*DM1*DM0*Dphi1+72/r*DM1*DM0*pow(Dphi1,3)*Dphi2+336/r*DM1*
DM0*pow(Dphi1,2)*Dphi3+696/r*DM1*DM0*Dphi1*pow(Dphi2,2)+240/r*pow(DM0,2)*Dphi1*
Dphi2*Dphi3+408/r*DM0*DM2*pow(Dphi1,2)*Dphi2-240*DM0*Dphi1*Dphi2*Dphi3-408/pow(r
,2)*DM1*DM0*pow(Dphi1,2)*Dphi2+72/r*DM1*DM0*Dphi1*Dphi4+408/r*DM1*DM0*Dphi2*
Dphi3+684/r*DM1*DM2*Dphi1*Dphi2+144/r*DM0*DM3*Dphi1*Dphi2+168/r*DM0*DM2*Dphi1*
Dphi3-552/pow(r,2)*DM1*DM0*Dphi1*Dphi3-684/pow(r,2)*DM0*DM2*Dphi1*Dphi2+1440/pow
(r,3)*DM1*DM0*Dphi1*Dphi2-12*pow(Dphi1,5));
}