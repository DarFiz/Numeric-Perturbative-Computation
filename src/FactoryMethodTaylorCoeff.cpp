#include <memory>
#include <string>
#include "ITaylorCoeff.h"
#include "TaylorCoeffClassic.h"
#include "TaylorCoeffFit.h"
#include "FactoryMethodTaylorCoeff.h"

namespace alg {

std::unique_ptr<ITaylorCoeff> FactoryMethodTaylorCoeff::factoryMethod(const std::string& typ_metody_obliczen,
                                                                      long double rmin, 
                                                                      long double rmax, 
                                                                      long double horizon) {
   
  if(typ_metody_obliczen == "classic"){
     return std::unique_ptr<ITaylorCoeff>(new TaylorCoeffClassic(rmin, rmax, horizon));
  }
  else if(typ_metody_obliczen == "fit"){
     return std::unique_ptr<ITaylorCoeff>(new TaylorCoeffFit(rmin, rmax, horizon));
  }
  else {
     return nullptr;
  }
   
}

}