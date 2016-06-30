#pragma once
#include <string>
#include <memory>
#include "ITaylorCoeff.h"
#include "IFactoryMethodTaylorCoeff.h"

namespace alg {

class FactoryMethodTaylorCoeff: public IFactoryMethodTaylorCoeff {
public:
std::unique_ptr<ITaylorCoeff> factoryMethod(const std::string& typ_metody_obliczen,      
                                            long double rmin, long double rmax, 
                                            long double horizon);
};

}