#pragma once
#include <memory>
#include <string>
#include "ITaylorCoeff.h"

namespace alg {

class IFactoryMethodTaylorCoeff {
public:
virtual std::unique_ptr<ITaylorCoeff> factoryMethod(const std::string& typ_metody_obliczen,
                                                    long double rmin, long double rmax, 
                                                    long double horizon) = 0;
virtual ~IFactoryMethodTaylorCoeff(){};
};

}