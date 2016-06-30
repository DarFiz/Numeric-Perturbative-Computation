#pragma once

class IPsi {
public:
  virtual ~IPsi(){};

  virtual void setPsiElem(int j, long double wartoscPsi) = 0;
  virtual long double getPsiElem(int j) = 0;
  virtual long double getDPsiElem(int j, int rzadpoch) = 0;
  virtual int getMaxj() = 0;
  virtual long double geth() = 0;
};
