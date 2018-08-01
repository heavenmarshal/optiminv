#ifndef __ANNEALFN_HPP__
#define __ANNEALFN_HPP__
#include "saoptim.hpp"

class annealFn{
public:
  virtual ~annealFn(){};
  virtual double annealing(int iter, saOptim* optim){return 0.0;};
};

// the default annealing function in R
class nlogAnneal: public annealFn{
public:
  nlogAnneal(double ti_, double offset_=E1):
    ti(ti_), offset(offset_){};
  double annealing(int iter, saOptim* optim)
  {
    return ti/log((double)iter + offset);
  }
private:
  double ti, offset;
};
#endif
