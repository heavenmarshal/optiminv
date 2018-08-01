#ifndef __ANNEALFN_HPP__
#define __ANNEALFN_HPP__
#include<cmath>
#define E1 1.7182818

class annealFn{
public:
  virtual ~annealFn(){};
  virtual double annealing(int iter, saOptim* optim);
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
}
#endif
