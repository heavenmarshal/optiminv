#ifndef __KERNELFN_HPP__
#define __KERNELFN_HPP__
#include<ctime>
#include<random>
#include"saoptim.hpp"

class kernelFn{
public:
  kernelFn(int nparam_): nparam(nparam_), generator(time(NULL)){};
  virtual ~kernelFn(){};
  virtual void propose(double *from, double *to, saOptim* optim){};
protected:
  int nparam;
  std::default_random_engine generator;
};

// normal proposal kernel
class normalKernel: public kernelFn{
public:
  normalKernel(int nparam_, double sigma_):
    kernelFn(nparam_), sigma(sigma_){};
  void propose(double *from, double *to, saOptim* optim)
  {
    int i;
    double scale = optim -> getT();
    scale *= sigma;
    std::normal_distribution<double> distribution(0.0,scale);
    for(i=0; i<nparam; ++i)
      to[i] = from[i] + distribution(generator);
  }
private:
  double sigma;
};
#endif
