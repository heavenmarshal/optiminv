#ifndef __SAOPTIM_HPP__
#define __SAOPTIM_HPP__
#include<ctime>
#include<random>
#include"fminfn.hpp"
#define E1 1.7182818
class kernelFn;
class annealFn;

class saOptim{
public:
  saOptim(int nparam_, int maxit_, int tmax_, double *xstart,
	  fminFn* mfun, annealFn *afun, kernelFn *kfun);
  ~saOptim();
  virtual void optim();
  void report(double* min, double *argmin,
	      double *fvals, double *samples);		// report the result for R iterface;
  double getT(){return temperature;};
private:
  int nparam;	// dimension of parameters
  int its;	// current number of iterations
  int maxit;	// maximum total number of iterations
  int tmax;			// maximum number of iteration at the same temperature
  double cmin;			// current minimum value
  double temperature;		// current temperature
  double *cargmin;		// current argmin
  double *fvaltr;		// history of min values
  double **sampletr;		// history of argmins
  fminFn *fmin;			// function (object) to be minimized
  annealFn *anneal;		// annealing function
  kernelFn *kernel;		// proposal kernel function
  std::default_random_engine generator;
};

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
};
class kernelFn{
public:
  kernelFn(int nparam_): nparam(nparam_), generator(time(NULL)){};
  virtual ~kernelFn(){};
  virtual void propose(double *from, double *to, saOptim* optim);
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
