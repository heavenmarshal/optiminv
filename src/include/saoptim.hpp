#ifndef __SAOPTIM_HPP__
#define __SAOPTIM_HPP__
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

#endif
