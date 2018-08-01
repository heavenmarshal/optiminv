#include <cstdlib>
#include <ctime>
#include "saoptim.hpp"
#include "kernelfn.hpp"
#include "annealfn.hpp"
#include "fminfn.hpp"
extern "C"{
  #include "matrix.h"
}
saOptim::saOptim(int nparam_, int maxit_, int tmax_, double *xstart,
		 fminFn *mfun, annealFn *afun, kernelFn *kfun):
  nparam(nparam_), its(0), maxit(maxit_), tmax(tmax_), fmin(mfun),
  anneal(afun), kernel(kfun), generator(time(NULL))
{
  cargmin = new_dup_vector(xstart,nparam);
  cmin = fmin->evaluate(cargmin);
  fvaltr = new_vector(maxit+1);
  sampletr = new_matrix(maxit+1,nparam);
  fvaltr[0] = cmin;
  dupv(sampletr[0],cargmin,nparam);
}
saOptim::~saOptim()
{
  free(cargmin);
  free(fvaltr);
  delete_matrix(sampletr);
}

void saOptim::optim()
{
  int k;
  double y, ytry, dy, logru;
  double *x, *xtry;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  x = new_dup_vector(cargmin,nparam);
  xtry = new_vector(nparam);
  y = cmin;
  its = 0;

  while(its < maxit)
  {
    temperature = anneal->annealing(its+1,this);
    k = 0;
    while((k < tmax) && (its < maxit))
    {
      kernel->propose(x, xtry, this);
      ytry = fmin -> evaluate(xtry);
      // todo: check finite
      dy = ytry - y;
      logru = distribution(generator);
      logru = log(logru);
      if((dy <= 0.0) || (logru < -dy/temperature))
      {
	dupv(x, xtry, nparam);
	y = ytry;
	if(y < cmin)
	{
	  cmin = y;
	  dupv(cargmin, x, nparam);
	}
      }
      its++; k++;
      fvaltr[its] = y;
      dupv(sampletr[its], x, nparam);
    }
  }
  free(x);
  free(xtry);
}

void saOptim::report(double *min, double *argmin,
		     double *fvals, double *samples)
{
  int nsps = maxit + 1;
  *min = cmin;
  dupv(argmin, cargmin, nparam);
  dupv(fvals, fvaltr, nsps);
  dupv(samples, *sampletr, nsps*nparam);
}
