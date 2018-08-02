#include <cstdlib>
#include "lasvdgpeifn.hpp"
extern "C"{
  #include "oeiinfo.h"
  #include "lasvdgpext.h"
  #include "matrix.h"
  #include "linalg.h"
}

lasvdgpEiFn::lasvdgpEiFn(int nparam_, unsigned int ndesign_, unsigned int tlen_, unsigned int n0_,
			 unsigned int nn_, unsigned int nfea_, unsigned int resvdThres_,
			 unsigned int every_, double frac_, double gstart_, double **design_,
			 double **resp_, double kmin_, double *xi_):
  lasvdgpFn(nparam_, ndesign_, tlen_, n0_,  nn_, nfea_, resvdThres_,
	    every_, frac_, gstart_, design_,  resp_), kmin(kmin_), xi(xi_){
  z2 = linalg_ddot(tlen, xi, 1, xi, 1);
}

double lasvdgpEiFn::evaluate(double *x)
{
  double ei;
  fitlasvdGP(x);
  ei = evalEI(lasvdgp, x, xi, z2, kmin);
  rmlasvdGP();
  return -ei;			// minimize negative ei
}
