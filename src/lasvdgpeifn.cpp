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
			 double **resp, double kmin_, double *xi_):
  lasvdgpFn(nparam_, ndesign_, tlen_, n0_,  nn_, nfea_, resvdThres_,
	    every_, frac_, gstart_, **design_,  **resp), kmin(kmin_), xi(xi_){
  z2 = linalg_ddot(tlen, xi, 1, xi, 1);
}

double lasvdgpEiFn::evaluate(double *x)
{
  int nbas, i;
  double ei, varres, iomemu2, bound, mumk;
  double *amat, *mub2star;

  fitlasvdGP(x);
  nbas = lasvdgp->nbas;
  amat = new_vector(nbas);
  mub2star = new_vector(nbas);
  lasvdgpEIhelp(lasvdgp, x, xi, &varres, &iomemu2, &bound, amat, mub2star, &mumk);
  iomemu2 += z2;
  mumk -= kmin;
  oeiinfo(1, nbas, tlen, varres, kmin, &iomemu2, &bound,
	  amat, mub2star, &mumk, &ei);
  ei -= mumk;
  rmlasvdGP();
  free(amat);
  free(mub2star);
  return -ei;			// minimize negative ei
}