#ifndef __LASVDGPEIFN_HPP__
#define __LASVDGPEIFN_HPP__
#include "lasvdgpfn.hpp"

class lasvdgpEiFn: public lasvdgpFn{
public:
  lasvdgpEiFn(int nparam_, unsigned int ndesign_, unsigned int tlen_, unsigned int n0_,
	      unsigned int nn_, unsigned int nfea_, unsigned int resvdThres_,
	      unsigned int every_, double frac_, double gstart_, double **design_,
	      double **resp, double kmin_, double *xi_);

  double evaluate(double *x);
private:
  double kmin, z2;
  double *xi;
};
#endif
