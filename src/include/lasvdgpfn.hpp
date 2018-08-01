#ifndef __LASVDGPFN_HPP__
#define __LASVDGPFN_HPP__
#include "fminfn.hpp"
extern "C"{
  #include "matrix.h"
  #include "lasvdgp.h"
}
class lasvdgpFn: public fminFn{
public:
  lasvdgpFn(int nparam_, unsigned int ndesign_, unsigned int tlen_, unsigned int n0_,
	    unsigned int nn_, unsigned int nfea_, unsigned int resvdThres_,
	    unsigned int every_, double frac_, double gstart_, double **design_,
	    double **resp_):
    fminFn(nparam_), ndesign(ndesign_), tlen(tlen_), n0(n0_), nn(nn_), nfea(nfea_),
    resvdThres(resvdThres_), every(every_), frac(frac_), gstart(gstart_), design(design_),
    resp(resp_), lasvdgp(NULL){};
  virtual ~lasvdgpFn()
  {
    rmlasvdGP();
  }
protected:
  unsigned int ndesign, tlen, n0, nn, nfea;
  unsigned int resvdThres, every;
  double frac, gstart;
  double **design, **resp;
  lasvdGP *lasvdgp;
  void fitlasvdGP(double *param);
  void rmlasvdGP();
};

class lasvdgpEsl2dFn: public lasvdgpFn{
public:
  lasvdgpEsl2dFn(int nparam_, unsigned int ndesign_, unsigned int tlen_, unsigned int n0_,
		 unsigned int nn_, unsigned int nfea_, unsigned int resvdThres_,
		 unsigned int every_, double frac_, double gstart_, double **design_,
		 double **resp_, double *xi_):
    lasvdgpFn(nparam_, ndesign_, tlen_, n0_,  nn_, nfea_, resvdThres_,
	      every_, frac_, gstart_, design_,  resp_),xi(xi_){};

  double evaluate(double *x);
private:
  double *xi;
};

#endif
