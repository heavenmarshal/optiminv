#include <cstdlib>
#include "lasvdgpfn.hpp"

extern "C"{
  #include "matrix.h"
  #include "linalg.h"
  #include "lasvdgp.h"
  #include "gp_sep.h"
}

static void vec_divide(double *vio, double* vdivid, int len)
{
  int i;
  for(i = 0; i < len; ++i)
    vio[i] /= vdivid[i];
}
void lasvdgpFn::fitlasvdGP(double *param)
{
  lasvdgp = newlasvdGP(param, design, resp, ndesign, nparam, tlen, nn, n0,
		       nfea, nn, 1, frac, gstart);
  jmlelasvdGP(lasvdgp, 100, 0);
  iterlasvdGP(lasvdgp, resvdThres, every, 100, 0);
}
void lasvdgpFn::rmlasvdGP()
{
  if(lasvdgp) deletelasvdGP(lasvdgp);
  lasvdgp = NULL;
}

double lasvdgpEsl2dFn::evaluate(double *x)
{
  int nbas, i;
  double cmean, cs2, cdf;
  double esl2d, *hatcxi;

  fitlasvdGP(x);		// todo: do not fit model at every input, reuse model
  nbas = lasvdgp->nbas;
  hatcxi = new_vector(nbas);
  linalg_dgemv(CblasTrans, tlen, nbas, 1.0, &(lasvdgp->basis), tlen, xi, 1,
	       0.0, hatcxi, 1);
  vec_divide(hatcxi, lasvdgp->reds, nbas);
  esl2d = 0.0;
  for(i = 0; i < nbas; ++i)
  {
    predGPsep_lite(lasvdgp->gpseps[i], 1, &x, &cmean, &cs2, &cdf,NULL);
    esl2d += sq(lasvdgp->reds[i])*(sq(cmean - hatcxi[i])+cs2);
  }
  free(hatcxi);
  rmlasvdGP();
  return esl2d;
}
