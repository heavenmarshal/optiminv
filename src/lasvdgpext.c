#include <stdlib.h>
#include "lasvdgp.h"
#include "gp_sep.h"
#include "matrix.h"
#include "linalg.h"

void predlasvdGPext(lasvdGP *lasvdgp, double *x, double *cmean,
		    double *cs2, double *sig2hat)
{
  int i, n0, tlen, nbas, reslen;
  double cdf;
  double **resid, **coeff;
  GPsep **gpseps;
  gpseps = lasvdgp->gpseps;
  n0 = lasvdgp -> n0;
  tlen = lasvdgp -> tlen;
  nbas = lasvdgp -> nbas;
  coeff = new_zero_matrix(nbas, n0);

  for(i=0; i<nbas; ++i)
  {
    linalg_daxpy(n0, lasvdgp->reds[i], gpseps[i]->Z, 1, coeff[i], 1);
    predGPsep_lite(gpseps[i], 1, &x, cmean+i, cs2+i, &cdf, NULL);
  }
  resid = new_p_submatrix_rows(lasvdgp->feaidx, lasvdgp->resp, n0, tlen, 0);
  linalg_dgemm(CblasNoTrans, CblasTrans, tlen, n0, nbas, -1.0, &(lasvdgp->basis), tlen,
	       coeff, n0, 1.0, resid, tlen);
  reslen = n0*tlen;
  *sig2hat = linalg_ddot(reslen,*resid,1,*resid,1);
  *sig2hat /= (reslen+2);
  delete_matrix(coeff);
  delete_matrix(resid);
}

void lasvdgpEIhelp(lasvdGP *lasvdgp, double* x, double *xi, double *varres,
		   double *iomemu2, double *bound, double *amat, double *mub2star,
		   double *mumk)
{
  int i, tlen, nbas;
  double *cmean, *cs2, *bz;
  double d, d2, sand, bilform, aval, maxval;
  nbas = lasvdgp->nbas;
  tlen = lasvdgp->tlen;
  cmean = new_vector(nbas);
  cs2 = new_vector(nbas);
  bz = new_vector(nbas);
  predlasvdGPext(lasvdgp, x, cmean, cs2, varres);
  linalg_dgemv(CblasTrans, tlen, nbas, 1.0, &(lasvdgp->basis), tlen, xi, 1,
	       0.0, bz, 1);
  *iomemu2 = 0.0;
  maxval = -1.0;
  for(i=0; i<nbas; ++i)
  {
    d = lasvdgp->reds[i];
    d2 = sq(d);
    aval = cs2[i]*d2+ (*varres);
    if(aval > maxval) maxval = aval;
    sand = cs2[i]/aval;
    bilform = sq(cmean[i])*d2-2.0*cmean[i]*bz[i]*d;
    *iomemu2 +=  bilform - d2*sand*sq(bz[i]);
    *iomemu2 += 2.0*sand*d2*d*cmean[i]*bz[i] - sq(cmean[i])*sq(d2)*sand;
    mub2star[i] = sand * sq(bz[i]*d-d2*cmean[i]);
    amat[i] = aval;
    *mumk += bilform + d2*cs2[i];
  }
  *iomemu2 /= (*varres);
  *bound = 0.5/maxval;
  *mumk += (*varres)*(double)tlen;
  free(cmean);
  free(cs2);
  free(bz);
}
