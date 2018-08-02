#include "fminfn.hpp"
#include "saoptim.hpp"
#include "kernelfn.hpp"
#include "annealfn.hpp"
#include "lasvdgpfn.hpp"
#include "lasvdgpeifn.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
extern "C"{
  #include "matrix.h"
  #include "linalg.h"
  #include "rhelp.h"
  #include "lasvdgp.h"
  #include "lasvdgpext.h"
  #include "oeiinfo.h"
}

extern "C"{
  void lasvdgpEsl2d(int *nparam_, int *maxit_, int *tmax_, unsigned int *ndesign_,
		    unsigned int *tlen_, unsigned int *n0_, unsigned int *nn_,
		    unsigned int *nfea_, unsigned int *resvdThres_, unsigned int *every_,
		    unsigned int *nthread_, unsigned int* nstarts_,double *frac_,
		    double *gstart_, double *ti_, double *kersig_, double *xstarts_,
		    double *xi_, double *design_, double *resp_, double *mins,
		    double *argmins, double *fvals, double *samples)
  {
    int nsamples, slen;
    unsigned int mxth, ndesign, nparam, tlen;
    double **design, **resp, **xstarts;
    ndesign = *ndesign_;
    nparam = *nparam_;
    tlen = *tlen_;
    design = new_matrix_bones(design_, ndesign, nparam);
    resp = new_matrix_bones(resp_, ndesign, tlen);
    xstarts = new_matrix_bones(xstarts_, *nstarts_, nparam);
    nsamples = (*maxit_)+1;
    slen = nparam * nsamples;
#ifdef _OPENMP
    mxth = omp_get_max_threads();
#else
    mxth = 1;
#endif
    if(*nthread_>mxth)
    {
      MYprintf(MYstdout, "NOTE: omp.threads(%d) > max(%d), using %d\n",
	       *nthread_, mxth, mxth);
      *nthread_ = mxth;
    }
#ifdef _OPENMP
#pragma omp parallel num_threads(*nthread_)
    {
      unsigned int i, start, step;
      start = omp_get_thread_num();
      step = *nthread_;
#else
      unsigned int i, start, step;
      start = 0; step = 1;
#endif
      for(i = start; i < *nstarts_; i += step)
      {
	lasvdgpEsl2dFn esl2dfn(nparam, ndesign, tlen, *n0_, *nn_,
			       *nfea_, *resvdThres_, *every_, *frac_,
			       *gstart_, design, resp, xi_);
	nlogAnneal anneal(*ti_);
	normalKernel kernel(nparam, *kersig_);
	saOptim saoptim(nparam, *maxit_, *tmax_, xstarts[i], &esl2dfn,
			&anneal, &kernel);
	saoptim.optim();
	saoptim.report(mins+i, argmins+i*nparam,
		       fvals+i*nsamples, samples+i*slen);
      }
#ifdef _OPENMP
    }
#endif
    free(xstarts);
    free(design);
    free(resp);
  }
  // optimize EI using the simulated annealing algorithm
  void lasvdgpEI(int *nparam_, int *maxit_, int *tmax_, unsigned int *ndesign_,
		 unsigned int *tlen_, unsigned int *n0_, unsigned int *nn_,
		 unsigned int *nfea_, unsigned int *resvdThres_, unsigned int *every_,
		 double *frac_, double *gstart_, double *ti_, double *kersig_,
		 double *xstart_, double *kmin_, double *xi_, double *design_,
		 double *resp_, double *mins, double *argmins, double *fvals,
		 double *samples)
  {
    int nparam;
    unsigned int ndesign, tlen;
    double **design, **resp;
    ndesign = *ndesign_;
    nparam = *nparam_;
    tlen = *tlen_;
    design = new_matrix_bones(design_, ndesign, nparam);
    resp = new_matrix_bones(resp_, ndesign, tlen);

    lasvdgpEiFn eifn(nparam, ndesign, tlen, *n0_, *nn_,
		     *nfea_, *resvdThres_, *every_, *frac_,
		     *gstart_, design, resp, *kmin_, xi_);
    nlogAnneal anneal(*ti_);
    normalKernel kernel(nparam, *kersig_);
    saOptim saoptim(nparam, *maxit_, *tmax_, xstart_, &eifn,
		    &anneal, &kernel);
    saoptim.optim();
    saoptim.report(mins, argmins, fvals, samples);
    free(design);
    free(resp);
  }
  // evaluate EI on a candidate set
  void evalLasvdgpEI(unsigned int *nparam_, unsigned int *ndesign_, unsigned int *tlen_,
		     unsigned int *ncand_, unsigned int *n0_, unsigned int *nn_,
		     unsigned int *nfea_, unsigned int *resvdThres_, unsigned int *every_,
		     unsigned int *nthread_, double *frac_, double *gstart_, double *kmin_,
		     double *xi_, double *design_, double *resp_, double *candid_,
		     double *eival)
  {
    unsigned int nparam, ndesign, tlen, ncand, mxth;
    double **design, **resp, **candid;
    ndesign = *ndesign_;
    nparam = *nparam_;
    tlen = *tlen_;
    ncand = *ncand_;
    design = new_matrix_bones(design_, ndesign, nparam);
    resp = new_matrix_bones(resp_, ndesign, tlen);
    candid = new_matrix_bones(candid_, ncand, nparam);
#ifdef _OPENMP
    mxth = omp_get_max_threads();
#else
    mxth = 1;
#endif
    if(*nthread_ > mxth)
    {
      MYprintf(MYstdout, "NOTE: omp.threads(%d) > max(%d), using %d\n",
	       *nthread_, mxth, mxth);
      *nthread_ = mxth;
    }
#ifdef _OPENMP
#pragma omp parallel num_threads(*nthread_)
    {
      unsigned int i, start, step;
      start = omp_get_thread_num();
      step = *nthread_;
#else
      unsigned int i, start, step;
      start = 0; step = 1;
#endif
      double z2;
      lasvdGP* lasvdgp;
      z2 = linalg_ddot(tlen, xi_, 1, xi_, 1);
      for(i = start; i < ncand; i+=step)
      {
	lasvdgp = newlasvdGP(candid[i], design, resp, ndesign, nparam, tlen,
			     *nn_, *n0_, *nfea_, *nn_, 1, *frac_, *gstart_);
	jmlelasvdGP(lasvdgp, 100, 0);
	iterlasvdGP(lasvdgp, *resvdThres_, *every_, 100, 0);
#ifdef _OPENMP
#pragma omp critical
	{
#endif
	  eival[i] = evalEI(lasvdgp, candid[i], xi_, z2, *kmin_);
#ifdef _OPENMP
	}
#endif
	deletelasvdGP(lasvdgp);
      }
#ifdef _OPENMP
    }
#endif
    free(design);
    free(resp);
    free(candid);
  }
}
