#ifndef __FMINFN_HPP__
#define __FMINFN_HPP__
class fminFn{
public:
  fminFn(int nparam_): nparam(nparam_){};
  virtual ~fminFn(){};
  virtual double evaluate(double *x){return 0.0;};
  virtual void gradient(double* x, double* grad){};
  virtual void hessian(double* x, double** hessian){};
protected:
  int nparam;
};
#endif
