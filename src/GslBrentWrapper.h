
#ifndef TFDH_GSL_BRENT_WRAPPER_H
#define TFDH_GSL_BRENT_WRAPPER_H


// simple base class that wraps the gsl_function interface.
// additional parameters to the function are to be held in
// member variables of the derived classes.
class GslFunction {
  public:
  virtual ~GslFunction() = default;
  virtual double f(double x) const = 0;
};

double GslBrent(const GslFunction& func, double x1, double x2, double dx_abs);


#endif // TFDH_GSL_BRENT_WRAPPER_H
