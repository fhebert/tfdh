
#ifndef TFDH_GSL_WRAPPERS_H
#define TFDH_GSL_WRAPPERS_H


// simple base class that wraps the gsl_function interface.
// additional parameters to the function are to be held in
// member variables of the derived classes.
class GslFunction {
  public:
  virtual ~GslFunction() = default;
  virtual double f(double x) const = 0;
};

// helper function
// TODO: documentation for the cleverness
double GslFunction_Unpacker(double x, void *params);



double gslBrent(const GslFunction& func, double x1, double x2, double dx_abs);

double gslQuadratureNG(const GslFunction& func, double x1, double x2,
    double eps_abs, double eps_rel);


#endif // TFDH_GSL_WRAPPERS_H
