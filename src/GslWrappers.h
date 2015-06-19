
#ifndef TFDH_GSL_WRAPPERS_H
#define TFDH_GSL_WRAPPERS_H

#include <gsl/gsl_spline.h>
#include <vector>


namespace GSL {

  // simple functor class to help interfacing with gsl_function.
  // derived classes representing a function f(x) should hold any needed
  // parameters as member variables, and operator()(double x) should evaluate
  // the function f(x).
  //
  // NOTE: the GSL wrappers below could in principle take C++11 lambdas instead
  //       of this custom function object, however i did not find an elegant
  //       way to call these lambdas using GSL's C interfaces...
  class FunctionObject {
    public:
      virtual ~FunctionObject() = default;
      virtual double operator()(double x) const = 0;
  };

  // simple wrapper class around GSL splines
  class Spline {
    private:
      gsl_interp_accel* acc;
      gsl_spline* spline;
    public:
      Spline(const std::vector<double>& x, const std::vector<double>& f);
      ~Spline();
      double eval(double r) const;
  };

  // find a root of func in the interval x1 to x2
  // one of (x1,x2) or (x2,x1) MUST bracket a root
  // implemented using GSL's Brent rootfinder method
  double findRoot(const GSL::FunctionObject& func, double x1, double x2,
      double eps_abs, double eps_rel);

  // definite integral of func from x1 to x2
  // implemented using GSL's adaptive Gauss quadrature
  double integrate(const GSL::FunctionObject& func, double x1, double x2,
      double eps_abs, double eps_rel);

}


#endif // TFDH_GSL_WRAPPERS_H
