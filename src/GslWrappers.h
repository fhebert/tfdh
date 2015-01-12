
#ifndef TFDH_GSL_WRAPPERS_H
#define TFDH_GSL_WRAPPERS_H

#include <gsl/gsl_spline.h>

struct RadialFunction;


namespace GSL {

  // simple wrapper around GSL splines
  class Spline {
    private:
      gsl_interp_accel* acc;
      gsl_spline* spline;
    public:
      Spline(const RadialFunction& data);
      ~Spline();
      double eval(double r) const;
  };


  // simple base class to help interfacing with gsl_function.
  //
  // derived classes should contain any necessary parameters as member
  // variables, and operator() should evaluate the function at the given
  // value of the parameter.
  //
  // NOTE: this is kind of like a hacky lambda. this is used instead of a
  // real c++11 lambda because i couldn't find an elegant way to mesh such
  // a lambda with the gsl_function interface.
  class FunctionObject {
    public:
      virtual ~FunctionObject() = default;
      virtual double operator()(double x) const = 0;
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
