
#ifndef TFDH_GSL_WRAPPERS_H
#define TFDH_GSL_WRAPPERS_H


namespace GSL {

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
      virtual double operator()(double x) const = 0;
  };

  // find a root of func in the interval x1 to x2
  // interval (x1,x2) MUST bracket a root and x1<x2
  // implemented using GSL's Brent rootfinder method
  double findRoot(const GSL::FunctionObject& func, double x1, double x2,
      double dx_abs);

  // definite integral of func from x1 to x2
  // implemented using GSL's adaptive Gauss quadrature
  double integrate(const GSL::FunctionObject& func, double x1, double x2,
      double eps_abs, double eps_rel);

}


#endif // TFDH_GSL_WRAPPERS_H
