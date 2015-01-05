
#ifndef TFDH_GSL_WRAPPERS_H
#define TFDH_GSL_WRAPPERS_H


namespace GSL {

  // simple base class that wraps the gsl_function interface.
  // additional parameters to the function are to be held in
  // member variables of the derived classes.
  class FunctionObject {
    public:
      virtual ~FunctionObject() = default;
      virtual double f(double x) const = 0;
  };

  // helper function
  // TODO: documentation for the cleverness
  double callFunctionFromObject(double x, void *params);


  // object to wrap a GSL ODE integrator
  // implements an adaptive stepper
  //class OdeIntegrator {
  //  public:
  //    virtual ~OdeIntegrator() = default; // ?
  //    virtual int step() const = 0;
  //};

  // find a real foot of func in the interval x1 to x2
  //   interval x1-x2 must bracket a root!
  // implemented using GSL's Brent rootfinder method
  double findRoot(const GSL::FunctionObject& func, double x1, double x2,
      double dx_abs);

  // definite integral of func from x1 to x2
  // implemented using GSL's adaptive Gauss quadrature
  double integrate(const GSL::FunctionObject& func, double x1, double x2,
      double eps_abs, double eps_rel);

}


#endif // TFDH_GSL_WRAPPERS_H
