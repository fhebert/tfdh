
#ifndef TFDH_RADIAL_INTEGRAL_H
#define TFDH_RADIAL_INTEGRAL_H

#include "GslWrappers.h"

#include <cmath>
#include <vector>


template <typename T>
double integrateOverRadius(const T& func, const double rmin, const double rmax)
{
  class MultiplyByJacobian : public GSL::FunctionObject {
    private:
      const T& f;
    public:
      MultiplyByJacobian(const T& func) : f(func) {}
      double operator()(const double r) const override {return 4*M_PI*r*r*f(r);};
  };
  const auto integrand = MultiplyByJacobian(func);
  const double eps = 1.e-6;
  const double eps_abs = eps * (rmax-rmin) * (fabs(integrand(rmax))+fabs(integrand(rmin))) / 2;
  const double eps_rel = eps;
  return GSL::integrate(integrand, rmin, rmax, eps_abs, eps_rel);
}


#endif // TFDH_RADIAL_INTEGRAL_H

