
#include "IntegrateOverRadius.h"

#include "GslWrappers.h"
#include "RadialFunction.h"

#include <cmath>


namespace {

  RadialFunction weighByVolumeElement(const RadialFunction& f) {
    std::vector<double> weighted(f.data);
    for (size_t r=0; r<weighted.size(); ++r) {
      weighted[r] *= 4.0*M_PI * f.radii[r] * f.radii[r];
    }
    return RadialFunction(f.radii, weighted);
  }

  class RadialInterpFunction : public GSL::FunctionObject {
    private:
      GSL::Spline spline;
    public:
      RadialInterpFunction(const RadialFunction& f) : spline(f) {}
      double operator()(double x) const {return spline.eval(x);}
  };
} // helper namespace


double integrateOverRadius(const RadialFunction& input) {

  // multiply data values by volume element (4 pi r^2) before integration
  const RadialFunction integrand = weighByVolumeElement(input);

  // construct a GSL cubic spline approximation to the data
  // TODO: see if GSL interfacing can be cleaned up
  RadialInterpFunction interp_integrand(integrand);

  // perform integration with GSL calls
  // TODO:
  //   consider other algorithm possibilities:
  //   1- use the spline's exact integral
  //   2- use a spectral type approach
  //   3- ?
  const double r_min = integrand.radii.front();
  const double r_max = integrand.radii.back();
  const double eps_abs = 1.0e-6;
  const double eps_rel = eps_abs;
  return GSL::integrate(interp_integrand, r_min, r_max, eps_abs, eps_rel);
}


// TODO: move this to a proper home
RadialFunction accumulateOverRadius(const RadialFunction& input) {
  const RadialFunction integrand = weighByVolumeElement(input);
  RadialInterpFunction interp_integrand(integrand);

  const double eps_abs = 1.0e-6;
  const double eps_rel = eps_abs;
  const double r_min = integrand.radii[0];

  std::vector<double> cumm(input.data.size(), 0);
  for (size_t i=1; i<cumm.size(); ++i) {
    const double r_max = integrand.radii[i];
    cumm[i] = integrate(interp_integrand, r_min, r_max, eps_abs, eps_rel);
  }
  return RadialFunction(input.radii, cumm);
}

