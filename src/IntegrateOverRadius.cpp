
#include "IntegrateOverRadius.h"

#include "GslWrappers.h"
#include "RadialFunction.h"

#include <cmath>
#include <gsl/gsl_spline.h>


// helpers
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
      gsl_interp_accel* acc;
      gsl_spline* spline;
    public:
      RadialInterpFunction(const RadialFunction& f) {
        acc = gsl_interp_accel_alloc();
        spline = gsl_spline_alloc(gsl_interp_cspline, f.radii.size());
        gsl_spline_init(spline, f.radii.data(), f.data.data(), f.radii.size());
      }
      ~RadialInterpFunction() {
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
      }
      double operator()(double x) const {
        return gsl_spline_eval(spline, x, acc);
      }
  };
}



double integrateOverRadius(const RadialFunction& input) {

  // multiply data values by volume element (4 pi r^2) before integration
  const RadialFunction integrand = weighByVolumeElement(input);

  // construct a GSL cubic spline approximation to the data
  // TODO: see if GSL intefacing can be cleaned up
  RadialInterpFunction interp_integrand(integrand);

  // perform integration with GSL calls
  // use simplest non-adaptive GSL algorithm until a need for more arises
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

