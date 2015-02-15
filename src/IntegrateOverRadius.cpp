
#include "IntegrateOverRadius.h"

#include "GslWrappers.h"

#include <cassert>
#include <cmath>


namespace {

  std::vector<double> weighByVolumeElement(const std::vector<double>& r,
      const std::vector<double>& f)
  {
    std::vector<double> weighted(f);
    for (size_t pt=0; pt<weighted.size(); ++pt) {
      weighted[pt] *= 4.0*M_PI * r[pt] * r[pt];
    }
    return weighted;
  }

  class RadialInterpFunction : public GSL::FunctionObject {
    private:
      GSL::Spline spline;
    public:
      RadialInterpFunction(const std::vector<double>& r, const std::vector<double>& f) : spline(r, f) {}
      double operator()(double x) const override {return spline.eval(x);}
  };
} // helper namespace


double integrateOverRadius(const std::vector<double>& r, const std::vector<double>& data) {
  assert(r.size()==data.size());

  // multiply data values by volume element (4 pi r^2) before integration
  const std::vector<double> integrand = weighByVolumeElement(r, data);

  // construct a GSL cubic spline approximation to the data
  // TODO: see if GSL interfacing can be cleaned up
  RadialInterpFunction interp_integrand(r, integrand);

  // perform integration with GSL calls
  // TODO:
  //   consider other algorithm possibilities:
  //   1- use the spline's exact integral
  //   2- use a spectral type approach
  //   3- ?
  const double r_min = r.front();
  const double r_max = r.back();
  const double eps_abs = 1.0e-6;
  const double eps_rel = eps_abs;
  return GSL::integrate(interp_integrand, r_min, r_max, eps_abs, eps_rel);
}


// TODO: move this to a proper home
std::vector<double> accumulateOverRadius(const std::vector<double>& r, const std::vector<double>& data) {
  assert(r.size()==data.size());
  const std::vector<double> integrand = weighByVolumeElement(r, data);
  RadialInterpFunction interp_integrand(r, integrand);

  const double eps_abs = 1.0e-6;
  const double eps_rel = eps_abs;
  const double r_min = r[0];

  std::vector<double> cumm(r.size(), 0);
  for (size_t pt=1; pt<r.size(); ++pt) {
    const double r_max = r[pt];
    cumm[pt] = integrate(interp_integrand, r_min, r_max, eps_abs, eps_rel);
  }
  return cumm;
}

