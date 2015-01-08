
#include "TfdhFunctions.h"

#include "GslWrappers.h"
#include "PlasmaFunctions.h"
#include "PlasmaState.h"
#include "PhysicalConstants.h"
#include "RadialFunction.h"

#include <cmath>


namespace {

  class FermiDiracDistribution : public GSL::FunctionObject {
    private:
      const double xi;
      const PlasmaState p;
    public:
      FermiDiracDistribution(const double xi, const PlasmaState& p) : xi(xi), p(p) {}
      double operator()(const double x) const {
        const double val = (1.0 + p.tau*x) * sqrt(x + p.tau * x*x/2.0) / (1.0 + exp(x-p.chi-xi));
        return val;
      }
  };

  double boundFermiDiracIntegral(const double xi, const PlasmaState& p) {
    FermiDiracDistribution fd(xi, p);
    const double eps = 1.e-6;
    return GSL::integrate(fd, 0, xi, eps, eps);
  }

}


double TFDH::neBound(const double phi, const PlasmaState& p)
{
  const double xi = phi/p.kt;
  if (xi < 0) {
    return 0;
  } else {
    const double& NePrefactor = PhysicalConstantsCGS::NePrefactor;
    return NePrefactor * pow(p.kt, 1.5) * boundFermiDiracIntegral(xi, p);
  }
}


RadialFunction TFDH::neBound(const RadialFunction& tfdh, const PlasmaState& p)
{
  std::vector<double> nebs(tfdh.data.size());
  for (size_t i=0; i<tfdh.radii.size(); ++i) {
    nebs[i] = neBound(tfdh.data[i], p);
  }
  return RadialFunction(tfdh.radii, nebs);
}

