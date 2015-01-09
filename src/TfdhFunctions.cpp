
#include "TfdhFunctions.h"

#include "Element.h"
#include "GslWrappers.h"
#include "IntegrateOverRadius.h"
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

  class EnergyDiff : public GSL::FunctionObject {
    private:
      const double kt;
      const int zion;
      const GSL::Spline& spline;
    public:
      EnergyDiff(const double kt, const int zion, const GSL::Spline& spline)
        : kt(kt), zion(zion), spline(spline) {}
      double operator()(const double r) const {
        return zion * spline.eval(r) - kt;
      }
  };
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


double TFDH::boundElectrons(const RadialFunction& tfdh, const PlasmaState& p)
{
  const RadialFunction& neBound = TFDH::neBound(tfdh, p);
  return integrateOverRadius(neBound);
}


std::vector<double> TFDH::exclusionRadii(const RadialFunction& tfdh,
    const Element& e, const PlasmaState& p)
{
  // for each ion species, find radius where:
  //   E_thermal = E_electrostatic  =>  kt = Zion phi
  const GSL::Spline spline(tfdh);
  std::vector<double> result(p.ni);
  for (size_t elem=0; elem<result.size(); ++elem) {
    const EnergyDiff delta(p.kt, p.comp.abundances[elem].element.Z, spline);
    const double eps_abs = 1e-6 * Plasma::radiusWignerSeitz(e, p);
    const double eps_rel = 1e-6;
    result[elem] = GSL::findRoot(delta, tfdh.radii.front(), tfdh.radii.back(), eps_abs, eps_rel);
  }
  return result;
}
