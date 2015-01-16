
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



double TFDH::boundElectrons(const RadialFunction& tfdh, const PlasmaState& p)
{
  std::vector<double> nebs(tfdh.data.size());
  for (size_t i=0; i<nebs.size(); ++i) {
    nebs[i] = Plasma::neBound(tfdh.data[i], p);
  }
  return integrateOverRadius(RadialFunction(tfdh.radii, nebs));
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
