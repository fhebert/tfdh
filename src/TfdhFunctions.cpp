
#include "TfdhFunctions.h"

#include "Element.h"
#include "GslWrappers.h"
#include "IntegrateOverRadius.h"
#include "PhysicalConstants.h"
#include "PlasmaFunctions.h"
#include "PlasmaState.h"
#include "TfdhSolution.h"

#include <cmath>


namespace {
  class EnergyDiff : public GSL::FunctionObject {
    private:
      const TfdhSolution& tfdh;
      const double kt;
      const int zion;
    public:
      EnergyDiff(const TfdhSolution& tfdh, const double kt, const int zion)
        : tfdh(tfdh), kt(kt), zion(zion) {}
      double operator()(const double r) const override {
        return zion * tfdh(r) - kt;
      }
  };
} // helper namespace



double TFDH::boundElectrons(const TfdhSolution& tfdh, const PlasmaState& p, const double cutoff)
{
  const auto f_ne_bound = [&] (const double r) -> double {return Plasma::neBound(tfdh(r), p, cutoff);};
  return integrateOverRadius(f_ne_bound, tfdh.r.front(), tfdh.r.back());
}


std::vector<double> TFDH::exclusionRadii(const TfdhSolution& tfdh,
    const Element& e, const PlasmaState& p)
{
  // for each ion species, find radius where:
  //   E_thermal == E_electrostatic  =>  kt == Zion phi
  std::vector<double> rexcl(p.ni);
  for (size_t elem=0; elem<rexcl.size(); ++elem) {
    const EnergyDiff delta(tfdh, p.kt, p.comp.species[elem].element.Z);
    const double eps_abs = 1e-6 * Plasma::radiusWignerSeitz(e, p);
    const double eps_rel = 1e-6;
    rexcl[elem] = GSL::findRoot(delta, tfdh.r.front(), tfdh.r.back(), eps_abs, eps_rel);
  }
  return rexcl;
}


TFDH::EnergyDeltas TFDH::embeddingEnergy(const TfdhSolution& tfdh,
    const Element& e, const PlasmaState& p)
{
  const auto f_dfi = [&] (const double r) -> double {return tfdh(r) * Plasma::totalIonChargeDensity(tfdh(r), p);};
  const double fi = integrateOverRadius(f_dfi, tfdh.r.front(), tfdh.r.back());

  const auto f_dfe = [&] (const double r) -> double {return - tfdh(r) * Plasma::ne(tfdh(r), p);};
  const double fe = integrateOverRadius(f_dfe, tfdh.r.front(), tfdh.r.back());

  const auto f_fce = [&] (const double r) -> double {
    const double& qe = PhysicalConstantsCGS::ElectronCharge;
    const double phi = tfdh(r);
    const double phi_ext = e.A * qe * qe / r;
    return 0.5 * (Plasma::totalIonChargeDensity(phi, p) - Plasma::ne(phi, p)) * (phi_ext - phi);
  };
  const double f2 = integrateOverRadius(f_fce, tfdh.r.front(), tfdh.r.back());

  const auto f_dni = [&] (const double r) -> double {
    double dni = 0;
    for (size_t i=0; i<p.ni.size(); ++i) {
      dni += Plasma::ni(tfdh(r), p)[i] - p.ni[i];
    }
    return dni;
  };
  const double ki = 1.5 * p.kt * integrateOverRadius(f_dni, tfdh.r.front(), tfdh.r.back());

  const auto f_dke = [&] (const double r) -> double {
    return Plasma::electronKineticEnergyDensity(tfdh(r), p) - Plasma::electronKineticEnergyDensity(0.0, p);};
  const double ke = integrateOverRadius(f_dke, tfdh.r.front(), tfdh.r.back());

  // TODO: are these even physically motivated?
  const double dni = ki;
  const auto f_dne = [&] (const double r) -> double {return Plasma::ne(tfdh(r), p) - p.ne;};
  const double dne = (Plasma::electronKineticEnergyDensity(0.0,p) / p.ne) * integrateOverRadius(f_dne, tfdh.r.front(), tfdh.r.back());

  return {fi, fe, f2, ki, ke, dni, dne, fi+fe+f2+ki+ke-dni-dne};
}


