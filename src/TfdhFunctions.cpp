
#include "TfdhFunctions.h"

#include "Element.h"
#include "GslWrappers.h"
#include "IntegrateOverRadius.h"
#include "PhysicalConstants.h"
#include "PlasmaFunctions.h"
#include "PlasmaState.h"
#include "RadialFunction.h"
#include "Species.h"

#include <cmath>
#include <iostream> // TODO: remove this once devel couts are removed


namespace {
  class EnergyDiff : public GSL::FunctionObject {
    private:
      const double kt;
      const int zion;
      const GSL::Spline& spline;
    public:
      EnergyDiff(const double kt, const int zion, const GSL::Spline& spline)
        : kt(kt), zion(zion), spline(spline) {}
      double operator()(const double r) const override {
        return zion * spline.eval(r) - kt;
      }
  };


  RadialFunction deltaIonFieldEnergy(const RadialFunction& tfdh, const PlasmaState& p) {
    std::vector<double> dife(tfdh.data.size(), 0);
    for (size_t i=0; i<dife.size(); ++i) {
      const double phi = tfdh.data[i];
      dife[i] = phi * Plasma::totalIonChargeDensity(phi, p);
    }
    return RadialFunction(tfdh.radii, dife);
  }

  RadialFunction deltaElectronFieldEnergy(const RadialFunction& tfdh, const PlasmaState& p) {
    std::vector<double> defe(tfdh.data.size(), 0);
    for (size_t i=0; i<defe.size(); ++i) {
      const double phi = tfdh.data[i];
      defe[i] = (-1) * phi * Plasma::ne(phi, p);
    }
    return RadialFunction(tfdh.radii, defe);
  }

  RadialFunction deltaFieldCountingEnergy(const RadialFunction& tfdh, const PlasmaState& p,
      const Element& e) {
    std::vector<double> dfce(tfdh.data.size(), 0);
    const double& qe = PhysicalConstantsCGS::ElectronCharge;
    for (size_t i=0; i<dfce.size(); ++i) {
      const double phi = tfdh.data[i];
      const double phi_ext = e.Z * qe * qe / tfdh.radii[i];
      dfce[i] = 0.5 * (Plasma::totalIonChargeDensity(phi, p) - Plasma::ne(phi, p)) * (phi_ext - phi);
    }
    return RadialFunction(tfdh.radii, dfce);
  }

  double deltaNumberIons(const RadialFunction& tfdh, const PlasmaState& p) {
    double result = 0;
    std::vector<double> deltaNi(tfdh.data.size(), 0);
    // TODO clean this up
    for (size_t i=0; i<p.ni.size(); ++i) {
      for (size_t pt=0; pt<deltaNi.size(); ++pt) {
        deltaNi[pt] = Plasma::ni(tfdh.data[pt], p)[i] - p.ni[i];
      }
      result += integrateOverRadius(RadialFunction(tfdh.radii, deltaNi));
    }
    return result;
  }

  RadialFunction deltaElectronKineticEnergy(const RadialFunction& tfdh, const PlasmaState& p) {
    std::vector<double> deke(tfdh.data.size(), 0);
    for (size_t i=0; i<deke.size(); ++i)
      deke[i] = Plasma::neKinetic(tfdh.data[i], p) - Plasma::neKinetic(0.0, p);
    return RadialFunction(tfdh.radii, deke);
  }

  double deltaNumberElectrons(const RadialFunction& tfdh, const PlasmaState& p) {
    std::vector<double> deltaNe(tfdh.data.size(), 0);
    for (size_t i=0; i<deltaNe.size(); ++i) {
      deltaNe[i] = Plasma::ne(tfdh.data[i], p) - p.ne;
    }
    return integrateOverRadius(RadialFunction(tfdh.radii, deltaNe));
  }

} // helper namespace



RadialFunction TFDH::boundElectronDensity(const RadialFunction& tfdh, const PlasmaState& p, const double cutoff)
{
  std::vector<double> nebs(tfdh.data.size());
  for (size_t i=0; i<nebs.size(); ++i) {
    nebs[i] = Plasma::neBound(tfdh.data[i], p, cutoff);
  }
  return RadialFunction(tfdh.radii, nebs);
}


double TFDH::boundElectrons(const RadialFunction& tfdh, const PlasmaState& p, const double cutoff)
{
  std::vector<double> nebs(tfdh.data.size());
  for (size_t i=0; i<nebs.size(); ++i) {
    nebs[i] = Plasma::neBound(tfdh.data[i], p, cutoff);
  }
  return integrateOverRadius(RadialFunction(tfdh.radii, nebs));
}


std::vector<double> TFDH::exclusionRadii(const RadialFunction& tfdh,
    const Element& e, const PlasmaState& p)
{
  // for each ion species, find radius where:
  //   E_thermal == E_electrostatic  =>  kt == Zion phi
  const GSL::Spline spline(tfdh);
  std::vector<double> rexcl(p.ni);
  for (size_t elem=0; elem<rexcl.size(); ++elem) {
    const EnergyDiff delta(p.kt, p.comp.species[elem].element.Z, spline);
    const double eps_abs = 1e-6 * Plasma::radiusWignerSeitz(e, p);
    const double eps_rel = 1e-6;
    rexcl[elem] = GSL::findRoot(delta, tfdh.radii.front(), tfdh.radii.back(), eps_abs, eps_rel);
  }
  return rexcl;
}


double TFDH::embeddingEnergy(const RadialFunction& tfdh, const PlasmaState& p, const Element& e)
{
  const double fi = integrateOverRadius(deltaIonFieldEnergy(tfdh, p));
  const double fe = integrateOverRadius(deltaElectronFieldEnergy(tfdh, p));
  const double f2 = integrateOverRadius(deltaFieldCountingEnergy(tfdh, p, e));
  const double ki = 1.5 * p.kt * deltaNumberIons(tfdh, p);
  const double ke = integrateOverRadius(deltaElectronKineticEnergy(tfdh, p));

  // TODO: are these even physically motivated?
  const double dni = ki;
  const double dne = (Plasma::neKinetic(0.0,p) / p.ne) * deltaNumberElectrons(tfdh, p);

  // TODO: remove these devel couts. information should be dumped to file instead
  std::cout << "\nembedding energies in units of kT:\n";
  std::cout << "ion field energy:             " << fi/p.kt << "\n";
  std::cout << "e- field energy:              " << fe/p.kt << "\n";
  std::cout << "overcounting of field energy: " << f2/p.kt << "\n";
  std::cout << "change in ion kinetic energy: " << ki/p.kt << "\n";
  std::cout << "change in e- kinetic energy:  " << ke/p.kt << "\n";
  std::cout << "energy from exchanging ions:  " << dni/p.kt << "\n";
  std::cout << "energy from exchanging e-'s:  " << dne/p.kt << "\n";
  std::cout << "sum of everything        : " << (fi+fe+f2+ki+ke-dni-dne)/p.kt << "\n";
  std::cout << "sum of the ones i like.. : " << (fi+fe+f2+ki+ke)/p.kt << "\n";
  return fi+f2+fe+ki+ke;
}


