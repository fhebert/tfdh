
#include "PlasmaFunctions.h"

#include "Composition.h"
#include "Element.h"
#include "Gfdi.h"
#include "GslWrappers.h"
#include "PhysicalConstants.h"
#include "PlasmaState.h"

#include <cmath>


namespace {

  const double me = PhysicalConstantsCGS::ElectronMass;
  const double hb = PhysicalConstantsCGS::Hbar;
  const double NePrefactor = pow(2.0*me, 1.5) / (2.0 * M_PI*M_PI * hb*hb*hb);

  class FermiDiracDistribution : public GSL::FunctionObject {
    private:
      const double xi;
      const PlasmaState p;
    public:
      FermiDiracDistribution(const double xi, const PlasmaState& p) : xi(xi), p(p) {}
      double operator()(const double x) const override {
        const double val = (1.0 + p.tau*x) * sqrt(x + p.tau * x*x/2.0) / (1.0 + exp(x-p.chi-xi));
        return val;
      }
  };
}



double Plasma::rhoFromNe(const double ne, const Composition& comp) {
  const double& mp = PhysicalConstantsCGS::ProtonMass;
  return ne * mp * comp.meanMolecularWeightPerElectron;
}



double Plasma::ne(const double chi, const double kt, const double tau) {
  const double i12 = gfdi(GFDI::Order12, chi, tau);
  const double i32 = gfdi(GFDI::Order32, chi, tau);
  return NePrefactor * pow(kt, 1.5) * (i12 + tau*i32);
}

double Plasma::ne(const double phi, const PlasmaState& p) {
  const double xi = fmax(0, phi/p.kt);
  return ne(p.chi + xi, p.kt, p.tau);
}

double Plasma::neBound(const double phi, const PlasmaState& p, const double cutoff) {
  const double xi = phi/p.kt;
  // conditions are:
  // * potential attractive ==> xi > 0
  // * integration bounds valid ==> (xi-cutoff) > 0
  if (xi > 0 and xi > cutoff) {
    FermiDiracDistribution fd(xi, p);
    const double eps = 1.e-6;
    return NePrefactor * pow(p.kt, 1.5) * GSL::integrate(fd, 0, xi-cutoff, eps, eps);
  } else {
    return 0;
  }
}

std::vector<double> Plasma::ni(const double phi, const PlasmaState& p) {
  const double xi = fmax(0, phi/p.kt);
  std::vector<double> ni = p.ni;
  for (size_t elem=0; elem<ni.size(); ++elem) {
    ni[elem] *= exp(-xi * p.comp.species[elem].element.Z);
  }
  return ni;
}



double Plasma::electronKineticEnergyDensity(const double phi, const PlasmaState& p) {
  const double xi = fmax(0, phi/p.kt);
  const double i32 = gfdi(GFDI::Order32, p.chi+xi, p.tau);
  const double i52 = gfdi(GFDI::Order52, p.chi+xi, p.tau);
  // TODO: check prefactor is the same here
  return NePrefactor * pow(p.kt, 2.5) * (i32 + p.tau*i52);
}

double Plasma::totalIonChargeDensity(const double phi, const PlasmaState& p) {
  const std::vector<double> ni = Plasma::ni(phi, p);
  double chargeDensity = 0;
  for (size_t elem=0; elem<ni.size(); ++elem) {
    chargeDensity += ni[elem] * p.comp.species[elem].element.Z;
  }
  return chargeDensity;
}



double Plasma::radiusWignerSeitz(const Element& e, const PlasmaState &p) {
  return pow((3*e.Z)/(4*M_PI*p.ne), 1.0/3.0);
}

double Plasma::energyWignerSeitz(const Element& e, const PlasmaState &p) {
  const double& qe = PhysicalConstantsCGS::ElectronCharge;
  return 0.9 * (qe*e.Z)*(qe*e.Z) / radiusWignerSeitz(e,p);
}

