
#include "PlasmaFunctions.h"

#include "Composition.h"
#include "Element.h"
#include "Gfdi.h"
#include "GslWrappers.h"
#include "PhysicalConstants.h"
#include "PlasmaState.h"

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
}



double Plasma::rhoFromNe(const double ne, const Composition& comp) {
  const double& mp = PhysicalConstantsCGS::ProtonMass;
  return ne * mp * comp.meanMolecularWeightPerElectron;
}


double Plasma::ne(const double chi, const double kt, const double tau) {
  const double i12 = gfdi(GFDI::Order12, chi, tau);
  const double i32 = gfdi(GFDI::Order32, chi, tau);
  const double& NePrefactor = PhysicalConstantsCGS::NePrefactor;
  return NePrefactor * pow(kt, 1.5) * (i12 + tau*i32);
}

double Plasma::ne(const double phi, const PlasmaState& p) {
  const double xi = fmax(0, phi/p.kt);
  return ne(p.chi + xi, p.kt, p.tau);
}


double Plasma::neBound(const double phi, const PlasmaState& p) {
  const double xi = phi/p.kt;
  if (xi < 0) {
    return 0;
  } else {
    const double& NePrefactor = PhysicalConstantsCGS::NePrefactor;
    FermiDiracDistribution fd(xi, p);
    const double eps = 1.e-6;
    return NePrefactor * pow(p.kt, 1.5) * GSL::integrate(fd, 0, xi, eps, eps);
  }
}


double Plasma::neKinetic(const double phi, const PlasmaState& p) {
  const double xi = fmax(0, phi/p.kt);
  const double i32 = gfdi(GFDI::Order32, p.chi+xi, p.tau);
  const double i52 = gfdi(GFDI::Order52, p.chi+xi, p.tau);
  const double& NePrefactor = PhysicalConstantsCGS::NePrefactor;
  return NePrefactor * pow(p.kt, 2.5) * (i32 + p.tau*i52);
}


std::vector<double> Plasma::ni(const double phi, const PlasmaState& p) {
  const double xi = (phi >= 0.0) ? phi/p.kt : 0.0;
  std::vector<double> nis(p.ni);
  for (size_t elem=0; elem<nis.size(); ++elem) {
    nis[elem] *= exp(-xi * p.comp.abundances[elem].element.Z);
  }
  return nis;
}

double Plasma::totalIonChargeDensity(const double phi, const PlasmaState& p) {
  const std::vector<double> nis = ni(phi, p);
  double chargeDensity = 0;
  for (size_t elem=0; elem<nis.size(); ++elem) {
    chargeDensity += nis[elem] * p.comp.abundances[elem].element.Z;
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

