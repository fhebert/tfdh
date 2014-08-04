
#include "PlasmaFunctions.h"

#include "Composition.h"
#include "Element.h"
#include "Gfdi.h"
#include "PhysicalConstants.h"
#include "PlasmaState.h"

#include <cmath>


namespace {
  const double& me = PhysicalConstantsCGS::ElectronMass;
  const double& hb = PhysicalConstantsCGS::Hbar;
  const double PrefactorNe = pow(2.0 * me, 1.5) / (2.0 * M_PI*M_PI * hb*hb*hb);
}


double Plasma::rhoFromNe(const double ne, const Composition& comp) {
  const double& mp = PhysicalConstantsCGS::ProtonMass;
  return ne * mp * comp.meanMolecularWeightPerElectron;
}


double Plasma::ne(const double chi, const double xi, const double kt,
    const double tau) {
  if (xi < 0) {
    return 0.0;
  }
  else {
    const double i12 = gfdi(GFDI::Order12, chi + xi, tau);
    const double i32 = gfdi(GFDI::Order32, chi + xi, tau);
    return PrefactorNe * pow(kt, 1.5) * (i12 + tau*i32);
  }
}

double Plasma::ne(const double phi, const PlasmaState& p) {
  return ne(p.chi, phi/p.kt, p.kt, p.tau);
}



std::vector<double> Plasma::ni(const double phi, const PlasmaState& p) {
  const double xi = phi/p.kt;
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

