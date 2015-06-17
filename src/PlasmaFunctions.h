
#ifndef TFDH_PLASMA_FUNCTIONS_H
#define TFDH_PLASMA_FUNCTIONS_H

#include <vector>

class Composition;
class Element;
class PlasmaState;


namespace Plasma {

  // for initializing a PlasmaState from an ne instead of a rho
  double rhoFromNe(double ne, const Composition& comp);

  // number densities
  double ne(double chi, double kt, double tau);
  double ne(double phi, const PlasmaState& p);
  double neBound(double phi, const PlasmaState& p, double cutoff=0);
  std::vector<double> ni(double phi, const PlasmaState& p);

  // energy/charge densities
  double electronKineticEnergyDensity(double phi, const PlasmaState& p);
  double totalIonChargeDensity(double phi, const PlasmaState& p);

  double radiusWignerSeitz(const Element& e, const PlasmaState &p);
  double energyWignerSeitz(const Element& e, const PlasmaState &p);

}


#endif // TFDH_PLASMA_FUNCTIONS_H
