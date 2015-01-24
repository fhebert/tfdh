
#ifndef TFDH_PLASMA_FUNCTIONS_H
#define TFDH_PLASMA_FUNCTIONS_H

#include <vector>

struct Composition;
struct Element;
struct PlasmaState;


namespace Plasma {

  double rhoFromNe(double ne, const Composition& comp);

  double ne(double chi, double kt, double tau);
  double ne(double phi, const PlasmaState& p);

  double neBound(double phi, const PlasmaState& p, double cutoff=0);

  // electron kinetic energy density
  double neKinetic(double phi, const PlasmaState& p);

  std::vector<double> ni(double phi, const PlasmaState& p);
  double totalIonChargeDensity(double phi, const PlasmaState& p);

  double radiusWignerSeitz(const Element& e, const PlasmaState &p);
  double energyWignerSeitz(const Element& e, const PlasmaState &p);

}


#endif // TFDH_PLASMA_FUNCTIONS_H
