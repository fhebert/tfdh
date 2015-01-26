
#ifndef TFDH_TFDH_FUNCTIONS_H
#define TFDH_TFDH_FUNCTIONS_H

#include <array>
#include <vector>

struct Element;
struct PlasmaState;
struct RadialFunction;


namespace TFDH {

  RadialFunction boundElectronDensity(const RadialFunction& tfdh, const PlasmaState& p, double cutoff=0);
  double boundElectrons(const RadialFunction& tfdh, const PlasmaState& p, double cutoff=0);
  std::vector<double> exclusionRadii(const RadialFunction& tfdh,
      const Element& e, const PlasmaState& p);

  double embeddingEnergy(const RadialFunction& tfdh, const PlasmaState& p,
      const Element& e);

}


#endif // TFDH_TFDH_FUNCTIONS_H
