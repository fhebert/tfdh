
#ifndef TFDH_TFDH_FUNCTIONS_H
#define TFDH_TFDH_FUNCTIONS_H

#include <array>
#include <vector>

struct Element;
struct PlasmaState;
struct RadialFunction;


namespace TFDH {

  double boundElectrons(const RadialFunction& tfdh, const PlasmaState& p);
  std::vector<double> exclusionRadii(const RadialFunction& tfdh,
      const Element& e, const PlasmaState& p);

  double embeddingEnergy(const RadialFunction& tfdh, const PlasmaState& p,
      const Element& e);

}


#endif // TFDH_TFDH_FUNCTIONS_H
