
#ifndef TFDH_TFDH_FUNCTIONS_H
#define TFDH_TFDH_FUNCTIONS_H

#include <vector>

struct Element;
struct PlasmaState;
struct RadialFunction;


namespace TFDH {

  double neBound(double phi, const PlasmaState& p);
  RadialFunction neBound(const RadialFunction& tfdh, const PlasmaState& p);

  double boundElectrons(const RadialFunction& tfdh, const PlasmaState& p);

  std::vector<double> exclusionRadii(const RadialFunction& tfdh,
      const Element& e, const PlasmaState& p);
}


#endif // TFDH_TFDH_FUNCTIONS_H
