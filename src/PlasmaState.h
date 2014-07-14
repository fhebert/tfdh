
#ifndef TFDH_PLASMA_STATE_H
#define TFDH_PLASMA_STATE_H

#include "Composition.h"

#include <vector>


struct PlasmaState {

  PlasmaState(double rho, double kt, const Composition& comp, bool isRel);

  // these "primary" variables are sufficient to define the state uniquely
  const double rho;
  const double kt;
  const Composition comp;
  const bool isRelativistic;

  // these "secondary" variables are computed from the primary ones
  const double ne;
  const std::vector<double> ni;
  const double tau;
  const double chi;
};


#endif // TFDH_PLASMA_STATE_H
