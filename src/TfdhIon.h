
#ifndef TFDH_TFDH_ION_H
#define TFDH_TFDH_ION_H

#include "Element.h"
#include "PlasmaState.h"
#include "TfdhOdeSolve.h"
#include "TfdhFunctions.h"
#include "TfdhSolution.h"

#include <vector>


struct TfdhIon {

  const PlasmaState ps;
  const Element e;
  const TfdhSolution tfdh;
  const double numberBoundElectrons;
  const TFDH::EnergyDeltas embeddingEnergies;
  const std::vector<double> exclusionRadii;

  TfdhIon(const PlasmaState& plasmaState, const Element& element)
  : ps(plasmaState), e(element),
    tfdh(TFDH::solve(e, ps)),
    numberBoundElectrons(TFDH::boundElectrons(tfdh, ps)),
    embeddingEnergies(TFDH::embeddingEnergy(tfdh, ps, e)),
    exclusionRadii(TFDH::exclusionRadii(tfdh, e, ps))
  {}

};


#endif // TFDH_TFDH_ION_H
