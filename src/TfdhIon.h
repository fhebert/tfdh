
#ifndef TFDH_TFDH_ION_H
#define TFDH_TFDH_ION_H

#include "Element.h"
#include "PlasmaState.h"
#include "TfdhFunctions.h"
#include "TfdhOdeSolve.h"
#include "TfdhSolution.h"

#include <vector>


class TfdhIon {
  public:
    TfdhIon(const PlasmaState& plasmaState, const Element& element);

    void printSummaryToFile(const std::string& filename, const std::string& time="<no time given>") const;
    void printRadialProfileToFile(const std::string& filename, const std::string& time="<no time given>") const;

    const PlasmaState ps;
    const Element e;
    const TfdhSolution tfdh;
    const double numberBoundElectrons;
    const TFDH::EnergyDeltas embeddingEnergies;
    const std::vector<double> exclusionRadii;
};


#endif // TFDH_TFDH_ION_H
