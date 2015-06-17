
#ifndef TFDH_TFDH_FUNCTIONS_H
#define TFDH_TFDH_FUNCTIONS_H

#include <vector>

class Element;
class PlasmaState;
class TfdhSolution;


namespace TFDH {

  struct EnergyDeltas {
    const double fi;
    const double fe;
    const double f2;
    const double ki;
    const double ke;
    const double ni;
    const double ne;
    const double total;
  };

  std::vector<double> boundElectronDensity(const TfdhSolution& tfdh, const PlasmaState& p, double cutoff=0);
  double boundElectrons(const TfdhSolution& tfdh, const PlasmaState& p, double cutoff=0);
  std::vector<double> exclusionRadii(const TfdhSolution& tfdh, const Element& e, const PlasmaState& p);

  EnergyDeltas embeddingEnergy(const TfdhSolution& tfdh, const PlasmaState& p, const Element& e);

}


#endif // TFDH_TFDH_FUNCTIONS_H
