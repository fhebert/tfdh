
#ifndef TFDH_TFDH_FUNCTIONS_H
#define TFDH_TFDH_FUNCTIONS_H

#include <vector>

struct Element;
struct PlasmaState;
struct TfdhSolution;


namespace TFDH {

  std::vector<double> boundElectronDensity(const TfdhSolution& tfdh, const PlasmaState& p, double cutoff=0);
  double boundElectrons(const TfdhSolution& tfdh, const PlasmaState& p, double cutoff=0);
  std::vector<double> exclusionRadii(const TfdhSolution& tfdh, const Element& e, const PlasmaState& p);

  double embeddingEnergy(const TfdhSolution& tfdh, const PlasmaState& p, const Element& e);

}


#endif // TFDH_TFDH_FUNCTIONS_H
