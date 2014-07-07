
#ifndef TFDH_COMPOSITION_H
#define TFDH_COMPOSITION_H

#include "Abundance.h"

#include <vector>


class Composition {

  public:
  Composition(const std::vector<Abundance>& abundances);

  const std::vector<Abundance> abundances;
  const double meanMolecularWeightPerElectron;
};


#endif // TFDH_COMPOSITION_H
