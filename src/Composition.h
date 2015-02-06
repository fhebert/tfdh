
#ifndef TFDH_COMPOSITION_H
#define TFDH_COMPOSITION_H

#include "Element.h"
#include "Species.h"

#include <vector>


struct Composition {

  Composition(const Element& element); // simple c'tor for one-component plasmas
  Composition(const std::vector<Species>& species);

  const std::vector<Species> species;
  const double meanMolecularWeightPerElectron;
};


#endif // TFDH_COMPOSITION_H
