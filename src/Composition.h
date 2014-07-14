
#ifndef TFDH_COMPOSITION_H
#define TFDH_COMPOSITION_H

#include "Abundance.h"
#include "Element.h"

#include <vector>


struct Composition {

  Composition(const Element& element); // simple c'tor for one-component plasmas
  Composition(const std::vector<Abundance>& abundances);

  const std::vector<Abundance> abundances;
  const double meanMolecularWeightPerElectron;
};


#endif // TFDH_COMPOSITION_H
