
#ifndef TFDH_COMPOSITION_H
#define TFDH_COMPOSITION_H

#include <vector>

struct Element;
struct Species;


struct Composition {
  const std::vector<Species> species;
  const double meanMolecularWeightPerElectron;

  Composition(const Element& element); // simple c'tor for one-component plasmas
  Composition(const std::vector<Species>& species);
};


#endif // TFDH_COMPOSITION_H
