
#ifndef TFDH_COMPOSITION_H
#define TFDH_COMPOSITION_H

#include "Element.h"

#include <ostream>
#include <vector>


// a "helper" struct to pair information used in Composition
struct Species {
  const double massFraction;
  const Element element;
};


struct Composition {
  const std::vector<Species> species;
  const double meanMolecularWeightPerElectron;

  Composition(const Element& element); // simple c'tor for one-component plasmas
  Composition(const std::vector<Species>& species);
};


std::ostream& operator<<(std::ostream& s, const Species& sp);
std::ostream& operator<<(std::ostream& s, const Composition& c);


#endif // TFDH_COMPOSITION_H
