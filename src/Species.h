
#ifndef TFDH_SPECIES_H
#define TFDH_SPECIES_H

#include "Element.h"

#include <cassert>
#include <string>


struct Species {
  const double massFraction;
  const Element element;

  Species(const double m, const Element& e) : massFraction(m), element(e) {
    assert(massFraction > 0 and massFraction <= 1 and "Mass fraction must be in (0,1)!");
  }
};


#endif // TFDH_SPECIES_H
