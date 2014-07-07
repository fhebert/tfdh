
#ifndef TFDH_ABUNDANCE_H
#define TFDH_ABUNDANCE_H

#include "Element.h"

#include <cassert>
#include <string>


struct Abundance {
  const double massFraction;
  const Element element;

  Abundance(const double m, const Element& e) : massFraction(m), element(e) {
    assert(massFraction > 0 && massFraction <= 1 &&
        "Mass fraction must be in (0,1)!");
  }
};


#endif // TFDH_ABUNDANCE_H
