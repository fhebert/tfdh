
#ifndef TFDH_ELEMENT_H
#define TFDH_ELEMENT_H

#include <cassert>
#include <string>


struct Element {
  const unsigned A;
  const unsigned Z;
  const std::string name;

  Element(unsigned a, unsigned z, const std::string n) : A(a), Z(z), name(n) {
    assert(A>=Z and "Elements have at least as many nucleons as protons.");
    assert(A<200 and Z<200 and "Input lies way outside the periodic table.");
  }
};


#endif // TFDH_ELEMENT_H
