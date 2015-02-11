
#ifndef TFDH_ELEMENTS_H
#define TFDH_ELEMENTS_H

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


inline std::string toString(const Element& e) {
  return "element: " + e.name + " with A=" + std::to_string(e.A)
    + ", Z=" + std::to_string(e.Z);
}


#endif // TFDH_ELEMENTS_H
