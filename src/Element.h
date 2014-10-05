
#ifndef TFDH_ELEMENTS_H
#define TFDH_ELEMENTS_H

#include <cassert>
#include <string>


struct Element {
  const unsigned A;
  const unsigned Z;
  const std::string name;

  Element(unsigned a, unsigned z, const std::string n) : A(a), Z(z), name(n) {
    assert(A>=Z && "Elements have at least as many nucleons as protons!");
  }
};


inline std::string ToString(const Element& e) {
  return "element: " + e.name + " with A=" + std::to_string(e.A)
    + ", Z=" + std::to_string(e.Z);
}


#endif // TFDH_ELEMENTS_H
