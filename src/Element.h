
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


namespace Elements {
  const Element H(1, 1, "Hydrogen");
  const Element He(4, 2, "Helium");
  const Element C(12, 6, "Carbon");
  const Element O(16, 8, "Oxygen");
  const Element Fe56(56, 26, "Iron-56");
} // namespace Elements


#endif // TFDH_ELEMENTS_H
