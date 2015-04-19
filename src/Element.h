
#ifndef TFDH_ELEMENT_H
#define TFDH_ELEMENT_H

#include <cassert>
#include <iostream>
#include <string>


class Element {
 public:
  Element(const unsigned a, const unsigned z, const std::string n)
    : A(a), Z(z), name(n)
  {
    assert(A>=Z and "Input has fewer nucleons than protons, is unphysical.");
    assert(A<300 and Z<120 and "Input lies outside the common periodic table.");
  }

  const unsigned A;
  const unsigned Z;
  const std::string name;
};


// Inline is primarily used to avoid "duplicate symbol" warnings arising from
// placing this operator definition in the header... I just don't want to make
// an extra .cpp file just to hold this two-liner.
inline std::ostream& operator<<(std::ostream& s, const Element& e) {
  s << e.name;
  return s;
}


#endif // TFDH_ELEMENT_H
