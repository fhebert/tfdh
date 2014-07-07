

#include "Element.h"
#include "Composition.h"
#include "PlasmaState.h"

#include <string>
#include <iostream>

int main() {
  std::cout << ToString(Elements::O) << "\n";
  std::cout << ToString(Element(6,5,"weird isotope?")) << "\n";

  const double rho = 1e6;
  const double t = 1e7;

  std::vector<Abundance> ab;
  ab.push_back(Abundance(0.5, Elements::H));
  ab.push_back(Abundance(0.5, Elements::He));
  const Composition c(ab);
  const PlasmaState ps(t, rho, c, false);

  std::cout << "plasma rho = " << ps.rho << "\n";
  std::cout << "plasma kt = " << ps.kt << "\n";
  std::cout << "plasma chi = " << ps.chi << "\n";
  std::cout << "plasma tau = " << ps.tau << "\n";
  std::cout << "plasma ne = " << ps.ne << "\n";
  for (double ni : ps.ni)
    std::cout << "plasma ni = " << ni << "\n";

  return 0;
}
