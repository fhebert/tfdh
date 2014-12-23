

#include "Element.h"
#include "Composition.h"
#include "IntegrateTFDH.h"
#include "PlasmaState.h"
#include "PhysicalConstants.h"
#include "RadialFunction.h"

#include <string>
#include <iostream>



// helper namespace to hold some useful elements
namespace Elements {
  const Element H(1, 1, "Hydrogen");
  const Element He(4, 2, "Helium");
  const Element C(12, 6, "Carbon");
  const Element O(16, 8, "Oxygen");
  const Element Fe56(56, 26, "Iron-56");
} // namespace Elements



int main() {
  std::cout << ToString(Elements::O) << "\n";
  std::cout << ToString(Element(6,5,"weird isotope?")) << "\n";

  const double rho = 1e6;
  const double t = 1e7;
  const double kt = t * PhysicalConstantsCGS::KBoltzmann;

  std::vector<Abundance> ab;
  ab.push_back(Abundance(0.5, Elements::H));
  ab.push_back(Abundance(0.5, Elements::He));
  const Composition c(ab);
  const PlasmaState ps(rho, kt, c, false);

  std::cout << "plasma rho = " << ps.rho << "\n";
  std::cout << "plasma kt = " << ps.kt << "\n";
  std::cout << "plasma chi = " << ps.chi << "\n";
  std::cout << "plasma tau = " << ps.tau << "\n";
  std::cout << "plasma ne = " << ps.ne << "\n";
  for (double ni : ps.ni)
    std::cout << "plasma ni = " << ni << "\n";

  const RadialFunction tfdh = TFDH::solve(Elements::Fe56, ps);

  writeToFile(tfdh, "test_data.dat");

  return 0;
}
