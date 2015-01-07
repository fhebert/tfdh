

#include "Element.h"
#include "Composition.h"
#include "PlasmaState.h"
#include "PhysicalConstants.h"
#include "RadialFunction.h"
#include "TfdhOdeSolve.h"

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

  const double rho = 1e8;
  const double t = 1e5;
  const double kt = t * PhysicalConstantsCGS::KBoltzmann;

  std::vector<Abundance> ab;
  ab.push_back(Abundance(1.0, Elements::He));
  const Composition c(ab);
  const PlasmaState ps(rho, kt, c, false);

  std::cout << "plasma rho = " << ps.rho << "\n";
  std::cout << "plasma kt = " << ps.kt << "\n";
  std::cout << "plasma chi = " << ps.chi << "\n";
  std::cout << "plasma tau = " << ps.tau << "\n";
  std::cout << "plasma ne = " << ps.ne << "\n";
  for (double ni : ps.ni)
    std::cout << "plasma ni = " << ni << "\n";

  //const RadialFunction tfdh = TFDH::solve(Elements::Fe56, ps);
  //writeToFile(tfdh, "test_data.dat");

  const double bound_electrons = TFDH::boundElectrons(Elements::Fe56, ps);
  std::cout << "number of bound electrons = " << bound_electrons << "\n";
  std::cout << "=> effective Z_net = " << Elements::Fe56.Z - bound_electrons << "\n";

  return 0;
}
