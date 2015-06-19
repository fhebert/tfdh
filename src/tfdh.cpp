

#include "Element.h"
#include "Composition.h"
#include "IntegrateOverRadius.h"
#include "PhysicalConstants.h"
#include "PlasmaFunctions.h"
#include "PlasmaState.h"
#include "TfdhIon.h"
#include "TfdhFunctions.h"
#include "TfdhOdeSolve.h"
#include "TfdhSolution.h"

#include <cassert>
#include <ctime>
#include <ostream>
#include <string>
#include <vector>



// helper namespace to hold some useful elements
namespace Elements {
  const Element H(1, 1, "Hydrogen");
  const Element He(4, 2, "Helium");
  const Element C(12, 6, "Carbon");
  const Element O(16, 8, "Oxygen");
  const Element Fe56(56, 26, "Iron-56");
} // namespace Elements


// helper namespace for functions
namespace {
  std::string getTime() {
    std::time_t t = std::time(nullptr);
    char str[128];
    assert(std::strftime(str, sizeof(str), "%c", std::localtime(&t)));
    return std::string(str);
  }
} // anon namespace



int main() {

  std::string time = getTime();

  const double rho = 1e3;
  const double t = 1e8;
  const double kt = t * PhysicalConstantsCGS::KBoltzmann;

  const auto species = std::vector<Species> {{ {0.7, Elements::He}, {0.3, Elements::H} }};
  const PlasmaState ps(rho, kt, Composition(species), false);

  const TfdhIon ion(ps, Elements::Fe56);
  ion.printSummaryToFile("summary.data", time);
  ion.printRadialProfileToFile("profile.data", time);

  return 0;
}
