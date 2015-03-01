

#include "Element.h"
#include "Composition.h"
#include "IntegrateOverRadius.h"
#include "PhysicalConstants.h"
#include "PlasmaFunctions.h"
#include "PlasmaState.h"
#include "Species.h"
#include "TfdhIon.h"
#include "TfdhFunctions.h"
#include "TfdhOdeSolve.h"
#include "TfdhSolution.h"

#include <iostream>
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



int main() {

  const double rho = 1e8;
  const double t = 1e5;
  const double kt = t * PhysicalConstantsCGS::KBoltzmann;

  std::vector<Species> species;
  species.push_back(Species(1.0, Elements::He));
  const Composition c(species);
  const PlasmaState ps(rho, kt, c, false);

  std::cout << "plasma rho = " << ps.rho << "\n";
  std::cout << "plasma kt = " << ps.kt << "\n";
  std::cout << "plasma chi = " << ps.chi << "\n";
  std::cout << "plasma tau = " << ps.tau << "\n";
  std::cout << "plasma ne = " << ps.ne << "\n";
  for (double ni : ps.ni)
    std::cout << "plasma ni = " << ni << "\n";

  const TfdhIon ion(ps, Elements::Fe56);

  std::cout << "number of bound electrons = " << ion.numberBoundElectrons << "\n";
  std::cout << "=> effective Z_net = " << Elements::Fe56.Z - ion.numberBoundElectrons << "\n";

  const double rws = Plasma::radiusWignerSeitz(Elements::Fe56, ps);
  const double scale = 26/PhysicalConstantsCGS::BohrRadius;

  std::cout << "rws = " << rws << ", rws*Ztr/a0 = " << rws*scale << "\n";
  for (double rex : ion.exclusionRadii)
    std::cout << "rex = " << rex << ", rex*Ztr/a0 = " << rex*scale << "\n";

  std::cout << "\nembedding energies in units of kT:\n";
  std::cout << "ion field energy:             " << ion.embeddingEnergies.fi/ps.kt << "\n";
  std::cout << "e- field energy:              " << ion.embeddingEnergies.fe/ps.kt << "\n";
  std::cout << "overcounting of field energy: " << ion.embeddingEnergies.f2/ps.kt << "\n";
  std::cout << "change in ion kinetic energy: " << ion.embeddingEnergies.ki/ps.kt << "\n";
  std::cout << "change in e- kinetic energy:  " << ion.embeddingEnergies.ke/ps.kt << "\n";
  std::cout << "energy from exchanging ions:  " << ion.embeddingEnergies.ni/ps.kt << "\n";
  std::cout << "energy from exchanging e-'s:  " << ion.embeddingEnergies.ne/ps.kt << "\n";
  std::cout << "sum of everything        : " << ion.embeddingEnergies.total/ps.kt << "\n";
  //std::cout << "sum of the ones i like.. : " << (fi+fe+f2+ki+ke)/p.kt << "\n";

  return 0;
}
