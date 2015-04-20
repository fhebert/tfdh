
#include "Composition.h"

#include "Element.h"
#include "Utils.h"

#include <cassert>
#include <cmath>
#include <ostream>
#include <vector>


namespace {
  double compute_mu_e(const std::vector<Species>& species) {
    double inv_mu_e = 0;
    for (const Species& s : species) {
      // NOTE: This math assumes a totally ionized plasma, otherwise need to
      // multiply each term in the sum by the ionization fraction.
      inv_mu_e += s.massFraction * s.element.Z / static_cast<double>(s.element.A);
    }
    return 1.0/inv_mu_e;
  }
} // helper namespace



Composition::Composition(const Element& element)
: species({{1.0, element}}),
  meanMolecularWeightPerElectron(element.Z / static_cast<double>(element.A))
{}


Composition::Composition(const std::vector<Species>& species)
: species(species),
  meanMolecularWeightPerElectron(compute_mu_e(species))
{
  // sanity checks on the mass fractions:
  // each must be in (0,1) and they must sum to 1
  double totalMassFraction = 0;
  for (const Species& s : species) {
    assert(s.massFraction > 0 and s.massFraction <= 1);
    totalMassFraction += s.massFraction;
  }
  assert(fabs(totalMassFraction-1.0) < 1e-15);
}



std::ostream& operator<<(std::ostream& s, const Species& sp) {
  s << sp.massFraction << " by mass of " << sp.element;
  return s;
}


std::ostream& operator<<(std::ostream& s, const Composition& c) {
  s << c.species;
  return s;
}
