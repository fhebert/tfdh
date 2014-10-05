
#include "Abundance.h"
#include "Composition.h"
#include "Element.h"

#include <cassert>
#include <vector>


namespace {

  double compute_mu_e(const std::vector<Abundance>& abundances) {
    double inv_mu_e = 0;
    for (const Abundance& a : abundances) {
      // NOTE: This math assumes a totally ionized plasma, otherwise need to
      // multiply each term in the sum by the ionization fraction.
      inv_mu_e += a.massFraction * a.element.Z / (double) a.element.A;
    }
    return 1.0/inv_mu_e;
  }

} // helper namespace


Composition::Composition(const Element& element)
: abundances(std::vector<Abundance>(1, Abundance(1.0, element))),
  meanMolecularWeightPerElectron(element.Z / (double) element.A)
{}

Composition::Composition(const std::vector<Abundance>& abundances)
: abundances(abundances),
  meanMolecularWeightPerElectron(compute_mu_e(abundances))
{
  // sanity check on the mass fractions: they must add to 1
  double totalMassFraction = 0;
  for (const Abundance& a : abundances) {
    totalMassFraction += a.massFraction;
  }
  assert(totalMassFraction == 1.0);
}
