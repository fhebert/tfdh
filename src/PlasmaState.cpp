
#include "PlasmaState.h"
#include "Abundance.h"
#include "Composition.h"
#include "GslBrentWrapper.h"
#include "PhysicalConstants.h"
#include "PlasmaFunctions.h"

#include <cassert>
#include <limits>
#include <vector>


// helper functions used in initializing the plasma state
namespace {

  double computeNe(const double rho, const Composition& comp) {
    const double& mp = PhysicalConstantsCGS::ProtonMass;
    return rho / (mp * comp.meanMolecularWeightPerElectron);
  }

  std::vector<double> computeNi(double rho, const Composition& comp) {
    const std::vector<Abundance>& abundances = comp.abundances;
    const double& mp = PhysicalConstantsCGS::ProtonMass;
    std::vector<double> result(abundances.size());
    for (size_t elem=0; elem<abundances.size(); ++elem) {
      result[elem] = rho * abundances[elem].massFraction
                     / (mp * abundances[elem].element.A);
    }
    return result;
  }

  class ChiFunction : public GslFunction {
    const double ne_target_;
    const double kt_;
    const double tau_;
    public:
    ChiFunction(const double ne_t, const double kt, const double tau)
      : ne_target_(ne_t), kt_(kt), tau_(tau) {}
    double f(const double chi) const {
      // argument phi = 0 in this rootfind
      return ne_target_ - Plasma::ne(chi, 0.0, kt_, tau_);
    }
  };

  double invertForChi(const double ne, const double kt, const double tau) {
    // start by initializing some values
    const ChiFunction func(ne, kt, tau);

    // find a chi interval which brackets the root
    // ne(chi) is an INCREASING function, so ChiFunction(chi) is DECREASING
    double chiA = 0;
    double deltaA = func.f(chiA);
    if (deltaA == 0.0) return chiA;
    double chiB = (deltaA > 0) - (deltaA < 0); // this is the sign operator
    double deltaB = func.f(chiB);
    while (deltaA*deltaB > 0.0) {
      chiA = chiB;
      deltaA = deltaB;
      chiB = 2*chiB;
      deltaB = func.f(chiB);
    }

    // find the root
    const double chi_eps = std::numeric_limits<double>::epsilon();
    return GslBrent(func, chiA, chiB, chi_eps);
  }

} // helper namespace



PlasmaState::
PlasmaState(const double rho, const double kt, const Composition& comp,
    const bool isRel)
: rho(rho),
  kt(kt),
  comp(comp),
  isRelativistic(isRel),
  ne(computeNe(rho, comp)),
  ni(computeNi(rho, comp)),
  tau(isRel ? (kt / PhysicalConstantsCGS::MeCC) : 0.0),
  chi(invertForChi(ne, kt, tau))
{
  assert(kt>0);
  assert(rho>0);
}

