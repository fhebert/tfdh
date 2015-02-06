
#include "PlasmaState.h"

#include "Composition.h"
#include "GslWrappers.h"
#include "PhysicalConstants.h"
#include "PlasmaFunctions.h"
#include "Species.h"

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
    const double& mp = PhysicalConstantsCGS::ProtonMass;
    std::vector<double> result(comp.species.size());
    for (size_t elem=0; elem<result.size(); ++elem) {
      result[elem] = rho * comp.species[elem].massFraction
                     / (mp * comp.species[elem].element.A);
    }
    return result;
  }

  class NeErrorFromChi : public GSL::FunctionObject {
    private:
      const double ne0, kt, tau;
    public:
      NeErrorFromChi(const double ne, const double kt, const double tau)
        : ne0(ne), kt(kt), tau(tau) {}
      double operator()(const double chi) const {
        // argument phi = 0
        return Plasma::ne(chi, kt, tau) - ne0;
      }
  };

  double invertForChi(const double ne, const double kt, const double tau) {
    // function to find root of:
    const NeErrorFromChi deltaNe(ne, kt, tau);

    // find a chi interval which brackets the root
    double chiA = 0;
    double deltaA = deltaNe(chiA);
    if (deltaA == 0.0) return chiA;
    // ne(chi) is an INCREASING function, so deltaNe(chi) is too
    const double signA = (deltaA > 0) - (deltaA < 0);
    double chiB = - signA;
    double deltaB = deltaNe(chiB);
    while (deltaA*deltaB > 0.0) {
      chiA = chiB;
      deltaA = deltaB;
      chiB = 2*chiB;
      deltaB = deltaNe(chiB);
    }

    // find the root -- for this shooting problem we need max accuracy
    const double eps = std::numeric_limits<double>::epsilon();
    return GSL::findRoot(deltaNe, chiA, chiB, eps, eps);
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

