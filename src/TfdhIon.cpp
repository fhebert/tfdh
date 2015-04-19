
#include "TfdhIon.h"

#include "Element.h"
#include "PlasmaFunctions.h"
#include "PlasmaState.h"
#include "PhysicalConstants.h"
#include "Species.h"
#include "TfdhOdeSolve.h"
#include "TfdhFunctions.h"
#include "TfdhSolution.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>


namespace {

  std::ostream& operator<<(std::ostream& s, const Species& sp) {
    s << sp.massFraction << " by mass of " << sp.element;
    return s;
  }

  template <typename T>
  std::ostream& operator<<(std::ostream& s, const std::vector<T>& v) {
    // do nothing if v is empty
    if (v.size()==1) {
      s << v[0];
    }
    else if (v.size() > 1) {
      s << "(";
      for (size_t i=0; i<v.size(); ++i) {
        s << v[i];
        if (i != v.size()-1)
          s << ", ";
      }
      s << ")";
    }
    return s;
  }

  std::ostream& operator<<(std::ostream& s, const Composition& c) {
    s << c.species;
    return s;
  }
}



TfdhIon::TfdhIon(const PlasmaState& plasmaState, const Element& element)
: ps(plasmaState),
  e(element),
  tfdh(TFDH::solve(e, ps)),
  numberBoundElectrons(TFDH::boundElectrons(tfdh, ps)),
  embeddingEnergies(TFDH::embeddingEnergy(tfdh, ps, e)),
  exclusionRadii(TFDH::exclusionRadii(tfdh, e, ps))
{}


void TfdhIon::printSummaryToFile(const std::string& filename) const
{
  std::ofstream outfile(filename);
  assert(outfile and "couldn't open file");
  //outfile.precision(16);

  outfile << "summary of TFDH ion-in-plasma calculation results\n";
  outfile << "code run on TODO\n"; // TODO
  outfile << "code revision = TODO\n"; // TODO
  outfile << "\n";

  outfile << "central ion = " << e << "\n";
  outfile << "\n";

  outfile << "plasma parameters\n";
  outfile << "composition = " << ps.comp << "\n";
  outfile << "rho = " << ps.rho << "\n";
  outfile << "t = " << ps.kt / PhysicalConstantsCGS::KBoltzmann << "\n";
  outfile << "ne = " << ps.ne << "\n";
  outfile << "ni = " << ps.ni << "\n";
  outfile << "tau = " << ps.tau << "\n";
  outfile << "chi = " << ps.chi << "\n";
  outfile << "\n";

  const double rws = Plasma::radiusWignerSeitz(e, ps);
  const double scale = PhysicalConstantsCGS::BohrRadius / e.Z;
  outfile << "Wigner-Seitz estimated quantities\n";
  outfile << "rws = " << rws << "\n";
  //outfile << "  r_ws / (a_0/Z_tr) = " << rws / scale << "\n";


  outfile << "rws = " << rws << ", rws*Ztr/a0 = " << rws*scale << "\n";
  for (double rex : exclusionRadii)
    outfile << "rex = " << rex << ", rex*Ztr/a0 = " << rex*scale << "\n";
  outfile << "\n";

  outfile << "TFDH global quantities\n";
  outfile << "number of bound electrons = " << numberBoundElectrons << "\n";
  outfile << "Z_net = " << e.Z - numberBoundElectrons << "\n";
  outfile << "embedding energy = " << embeddingEnergies.total/ps.kt << "\n";
  outfile << "\n";

  outfile << "embedding energy breakdown\n";
  outfile << "ion field energy:             " << embeddingEnergies.fi/ps.kt << "\n";
  outfile << "e- field energy:              " << embeddingEnergies.fe/ps.kt << "\n";
  outfile << "overcounting of field energy: " << embeddingEnergies.f2/ps.kt << "\n";
  outfile << "change in ion kinetic energy: " << embeddingEnergies.ki/ps.kt << "\n";
  outfile << "change in e- kinetic energy:  " << embeddingEnergies.ke/ps.kt << "\n";
  outfile << "energy from exchanging ions:  " << embeddingEnergies.ni/ps.kt << "\n";
  outfile << "energy from exchanging e-'s:  " << embeddingEnergies.ne/ps.kt << "\n";

  return;
}

