
#include "TfdhIon.h"

#include "Element.h"
#include "IntegrateOverRadius.h"
#include "PhysicalConstants.h"
#include "PlasmaFunctions.h"
#include "PlasmaState.h"
#include "TfdhFunctions.h"
#include "TfdhOdeSolve.h"
#include "TfdhSolution.h"
#include "Utils.h"

#include <fstream>
#include <ostream>
#include <string>
#include <vector>


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


void TfdhIon::printRadialProfileToFile(const std::string& filename) const
{
  std::ofstream outfile(filename);
  assert(outfile and "couldn't open file");
  //outfile.precision(16);

  // print header:
  outfile << "# radial profile of TFDH ion-in-plasma calculation results\n";
  outfile << "# code run on TODO\n"; // TODO
  outfile << "# code revision = TODO\n"; // TODO
  outfile << "#\n";
  outfile << "# col #0 = radius\n";
  outfile << "# col #1 = potential [Zeff e^2 / r]\n";
  outfile << "# col #2 = electron density\n";
  outfile << "# col #3 = bound portion of electron density\n";
  outfile << "# col #4 = free portion of electron density (total - bound)\n";
  outfile << "# col #5 = total ion charge density (Zi * ni)\n";
  outfile << "# col #6 = enclosed net charge in units of q_e\n";
  outfile << "# col #7 = enclosed bound electron charge in units of q_e\n";

  // needed for cumulative distribution output
  const auto f_charge = [&] (const double r) -> double {
    return Plasma::totalIonChargeDensity(tfdh(r), ps) - Plasma::ne(tfdh(r), ps);
  };
  const auto f_nb = [&] (const double r) -> double {
    return Plasma::neBound(tfdh(r), ps);
  };
  double accumulate_charge = 0;
  double accumulate_numbound = 0;

  // print profiles:
  const std::string sep = "    ";
  for (size_t i=0; i<tfdh.r.size(); ++i) {
    // columns 0,1 -- the tfdh (r,phi) results
    outfile << tfdh.r[i] << sep << tfdh.phi[i];

    // columns 2,3,4 -- the electron densities
    const double ne = Plasma::ne(tfdh.phi[i], ps);
    const double neb = Plasma::neBound(tfdh.phi[i], ps);
    outfile << sep << ne << sep << neb << sep << ne-neb;

    // column 5 -- ion charge density
    outfile << sep << Plasma::totalIonChargeDensity(tfdh.phi[i], ps);

    // columns 6,7 -- the cumulative distributions
    if (i > 0) {
      accumulate_charge += integrateOverRadius(f_charge, tfdh.r[i-1], tfdh.r[i]);
      accumulate_numbound += integrateOverRadius(f_nb, tfdh.r[i-1], tfdh.r[i]);
    }
    outfile << sep << e.Z + accumulate_charge << sep << accumulate_numbound;
    outfile << "\n";
  }

  return;
}
