
#include "TfdhIon.h"

#include "Element.h"
#include "IntegrateOverRadius.h"
#include "PlasmaFunctions.h"
#include "PlasmaState.h"
#include "PhysicalConstants.h"
#include "TfdhOdeSolve.h"
#include "TfdhFunctions.h"
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
  outfile << "# col #4 = free portion of electron density (= 2-3)\n";
  outfile << "# col #5 = enclosed net charge in units of q_e\n";
  outfile << "# col #6 = enclosed bound electron charge in units of q_e\n";

  // intermediate quantities
  // TODO: see about cleaning this up
  std::vector<double> ne(tfdh.r.size(), 0);
  std::vector<double> ntot(tfdh.r.size(), 0);
  for (size_t i=0; i<ne.size(); ++i) {
    ne[i] = Plasma::ne(tfdh.phi[i], ps);
    ntot[i] = Plasma::totalIonChargeDensity(tfdh.phi[i], ps) - ne[i];
  }
  const std::vector<double> neb = TFDH::boundElectronDensity(tfdh, ps);

  const std::vector<double> enc_plasma_charge = accumulateOverRadius(tfdh.r, ntot);
  const std::vector<double> enc_bound_elec = accumulateOverRadius(tfdh.r, neb);

  // print profiles:
  const std::string sep = "    ";
  for (size_t i=0; i<tfdh.r.size(); ++i) {
    outfile << tfdh.r[i] << sep << tfdh.phi[i] << sep;
    outfile << ne[i] << sep << neb[i] << sep << ne[i]-neb[i] << sep;
    outfile << e.Z + enc_plasma_charge[i] << sep << enc_bound_elec[i];
    outfile << "\n";
  }

  return;
}
