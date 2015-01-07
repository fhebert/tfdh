
#include "IntegrateTFDH.h"

#include "Element.h"
#include "IntegrateOverRadius.h"
#include "PlasmaFunctions.h"
#include "PlasmaState.h"
#include "PhysicalConstants.h"
#include "RadialFunction.h"
#include "TfdhFunctions.h"

#include <cassert>
#include <cmath>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>


namespace {

  struct TfdhParams {
    const Element& e;
    const PlasmaState& p;
    TfdhParams(const Element& re, const PlasmaState& rp) : e(re), p(rp) {}
  };

  int tfdhOde(double r, const double f[], double dfdr[], void *params) {
    const double& qe = PhysicalConstantsCGS::ElectronCharge;
    const double phi = qe*f[0]/r;
    const PlasmaState& p = static_cast<TfdhParams*>(params)->p;
    const double ne = Plasma::ne(phi, p);
    const double ionChargeDensity = Plasma::totalIonChargeDensity(phi, p);
    dfdr[0] = f[1];
    dfdr[1] = -4.0*M_PI*qe * r * (ionChargeDensity - ne);
    return GSL_SUCCESS;
  }
}


RadialFunction TFDH::solve(const Element& e, const PlasmaState& p)
{
  // set parameters here:
  const double rws = Plasma::radiusWignerSeitz(e,p);
  const double ri = 1e-5 * rws;
  // TODO: large outer radius set to avoid assert failures when rootfinding
  //       for the central potential. this means we're no longer checking that
  //       rf < rws
  //       should add this check back in somewhere
  const double rf = 1e3 * rws;

  const double dv0 = findPotentialRoot(e, p, ri, rf);
  return integrateODE(e, p, ri, rf, dv0);
}


RadialFunction TFDH::integrateODE(const Element& e, const PlasmaState& p,
    const double r_init, const double r_final, const double dv0)
{
  TfdhParams params(e,p);
  const size_t dim = 2;
  const double eps_abs = 1e-6;
  const double eps_rel = 0;
  double r = r_init;
  double dr = r_init;
  const double& qe = PhysicalConstantsCGS::ElectronCharge;
  double solution[dim] = {qe*e.Z + r_init*dv0, dv0};

  gsl_odeiv2_step* step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd, dim);
  gsl_odeiv2_control* ctrl = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
  gsl_odeiv2_evolve* ev = gsl_odeiv2_evolve_alloc(dim);
  gsl_odeiv2_system sys = {tfdhOde, nullptr, dim, &params};

  // TODO: fix units. too many quanities: solution / phi / xi / ...
  std::vector<double> radii, potentials;
  radii.push_back(r_init);
  potentials.push_back(qe*solution[0]/r_init);

  while (r < r_final) {
    const int status = gsl_odeiv2_evolve_apply(ev, ctrl, step, &sys, &r, r_final, &dr, solution);
    assert(status==GSL_SUCCESS);
    radii.push_back(r);
    potentials.push_back(qe*solution[0]/r);
    if (solution[0] <= 0 or (solution[1]-solution[0])/r > 0) break;
  }
  assert(r < r_final and "integrated ODE until final radius without terminating!");

  gsl_odeiv2_evolve_free(ev);
  gsl_odeiv2_control_free(ctrl);
  gsl_odeiv2_step_free(step);
  return RadialFunction(radii, potentials);
}


double TFDH::findPotentialRoot(const Element& e, const PlasmaState& p,
    const double r_init, const double r_final)
{
  // find interval that brackets correct potential
  double v_low = 0;
  double v_high = 0;
  {
    bool success = false;
    const double v_step = 100.0;
    const int bracket_attempts = 100;
    for (int i=0; i<bracket_attempts; ++i) {
      const auto& rf = integrateODE(e, p, r_init, r_final, v_low);
      if (rf.data.back() < 0.0) {
        success = true;
        break;
      }
      v_high = v_low;
      v_low -= v_step;
    }
    assert(success and "failed to bracket potential root");
  }

  // find "root" within this bracket
  // we aren't looking for a traditional root, just the boundary between
  // potentials that diverge to +infty or -infty... so do something simple
  double v_mid = 0.0;
  {
    bool success = false;
    const int root_attempts = 100;
    for (int i=0; i<root_attempts; ++i) {
      v_mid = (v_low + v_high)/2.0;
      if (v_mid==v_low or v_mid==v_high) { // underflow
        success = true;
        break;
      }

      const auto& rf = integrateODE(e, p, r_init, r_final, v_mid);
      ((rf.data.back() >= 0.0) ? v_high : v_low) = v_mid;
    }
    assert(success and "failed to find potential root within bracket");
  }

  // TODO: return the corresponding RadialFunction, to avoid computing it
  //       a second time!
  return v_mid;
}


double TFDH::boundElectrons(const Element& e, const PlasmaState& p)
{
  const RadialFunction& tfdh = solve(e, p);
  const RadialFunction& neBound = TFDH::neBound(tfdh, p);
  return integrateOverRadius(neBound);
}
