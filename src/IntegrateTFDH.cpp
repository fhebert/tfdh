
#include "IntegrateTFDH.h"

#include "Element.h"
#include "PlasmaFunctions.h"
#include "PlasmaState.h"
#include "PhysicalConstants.h"
#include "RadialFunction.h"

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
    return 0;
  }
}


RadialFunction TFDH::solve(const Element& e, const PlasmaState& p)
{
  // set parameters here:
  const double ri = 1e-10;
  const double rf = 1e-6;

  const double dv0 = rootfindPotential(e, p, ri, rf);
  return integrateODE(e, p, ri, rf, dv0);
}


RadialFunction TFDH::integrateODE(const Element& e, const PlasmaState& p,
    const double r_init, const double r_final, const double dv0)
{
  const double dr_start = r_init;
  const double eps_abs = 1e-6;
  const double eps_rel = 0;
  TfdhParams params(e,p);
  gsl_odeiv2_system sys = {tfdhOde, nullptr, 2, &params};
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys,
      gsl_odeiv2_step_rk8pd, dr_start, eps_abs, eps_rel);

  // TODO -- clean all this up
  const double& qe = PhysicalConstantsCGS::ElectronCharge;
  double r = r_init;
  double solution[2] = {qe*e.Z + r_init*dv0, dv0};

  std::vector<double> radii, potentials;

  const int numsteps = 100;
  for (int i=0; i<numsteps; ++i) {
    const double ri = r_init + (i+1)*(r_final-r_init)/numsteps;
    const int status = gsl_odeiv2_driver_apply(d, &r, ri, solution);
    assert(status==GSL_SUCCESS);
    radii.push_back(ri);
    potentials.push_back(solution[0]/ri);
  }

  gsl_odeiv2_driver_free(d);
  return RadialFunction(radii, potentials);
}


double TFDH::rootfindPotential(const Element& e, const PlasmaState& p,
    const double r_init, const double r_final)
{
  // find interval that brackets correct potential
  double v0 = 0;
  double v1 = 0;
  {
    bool success = false;
    const double dv = 1.0;
    const int bracket_attempts = 100;
    for (int i=0; i<bracket_attempts; ++i) {
      const auto& rf = integrateODE(e, p, r_init, r_final, v0);
      if (rf.data.back() < 0.0) {
        success = true;
        break;
      }
      v1 = v0;
      v0 -= dv;
    }
    assert(success and "failed to bracket potential root");
  }

  // iterate guess
  // avoid fancy rootfinders here, because we aren't dealing with a smooth
  // function and its root -- rather, we are trying to find the boundary
  // between two distinct sets: dv0 too large, dv0 too small
  double v_try = 0.0;
  {
    bool success = false;
    const int root_attempts = 100;
    for (int i=0; i<root_attempts; ++i) {
      v_try = (v0+v1)/2.0;
      if (v_try==v0 or v_try==v1) { // underflow
        success = true;
        break;
      }

      const auto& rf = integrateODE(e, p, r_init, r_final, v_try);
      if (rf.data.back() >= 0.0) {
        v1 = v_try;
      } else {
        v0 = v_try;
      }
    }
    assert(success and "failed to find potential root within bracket");
  }

  return v_try;
}


