
#include "TfdhOdeSolve.h"

#include "Element.h"
#include "PhysicalConstants.h"
#include "PlasmaFunctions.h"
#include "PlasmaState.h"
#include "TfdhSolution.h"

#include <cassert>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <vector>


namespace {

  struct IntegrationResults {
    const std::vector<double> rs;
    const std::vector<double> phis;
  };

  struct RhsParams {
    const PlasmaState& p;
  };


  int tfdhOdeRhs(const double r, const double f[], double dfdr[], void *params) {
    const double& qe = PhysicalConstantsCGS::ElectronCharge;
    const double phi = qe*f[0]/r;
    const PlasmaState& p = static_cast<RhsParams*>(params)->p;
    const double ne = Plasma::ne(phi, p);
    const double ionChargeDensity = Plasma::totalIonChargeDensity(phi, p);
    dfdr[0] = f[1];
    dfdr[1] = -4.0*M_PI*qe * r * (ionChargeDensity - ne);
    return GSL_SUCCESS;
  }


  IntegrationResults integrateODE(const Element& e, const PlasmaState& p,
      const double r_init, const double r_final, const double dv0)
  {
    const double eps_abs = 1e-6;
    const double eps_rel = 0;
    const double& qe = PhysicalConstantsCGS::ElectronCharge;
    const size_t dim = 2;
    double solution[dim] = {qe*e.Z + r_init*dv0, dv0};
    RhsParams params {p};

    gsl_odeiv2_evolve* ev = gsl_odeiv2_evolve_alloc(dim);
    gsl_odeiv2_control* ctrl = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
    gsl_odeiv2_step* step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd, dim);
    gsl_odeiv2_system sys = {tfdhOdeRhs, nullptr, dim, &params};

    // vectors in which to store (r,phi) at each step
    std::vector<double> rs = {r_init};
    std::vector<double> phis = {qe*solution[0]/r_init};

    double r = r_init;
    double dr = r_init;
    const double max_dr_over_r = 0.2;
    while (r < r_final) {
      dr = fmin(dr, max_dr_over_r * r); // prevent dr from being "too big"
      const int status = gsl_odeiv2_evolve_apply(ev, ctrl, step, &sys, &r, r_final, &dr, solution);
      assert(status==GSL_SUCCESS);
      rs.push_back(r);
      phis.push_back(qe*solution[0]/r);
      if (solution[0] <= 0 or (solution[1]-solution[0])/r > 0) break;
    }
    assert(r < r_final and "integrated ODE until final radius without terminating!");

    gsl_odeiv2_step_free(step);
    gsl_odeiv2_control_free(ctrl);
    gsl_odeiv2_evolve_free(ev);
    return {rs, phis};
  }


  double findPotentialRoot(const Element& e, const PlasmaState& p,
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
        const auto& tfdh = integrateODE(e, p, r_init, r_final, v_low);
        if (tfdh.phis.back() < 0.0) {
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

        const auto& tfdh = integrateODE(e, p, r_init, r_final, v_mid);
        ((tfdh.phis.back() >= 0.0) ? v_high : v_low) = v_mid;
      }
      assert(success and "failed to find potential root within bracket");
    }

    return v_mid;
  }
}



TfdhSolution TFDH::solve(const Element& e, const PlasmaState& p)
{
  // ODE integration bounds
  const double rws = Plasma::radiusWignerSeitz(e,p);
  const double ri = 1e-4 * rws;
  const double rf = 1e3 * rws;

  // NOTE: with this setup, the "correct" ODE is integrated twice -- first while
  // finding the correct potential, then again using the correct potential.
  // this should be a negligible cost, but could be optimized away if need be.
  const double dv0 = findPotentialRoot(e, p, ri, rf);
  const IntegrationResults& results = integrateODE(e, p, ri, rf, dv0);
  return TfdhSolution(results.rs, results.phis);
}

