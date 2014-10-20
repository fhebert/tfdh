
#include "IntegrateTFDH.h"

#include "Element.h"
#include "PlasmaFunctions.h"
#include "PlasmaState.h"
#include "PhysicalConstants.h"
#include "RadialFunction.h"

#include <cmath>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include <iostream>

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


RadialFunction TFDH::integrate(const Element& e, const PlasmaState& p)
{
  const double dr_start = 1e-6;
  const double eps_abs = 1e-6;
  const double eps_rel = 0;
  TfdhParams params(e,p);
  gsl_odeiv2_system sys = {tfdhOde, nullptr, 2, &params};
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys,
      gsl_odeiv2_step_rk8pd, dr_start, eps_abs, eps_rel);

  // TODO -- clean all this up
  const double r_init = 1e-10;
  const double r_final = 1e-6;
  double r = r_init;
  double solution[2] = {1,2};

  std::vector<double> radii, potentials;

  const int numsteps = 100;
  for (int i=0; i<numsteps; ++i) {
    const double ri = r_init + (i+1)*(r_final-r_init)/numsteps;
    const int status = gsl_odeiv2_driver_apply(d, &r, ri, solution);
    assert(status==GSL_SUCCESS);
    radii.push_back(ri);
    potentials.push_back(solution[0]/ri);
  }

  return RadialFunction(radii, potentials);
}
