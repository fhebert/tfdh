
#include "GslWrappers.h"

#include <cassert>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


namespace {
  // free function to interface between the gsl_function interface and the
  // more user-friendly FunctionObject class.
  //
  // f.params points to a FunctionObject, cast to void*, which provides the
  //     interface AND the data needed to compute f(x)
  // f.function calls this wrapper, which casts back to FunctionObject and
  //     calls the member function f(x)
  inline double callFunctionFromObject(const double x, void *params) {
    const GSL::FunctionObject* func = static_cast<const GSL::FunctionObject*>(params);
    return (*func)(x);
  }
}


GSL::Spline::Spline(const std::vector<double>& x, const std::vector<double>& f)
: acc(gsl_interp_accel_alloc()),
  spline(gsl_spline_alloc(gsl_interp_cspline, x.size()))
{
  assert(x.size()==f.size());
  gsl_spline_init(spline, x.data(), f.data(), x.size());
}

GSL::Spline::~Spline()
{
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}

double GSL::Spline::eval(const double r) const
{
  return gsl_spline_eval(spline, r, acc);
}


double GSL::findRoot(const GSL::FunctionObject& func, const double xa, const double xb,
    const double eps_abs, const double eps_rel)
{
  // check that interval brackets the root
  assert(func(xa)*func(xb) <= 0);
  if (xa==xb) return xa;

  // set up solver
  gsl_function f;
  f.params = const_cast<GSL::FunctionObject*>(&func);
  f.function = &callFunctionFromObject;

  gsl_root_fsolver* solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  gsl_root_fsolver_set(solver, &f, fmin(xa,xb), fmax(xa,xb));

  // iterate
  const int max_iter = 80;
  int status = GSL_CONTINUE;
  for (int iter=0; (status==GSL_CONTINUE) and (iter<max_iter); ++iter) {
    status = gsl_root_fsolver_iterate(solver);
    const double xlo = gsl_root_fsolver_x_lower(solver);
    const double xhi = gsl_root_fsolver_x_upper(solver);
    status = gsl_root_test_interval(xlo, xhi, eps_abs, eps_rel);
  }
  assert(status==GSL_SUCCESS);

  const double result = gsl_root_fsolver_root(solver);
  gsl_root_fsolver_free(solver);
  return result;
}


double GSL::integrate(const GSL::FunctionObject& func, const double xi, const double xf,
    const double eps_abs, const double eps_rel)
{
  // set up integrator
  gsl_function f;
  f.params = const_cast<GSL::FunctionObject*>(&func);
  f.function = &callFunctionFromObject;

  // NOTE: this choice has worked so far, but could be increased if necessary
  const size_t max_intervals = 100;
  gsl_integration_workspace* ws = gsl_integration_workspace_alloc(max_intervals);

  double result = 0.0;
  double abs_err = 0.0;
  gsl_integration_qag(&f, xi, xf, eps_abs, eps_rel,
      max_intervals, GSL_INTEG_GAUSS21, ws,
      &result, &abs_err);

  gsl_integration_workspace_free(ws);
  return result;
}

