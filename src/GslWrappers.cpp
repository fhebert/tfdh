
#include "GslWrappers.h"

#include <cassert>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


double GslFunction_Unpacker(const double x, void *params) {
  const GslFunction* object = static_cast<const GslFunction*>(params);
  return object->f(x);
}


double gslBrent(const GslFunction& func, const double xa, const double xb,
    const double dx_abs)
{
  // check that interval is non-empty
  // (this does not check that it is indeed bracketing a root! -- that job is
  // left to the caller of this function)
  assert(xa != xb);

  // cleverness to enable the rootfind function and data ('params' in GSL) to
  // be held together in a class (hat tip: Mark Scheel at Caltech)
  //
  // f.params points to a GslFunction object, cast to void*, which provides
  //     the interface AND the data needed to compute f(x). in this sense, it
  //     "packs" the object
  // f.function points to a wrapper which casts the params void* back into a
  //     GslFunction object and then calls the member function f(x) on it.
  //     in this sense it is "unpacking" the object.
  gsl_function f;
  f.params = const_cast<GslFunction*>(&func);
  f.function = &GslFunction_Unpacker;

  // set up solver
  gsl_root_fsolver* solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  gsl_root_fsolver_set(solver, &f, fmin(xa,xb), fmax(xa,xb));

  // iterate
  const int max_iter = 200;
  const double dx_rel = dx_abs; // will give smallest resolvable relative dx
  int status = GSL_CONTINUE;
  for (int iter=0; (status==GSL_CONTINUE) and (iter<max_iter); ++iter) {
    status = gsl_root_fsolver_iterate(solver);
    const double xlo = gsl_root_fsolver_x_lower(solver);
    const double xhi = gsl_root_fsolver_x_upper(solver);
    status = gsl_root_test_interval(xlo, xhi, dx_abs, dx_rel);
  }
  assert(status==GSL_SUCCESS);

  const double result = gsl_root_fsolver_root(solver);
  gsl_root_fsolver_free(solver);
  return result;
}


double gslQuadratureNG(const GslFunction& func, const double xi, const double xf,
    const double eps_abs, const double eps_rel)
{
  double result = 0.0;
  double abs_err = 0.0;
  size_t num_eval = 0;
  {
    gsl_function f;
    f.params = const_cast<GslFunction*>(&func);
    f.function = &GslFunction_Unpacker;
    gsl_integration_qng(&f, xi, xf, eps_abs, eps_rel, &result, &abs_err, &num_eval);
  }
  return result;
}


double gslQuadratureAG(const GslFunction& func, const double xi, const double xf,
    const double eps_abs, const double eps_rel)
{
  double result = 0.0;
  double abs_err = 0.0;
  // TODO: determine minimal sufficient size, or perhaps allocate a global
  //       workspace to reduce memory load
  const size_t max_intervals = 100;
  gsl_integration_workspace* ws = gsl_integration_workspace_alloc(max_intervals);
  {
    gsl_function f;
    f.params = const_cast<GslFunction*>(&func);
    f.function = &GslFunction_Unpacker;
    gsl_integration_qag(&f, xi, xf, eps_abs, eps_rel,
        max_intervals, GSL_INTEG_GAUSS21, ws,
        &result, &abs_err);
  }
  gsl_integration_workspace_free(ws);
  return result;
}

