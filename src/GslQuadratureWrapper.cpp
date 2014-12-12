
#include "GslQuadratureWrapper.h"
#include "GslFunction.h"

#include <gsl/gsl_integration.h>


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

