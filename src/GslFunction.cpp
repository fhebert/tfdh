
#include "GslFunction.h"

double GslFunction_Unpacker(const double x, void *params) {
  const GslFunction* object = static_cast<const GslFunction*>(params);
  return object->f(x);
}

