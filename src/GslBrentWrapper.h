
#ifndef TFDH_GSL_BRENT_WRAPPER_H
#define TFDH_GSL_BRENT_WRAPPER_H

#include "GslFunction.h"


double gslBrent(const GslFunction& func, double x1, double x2, double dx_abs);


#endif // TFDH_GSL_BRENT_WRAPPER_H
