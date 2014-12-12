
#ifndef TFDH_GSL_QUADRATURE_WRAPPER_H
#define TFDH_GSL_QUADRATURE_WRAPPER_H

#include "GslFunction.h"


double gslQuadratureNG(const GslFunction& func, double x1, double x2,
    double eps_aps, double eps_rel);


#endif // TFDH_GSL_QUADRATURE_WRAPPER_H
