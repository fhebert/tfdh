
#ifndef TFDH_ODE_SOLVE_H
#define TFDH_ODE_SOLVE_H

#include "TfdhSolution.h"

#include <vector>

class Element;
class PlasmaState;


namespace TFDH {
  TfdhSolution solve(const Element& e, const PlasmaState& p);
}


#endif // TFDH_ODE_SOLVE_H
