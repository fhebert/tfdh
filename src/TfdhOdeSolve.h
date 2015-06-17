
#ifndef TFDH_ODE_SOLVE_HPP
#define TFDH_ODE_SOLVE_HPP

#include "TfdhSolution.h"

#include <vector>

class Element;
class PlasmaState;


namespace TFDH {
  TfdhSolution solve(const Element& e, const PlasmaState& p);
}


#endif // TFDH_ODE_SOLVE_HPP
