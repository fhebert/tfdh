
#ifndef TFDH_ODE_SOLVE_HPP
#define TFDH_ODE_SOLVE_HPP

class Element;
struct PlasmaState;
struct TfdhSolution;


namespace TFDH {

  TfdhSolution solve(const Element& e, const PlasmaState& p);
  TfdhSolution integrateODE(const Element& e, const PlasmaState& p,
      double r_init, double r_final, double dv0);
  double findPotentialRoot(const Element& e, const PlasmaState& p,
      double r_init, double r_final);

}


#endif // TFDH_ODE_SOLVE_HPP
