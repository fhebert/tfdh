
#ifndef TFDH_TFDH_FUNCTIONS_H
#define TFDH_TFDH_FUNCTIONS_H

struct PlasmaState;
struct RadialFunction;


namespace TFDH {

  double neBound(double phi, const PlasmaState& p);
  RadialFunction neBound(const RadialFunction& tfdh, const PlasmaState& p);

}


#endif // TFDH_TFDH_FUNCTIONS_H
