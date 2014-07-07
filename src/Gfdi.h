
#ifndef TFDH_GFDI_H
#define TFDH_GFDI_H


// the int value associated with each enum item is the index used
// in the gfdi() function for accessing the tabulated coefficients
// used in a function evaluation at this order
enum class GFDI {Order12=0, Order32=1, Order52=2};


// The function gfdi() evaluates an analytic approximation to a
// "generalized Fermi-Dirac interal", from which comes the acronym.
// The analytic approximation used is taken from the paper:
//  "Equation of state of fully-ionized electron-ion plasmas"
//  Chabrier & Potekhin, 1998, Phys Rev E
//
// order - order of the generalized Fermi-Dirac integral
// chi - degeneracy parameter, mu/kT
// tau - kT/mcc
double gfdi(GFDI order, double chi, double tau);


#endif // TFDH_GFDI_H
