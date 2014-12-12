
#ifndef TFDH_PHYSICAL_CONSTANTS_CGS_H
#define TFDH_PHYSICAL_CONSTANTS_CGS_H

#include <cmath>


namespace PhysicalConstantsCGS
{
  // nuclear physics constants
  const double ProtonMass = 1.672621e-24; // gram
  const double ElectronMass = 9.1093826e-28; // cm
  const double ElectronCharge = 4.80320441e-10; // statvolt or esu

  // other constants
  const double HPlanck = 6.6260693e-27; // erg sec
  const double Hbar = 1.05457168e-27; // erg sec
  const double GNewton = 6.6742e-8; // 'cgs'
  const double KBoltzmann = 1.3806505e-16; // erg / kelvin
  const double SpeedLight = 2.99792458e10; // cm / s

  // useful combinations
  const double MeCC = ElectronMass * SpeedLight * SpeedLight;
  const double MpCC = ProtonMass * SpeedLight * SpeedLight;
  const double NePrefactor = pow(2.0*ElectronMass, 1.5)
                             / (2.0 * M_PI*M_PI * Hbar*Hbar*Hbar);

  // useful physical scales
  const double BohrRadius = 5.2917721e-9; // cm
  const double SolarMass = 1.98892e33; // gram
  const double SolarRadius = 6.961e10; // cm
}


#endif // TFDH_PHYSICAL_CONSTANTS_CGS_H
