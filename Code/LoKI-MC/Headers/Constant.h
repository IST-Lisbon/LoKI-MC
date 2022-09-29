#ifndef __Constant__
#define __Constant__
#include <cmath>

namespace Constant{
	const double boltzmann = 1.38064852e-23;           									// Boltzmann constant in J/K
	const double electronCharge = 1.6021766208e-19;    									// Electron charge in C
	const double electronMass = 9.10938356e-31;        									// Electron mass in Kg
	const double unifiedAtomicMass = 1.660539040e-27;  									// Unified Atomic Mass unit (UAM) in kg
	const double bohrRadius = 5.2917721067e-11;        									// Bohr radius in m
	const double vacuumPermittivity = 8.854187817e-12; 									// Vacuum permittivity in F/m
	const double vacuumPermeability = 1.25663706212e-6;									// Vacuum permeability in N A^(-2)
	const double planck = 6.626070040e-34;            									// Plank constant in J s
	const double planckReduced = planck/(2.0*M_PI);
	const double speedOfLight = 299792458;             									// Speed of light in vacuum in m/s
	const double atmosphereInPa = 101325;              									// Standard atmosphere in Pa
	const double atmosphereInTorr = 760;               									// Standard atmosphere in Torr (not obtained from NIST database)
	const double avogadro = 6.02214076e23;												// Avogadro's constant
	const double boltzmannInEV = boltzmann/electronCharge; 								// Boltzmann constant in eV/K
	const double planckReducedInEV = planck/(2.0*M_PI*electronCharge);
	const double idealGas = boltzmann*avogadro;											// ideal gas constant in J K^(-1) mol^(-1)

	const double NON_DEF = -123456789;
};

#endif
