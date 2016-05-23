#ifndef _atomic_constants_H
#define _atomic_constants_H

// Proton mass (g)
const double mp        = 1.67262158e-24;

// Electron mass (g)
const double me        =  9.10938291e-28;

// Speed of light (cm/s)
const double light     =  2.99792458e+10;

// electron volt in (cgs)
const double eV        =  1.60217653e-12;

// Boltzmann constant (cgs)
const double boltz     = 1.3806504e-16;

// Boltzmann constant (eV)
const double boltzEv   = 8.6173324e-5;

// Planck's constant (cgs)
const double planck    = 6.626075540e-27;

// Electron charge (cgs)
const double esu       = 4.80320427e-10;

// Atomic mass unit in grams
const double amu       = 1.660539e-24;

// Electron rest mass in MeV
const double mec2      = 0.51099891013;

// Parsec (cm)
const double pc        = 3.08567758e18;

// Bohr radius (cm)
const double a0        = 0.052917721092e-7;

// Year (sec)
const double year      = 365.242*24*3600;

// Solar mass (g)
const double msun      = 1.9891e33;

// Atomic masses
const std::vector<double> atomic_mass = {0.000549,  // 0  electron
					 1.00797,   // 1  H
					 4.00260,   // 2  He
					 6.941,     // 3  Li
					 9.01218,   // 4  Be
					 10.81,     // 5  B
					 12.011,    // 6  C
					 14.0067,   // 7  N
					 15.9994,   // 8  O
					 18.998403, // 9  F
					 20.179,    // 10 Ne
					 22.98977,  // 11 Na
					 24.305,    // 12 Mg
					 26.98154,  // 13 Al
					 28.0855,   // 14 Si
					 30.97376,  // 15 P
					 32.06,     // 16 S
					 35.453,    // 17 Cl
					 39.948,    // 18 Ar
					 39.0983,   // 19 K
					 40.08,     // 20 Ca
					 44.9559,   // 21 Sc
					 47.90,     // 22 Ti
					 50.9415,   // 23 V
					 51.996,    // 24 Cr
					 54.9380,   // 25 Mn
					 55.847,    // 26 Fe
					 58.9332,   // 27 Co
					 58.70,     // 28 Ni
					 63.546,    // 29 Cu
					 65.38 };   // 30 Zn
#endif
