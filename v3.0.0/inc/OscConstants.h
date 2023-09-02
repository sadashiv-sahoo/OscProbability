/**
 *
 *  @file  OscConstants.h
 *  @brief Introductions of Physical Constants
 *
 *  This contains the values of physical parameters
 *  used throughout this programme as Global variables.
 *  The values are considered in are taken from : <br>
 *  Particle Data Group <a href="https://journals.aps.org/prd/pdf/10.1103/PhysRevD.98.030001"> [Phys. Rev. D 98, 030001 (2018)] </a>
 *  Chapter 1. Physical Constants, Page-127
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #ifndef _OscConstants_H
    #define _OscConstants_H

    #include <cmath>
    #include <complex>
    #include <iostream>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//

    const static int       NoMix     = 1;
    const static double    GeV       = 1.0E9;
    const static double    deg       = 180./M_PI;                   // rad to deg
    const static double    rad       = M_PI/180.;                   // deg to rad
    const static double    h         = 6.626070040E-34;             // (j-s) Planck's Constant
    const static double    c         = 299792458.0;                 // (m/s) Speed of Light in vacuum 
    const static double    eV        = 1.6021766208E-19;            // (eV)
    const static double    GF        = 1.1663787E-23;               // (eV^-2)  Fermi Coupling Constant
    const static double    Na        = 6.022140857E23;              // (mol^-1) Avogrado Constant
    const static double    hbar      = h/(2.*M_PI);                 // h/2Pi;
    const static double    HbarC     = (1E9*hbar*c)/eV;             // (GeV) [(1.0x10^9 x hbar x c)/eV]
    const static double    Per_Sec   = (hbar/eV);                   // 1/s [6.58211899x10^(-16) eV]
    const static double    km2_eV    = 1E12/HbarC;                  // km (eV^-1)
    const static double    cm2_eV    = km2_eV/1E5;                  // cm (eV^-1)
    const static double    cc2_eV3   = pow(cm2_eV,3);               // cc (eV^-3)
    const static double    sqrt2GF   = (sqrt(2)*(GF)*Na)/(cc2_eV3); // (g/cc)*eV [7.63247E-14]
    const static double    GFbysqrt2 = sqrt2GF/2.;                  // (g/cc)*eV [3.81623E-14]
    const static double    sinsqrThw = 0.23122;                     // sin^2(theta_W) [Weinberg Angle]
    const static double    Cv_by_Ca  = 1.0 - 4*sinsqrThw;           // Cv/Ca = 1 - 4*sin^2(theta_W)
    const static double    rS        = 3.00000000;                  // (up_quark_No.dens)/(electron_No.dens)  [Valency Quarks]
    const static std::complex <double> img(0.,1.);                  // Imaginary Number i = sqrt(-1)

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 

  #endif






