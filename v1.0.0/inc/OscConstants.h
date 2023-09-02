
    //////////////////////////////////////////////////////////////////
    //   .........oooOOO 18th-Aug-2020 Last Update OOOooo........   //
    //                        SADASHIV SAHOO                        //
    //               India-based Neutrino Observatory               //
    //        Homi Bhabha National Institute, Mumbai, INDIA         // 
    //          @email: sadashiv.sahoo@tifr.res.in                  //
    //                : sadashiv.sahoo@iopb.res.in                  //
    //////////////////////////////////////////////////////////////////

    #ifndef _OscConstants_H
    #define _OscConstants_H

    #include <cmath>
    #include <complex>
    #include <iostream>

    using namespace std;

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//  

    //* Particle Data Group [Phys. Rev. D 98, 030001 (2018)] *//
    //* (Chapter 1.Physical Constants) Page-127 *//
    const static int       dim       = 3;                           // 3 Flavour Neutrinos
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
    const static double    GFbysqrt2 = sqrt2GF/2.0;                 // (g/cc)*eV [3.81623E-14]
    const static double    eS        = 0.500000000;                 // CC density fraction                    (Tentative Value)
    const static double    nS        = 0.500000000;                 // NC density fraction                    (Tentative Value)
    const static double    rS        = 3.000000000;                 // (up_quark_density)/(electron_density)  (Tentative Value) 
    const complex <double> img(0.,1.);                              // Imaginary Number i = sqrt(-1)

  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 

  #endif






