
    //////////////////////////////////////////////////////////////////
    //   .........oooOOO 18th-Aug-2020 Last Update OOOooo........   //
    //                        SADASHIV SAHOO                        //
    //               India-based Neutrino Observatory               //
    //        Homi Bhabha National Institute, Mumbai, INDIA         // 
    //          @email: sadashiv.sahoo@tifr.res.in                  //
    //                : sadashiv.sahoo@iopb.res.in                  //
    //////////////////////////////////////////////////////////////////

    #ifndef _NeutrinoOscProbDEC_H
    #define _NeutrinoOscProbDEC_H

    #include <OscProbability.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 
    //* Class Definations *//

    class NeutrinoOscProbDEC : virtual NeutrinoOscProbSTD {

    private:

    //* ....oooOOOOO Non Standard Neutrino Decay OOOOOooo.... *//
    //*  NS NuDecay Parameters *//
    complex <double> _a_11_;  // alpha_11-Mass Basis;
    complex <double> _a_12_;  // alpha_12-Mass Basis;
    complex <double> _a_13_;  // alpha_13-Mass Basis;
    complex <double> _a_21_;  // alpha_21-Mass Basis;
    complex <double> _a_22_;  // alpha_22-Mass Basis;
    complex <double> _a_23_;  // alpha_23-Mass Basis;
    complex <double> _a_31_;  // alpha_31-Mass Basis;
    complex <double> _a_32_;  // alpha_32-Mass Basis;
    complex <double> _a_33_;  // alpha_33-Mass Basis;

    public:

    //* ....oooOOOOO Non Standard Neutrino Decay  OOOOOooo.... *//
    //* NS NuDecay Set Functions *//                                                     
    void SetNuDecay( complex <double> a_11, complex <double> a_12, complex <double> a_13,        // | a11 | a12 | a13 |; NS NuDecay
                       complex <double> a_21, complex <double> a_22, complex <double> a_23,      // | a21 | a22 | a23 |; Parameters
                         complex <double> a_31, complex <double> a_32, complex <double> a_33 );  // | a31 | a32 | a33 |; (On Mass Basis)
    Matrix<complex<double>,dim,dim> SetNuDecayMatrix(double Erg);                                // NS NuDecay Matrix; (On Mass Basis)

    //* NS NuDecay Oscillation Probability Functions *//
    double EarthMatOscProbNuDecay(int ij, int jk, int type, double Erg);                             // LIV Earth Matter Oscilation Probabilty;
    double ProfiledMatOscProbNuDecay(int ij, int jk, int type, double Erg);                          // LIV Profiled Matter Oscilation Probabilty;
    double VaccumOscProbNuDecay(int ij, int jk, int type, double Erg, double length);                // LIV Vaccum Oscillation Probability;
    double ConstMatOscProbNuDecay(int ij, int jk, int type, double Erg, double dens, double length); // LIV Const-Dens-Mat Oscillation Probability;

  }; 

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 

    #endif

