
    //////////////////////////////////////////////////////////////////
    //   .........oooOOO 18th-Aug-2020 Last Update OOOooo........   //
    //                        SADASHIV SAHOO                        //
    //               India-based Neutrino Observatory               //
    //        Homi Bhabha National Institute, Mumbai, INDIA         // 
    //          @email: sadashiv.sahoo@tifr.res.in                  //
    //                : sadashiv.sahoo@iopb.res.in                  //
    //////////////////////////////////////////////////////////////////

    #ifndef _NeutrinoOscProbNSI_H
    #define _NeutrinoOscProbNSI_H

    #include <OscProbability.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 
    //* Non-Standard Interaction *//

    class NeutrinoOscProbNSI : virtual NeutrinoOscProbSTD {

    private:

    //* ....oooOOOOO NSI Settings OOOOOooo.... *//
    //* NSI Parameters for electron *//
    complex <double> nsiE_ee;        // NSI-Electron-ee
    complex <double> nsiE_mumu;      // NSI-Electron-MuMu
    complex <double> nsiE_tautau;    // NSI-Electron-TauTau
    complex <double> nsiE_emu;       // NSI-Electron-eMu
    complex <double> nsiE_etau;      // NSI-Electron-eTau
    complex <double> nsiE_mutau;     // NSI-Electron-MuTau

    //* NSI Parameters for u-quark *//
    complex <double> nsiU_ee;        // NSI-U-quark-ee
    complex <double> nsiU_mumu;      // NSI-U-quark-MuMu
    complex <double> nsiU_tautau;    // NSI-U-quark-TauTau
    complex <double> nsiU_emu;       // NSI-U-quark-eMu
    complex <double> nsiU_etau;      // NSI-U-quark-eTau
    complex <double> nsiU_mutau;     // NSI-U-quark-MuTau

    //* NSI Parameters for d-quark *//
    complex <double> nsiD_ee;        // NSI-D-quark-ee
    complex <double> nsiD_mumu;      // NSI-D-quark-MuMu
    complex <double> nsiD_tautau;    // NSI-D-quark-TauTau
    complex <double> nsiD_emu;       // NSI-D-quark-eMu
    complex <double> nsiD_etau;      // NSI-D-quark-eTau
    complex <double> nsiD_mutau;     // NSI-D-quark-MuTau

    public:

    //* ....oooOOOOO Non Standard Interactions  OOOOOooo.... *//
    //* NSI Set Functions *//
    void SetNSI_E( complex <double> E_ee_, complex <double> E_mumu_,
                       complex <double> E_tautau_, complex <double> E_emu_, 
                           complex <double> E_etau_, complex <double> E_mutau_ );       // NSI-with-Electron Parameters
    void SetNSI_U( complex <double> U_ee_, complex <double> U_mumu_,
                       complex <double> U_tautau_, complex <double> U_emu_, 
                           complex <double> U_etau_, complex <double> U_mutau_ );       // NSI-with-U-quark Parameters
    void SetNSI_D( complex <double> D_ee_, complex <double> D_mumu_,
                       complex <double> D_tautau_, complex <double> D_emu_, 
                           complex <double> D_etau_, complex <double> D_mutau_ );       // NSI-with-D-quark Parameters

    Matrix<complex<double>,dim,dim> SetNSIMatrixE();                                    // NSI-with-Electron
    Matrix<complex<double>,dim,dim> SetNSIMatrixU();                                    // NSI-with-U-quark
    Matrix<complex<double>,dim,dim> SetNSIMatrixD();                                    // NSI-with-D-quark

    //* NSI Get Functions *//
    //* NSI Modified Mass *//
    double GetEarthNSIModMass(int i, int j, int type, double Erg);                      // Earth Matter NSI Modified Mass;
    double GetProfileNSIModMass(int i, int j, int type, double Erg);                    // Profile Matter NSI Modified Mass;
    double GetConstMatNSIModMass(int i, int j, int type, double dens, double Erg);      // Constant Density Matter NSI Modified Mass;

    //* NSI Modified Mix-Angles *//
    double GetEarthNSIMixAngle(int i, int j, int type, double Erg);                     // Earth Matter NSI Modified Mixing Angles;
    double GetProfileNSIMixAngle(int i, int j, int type, double Erg);                   // Profile Matter NSI Modified Mixing Angles;
    double GetConstMatNSIMixAngle(int i, int j, int type, double dens, double Erg);     // Constant Density Matter NSI Modified Mixing Angles;

    //* NSI Oscillation Probability Functions *//
    double EarthMatOscProbNSI(int ij, int jk, int type, double Erg);                             // NSI Earth Matter Oscilation Probabilty;
    double ProfiledMatOscProbNSI(int ij, int jk, int type, double Erg);                          // NSI Profiled Matter Oscilation Probabilty;
    double ConstMatOscProbNSI(int ij, int jk, int type, double Erg, double dens, double length); // NSI Const-Dens-Mat Oscillation Probability;

  }; 

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 

    #endif






