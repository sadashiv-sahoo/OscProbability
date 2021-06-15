
    //////////////////////////////////////////////////////////////////
    //   .........oooOOO 18th-Aug-2020 Last Update OOOooo........   //
    //                        SADASHIV SAHOO                        //
    //               India-based Neutrino Observatory               //
    //        Homi Bhabha National Institute, Mumbai, INDIA         // 
    //          @email: sadashiv.sahoo@tifr.res.in                  //
    //                : sadashiv.sahoo@iopb.res.in                  //
    //////////////////////////////////////////////////////////////////

    #ifndef _NeutrinoOscProbSME_H
    #define _NeutrinoOscProbSME_H

    #include <OscProbability.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 
    //* Lorentz Invariance Violation in Standard Model Extension Framework *//

    class NeutrinoOscProbSME : virtual NeutrinoOscProbSTD {

    private:

    //* ....oooOOOOO Gauge Invariant LIV Settings OOOOOooo.... *//
    //* GILIV Control Parameter factors *//
    double aL;
    double cL;

    //* LIV Gauge-Symmetry Elements *//
    //* LIV ODD(a) *//
    complex <double> aX_ee;         // OscLIV in Electron to Electron Channel 
    complex <double> aX_mumu;       // OscLIV in Muon to Muon Channel 
    complex <double> aX_tautau;     // OscLIV in Tauon to Tauon Channel 
    complex <double> aX_emu;        // OscLIV in Electron to Muon Channel 
    complex <double> aX_etau;       // OscLIV in Electron to Tauon Channel 
    complex <double> aX_mutau;      // OscLIV in Muon to Tauon Channel 

    //* LIV EVEN(c) *//
    complex <double> cXX_ee;        // OscLIV in Electron to Electron Channel
    complex <double> cXX_mumu;      // OscLIV in Muon to Muon Channel
    complex <double> cXX_tautau;    // OscLIV in Tauon to Tauon Channel
    complex <double> cXX_emu;       // OscLIV in Electron to Muon Channel
    complex <double> cXX_etau;      // OscLIV in Electron to Tauon Channel
    complex <double> cXX_mutau;     // OscLIV in Muon to Tauon Channel

    //* ....oooOOOOO Gauge Violating LIV Settings OOOOOooo.... *//
    //* GVLIV Control Parameter factors *//
    double _h_;
    double _g_;
    complex <double> epsilon;

    //* GVLIV Gauge-Symmetry-Violating Elements *//
    //* GVLIV ODD(g) *//
    complex <double> gXXX_eebar;          // OscGVLIV in Electron to Electron Channel 
    complex <double> gXXX_mumubar;        // OscGVLIV in Muon to Anti-Muon Channel 
    complex <double> gXXX_tautaubar;      // OscGVLIV in Tauon to Anti-Tauon Channel 
    complex <double> gXXX_emubar;         // OscGVLIV in Electron to Anti-Muon Channel 
    complex <double> gXXX_etaubar;        // OscGVLIV in Electron to Anti-Tauon Channel 
    complex <double> gXXX_mutaubar;       // OscGVLIV in Muon to Anti-Tauon Channel 

    //* GVLIV EVEN(h) *//
    complex <double> hXX_eebar;           // OscGVLIV in Electron to Anti-Electron Channel
    complex <double> hXX_mumubar;         // OscGVLIV in Muon to Anti-Muon Channel
    complex <double> hXX_tautaubar;       // OscGVLIV in Tauon to Anti-Tauon Channel
    complex <double> hXX_emubar;          // OscGVLIV in Electron to Anti-Muon Channel
    complex <double> hXX_etaubar;         // OscGVLIV in Electron to Anti-Tauon Channel
    complex <double> hXX_mutaubar;        // OscGVLIV in Muon to Anti-Tauon Channel
        
    public:

    virtual ~NeutrinoOscProbSME() {};

    //* ....oooOOOOO Standard Model Extension Lorentz Invariace Violation OOOOOooo.... *//
    //* SME-GILIV Set Functions *//
    void SetGILIVControls(int _mu_, int _nu_);                                       // Momentum-Four-Vector indices;        
    void SetGILIVaX( complex <double> _ee, complex <double> _mumu,
                       complex <double> _tautau, complex <double> _emu, 
                           complex <double> _etau, complex <double> _mutau );        // CPT-Odd Parameters
    void SetGILIVcXX( complex <double> _ee_, complex <double> _mumu_,
                        complex <double> _tautau_, complex <double> _emu_, 
                            complex <double> _etau_, complex <double> _mutau_ );     // CPT-Even-LIV Parameters 
    Matrix<complex<double>,dim,dim> SetGILIVMatrixA();                               // LIV-ODD  Matrix  in GeV Unit 
    Matrix<complex<double>,dim,dim> SetGILIVMatrixC();                               // LIV-EVEN Matrix  Unit Less 
    Matrix<complex<double>,dim,dim> SetGILIVMatrixA(int type);                       // LIV-ODD  Matrix* for antiparticle
    Matrix<complex<double>,dim,dim> SetGILIVMatrixC(int type);                       // LIV-EVEN Matrix* for antiparticle

    //* SME-GILIV Get Mass Modifies Functions *//
    double GetVacGILIVModMass(int i, int j, int type, double Erg);                   // Vac LIV Modified Mass;
    double GetEarthGILIVModMass(int i, int j, int type, double Erg);                 // Earth Matter LIV Modified Mass;
    double GetProfileGILIVModMass(int i, int j, int type, double Erg);               // Profile Matter LIV Modified Mass;
    double GetConstMatGILIVModMass(int i, int j, int type, double dens, double Erg); // Const Density Matter LIV Modified Mass;

    //* SME-GILIV Get Mix-Angles Modifies Functions *//
    double GetVacGILIVMixAngle(int i, int j, int type, double Erg);                  // Vac LIV Modified Mixing Angles;
    double GetEarthGILIVMixAngle(int i, int j, int type, double Erg);                // Earth Matter LIV Modified Mixing Angles;
    double GetProfileGILIVMixAngle(int i, int j, int type, double Erg);              // Profile Matter LIV Modified Mixing Angles;
    double GetConstMatGILIVMixAngle(int i, int j, int type, double dens, double Erg);// Const Density Matter LIV Modified Mixing Angles;

    //* SME-GILIV Oscillation Probability Functions *//
    double EarthMatOscProbGILIV(int ij, int jk, int type, double Erg);                             // LIV Earth Matter Oscilation Probabilty;
    double ProfiledMatOscProbGILIV(int ij, int jk, int type, double Erg);                          // LIV Profiled Matter Oscilation Probabilty;
    double VaccumOscProbGILIV(int ij, int jk, int type, double Erg, double length);                // LIV Vaccum Oscillation Probability;
    double ConstMatOscProbGILIV(int ij, int jk, int type, double Erg, double dens, double length); // LIV Const-Dens-Mat Oscillation Probability;

    //* SME-GVLIV Set Functions *//
    void SetGVLIVControls(int _mu_, int _sigma_);                                     // Momentum-Four-Vector indices;        
    void SetGVLIVPolzation(complex <double> eps);                                     // GVLIV Vacuum Polarization;
    void SetGVLIVgXXX(complex <double> _eebar, complex <double> _mumubar,
                         complex <double> _tautaubar, complex <double> _emubar, 
                            complex <double> _etaubar, complex <double> _mutaubar);   // CPT-Odd Parameters
    void SetGVLIVhXX(complex <double> _eebar_, complex <double> _mumubar_,
                        complex <double> _tautaubar_, complex <double> _emubar_, 
                           complex <double> _etaubar_, complex <double> _mutaubar_);  // CPT-Even-GVLIV Parameters

    Matrix <complex<double>,dim,dim> SetGVLIVMatrixC();                               // Charge-Conjugate-Matrix;
    Matrix <complex<double>,dim,dim> SetGVLIVMatrixG();                               // GVLIV-ODD  Matrix  in GeV Unit 
    Matrix <complex<double>,dim,dim> SetGVLIVMatrixH();                               // GVLIV-EVEN Matrix  Unit Less 

    //* SME-GVLIV Get Functions *//
    Matrix <complex<double>,dim,dim> GetGVLIVHamiltonian(int ij, int jk, double Erg); // GVLIV Partial Hamitonian;

    //* SME-GVLIV Oscillation Probability Functions *//
    double VaccumOscProbGVLIV(int ij, int jk, double Erg, double length);             // GVLIV Vaccum Oscillation Probability;

  }; 

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 

    #endif






