
    //////////////////////////////////////////////////////////////////
    //   .........oooOOO 18th-Aug-2020 Last Update OOOooo........   //
    //                        SADASHIV SAHOO                        //
    //               India-based Neutrino Observatory               //
    //        Homi Bhabha National Institute, Mumbai, INDIA         // 
    //          @email: sadashiv.sahoo@tifr.res.in                  //
    //                : sadashiv.sahoo@iopb.res.in                  //
    //////////////////////////////////////////////////////////////////

    #ifndef _NeutrinoOscProbSTD_H
    #define _NeutrinoOscProbSTD_H

    #include <OscProbability.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 
    //* Class Definations *//

    class NeutrinoOscProbSTD {

    protected:

    //* ....oooOOOOO Default Global Settings OOOOOooo.... *// 
    //* Mass Square Differences *//
    double m21;              // Sol in eV scale;
    double m31;              // Atm in eV scale;

    //* Mixing Angles Variables *//
    double theta12;          // in degree scale;
    double theta23;          // in degree sacle;
    double theta13;          // in degree scale;
    double Delta_cp;         // in degree scale;
    complex <double> cpterm; // Calculates cos(Delta_cp) + i Sin(Delta_cp);

    //* Earth Variables *//                     // All Scale in KMs;
    double Depth;
    double AtmosHeight;
    double EarthRadius;

    //* Earth Density Variables *//
    vector < pair < double, double > > PREM;    // PREM Global Fixed   Data Store;
    vector < pair < double, double > > Out_In;  // PREM Local  Dynamic Data Store;
    vector < pair < double, double > > Profile; // User's      Fixed   Data Store;

    public:

    virtual ~NeutrinoOscProbSTD() {};
    void SetEarthDensityProfile();                                      // Loads Earth's PREM Profile for Atmospheric Expt.;
    void SetGivenDensityProfile(double profile[][2], size_t layers);    // Loads Fixed Baseline Expt.;

    //* ....oooOOOOO Default Global Functions OOOOOooo.... *//
    //* Mandatory Set Functions *//
    Matrix<complex<double>,dim,dim> SetPMNS();                           // Sets PMNS Matrix;
    Matrix<complex<double>,dim,dim> SetMassHamiltonian(double Erg);      // Hamiltonian in Mass-Basis; 
    void SetMassSqDiff(double Sol, double Atm);                          
    void SetMixAngles(double Th12, double Th23, double Th13, double CP); 

    //* EarthMatOscProb Set Functions *//
    void SetZenith(double CosZ);                                          
    void SetDetDepth(double depth);                                       // Includes Detector Depth;
    void SetEarthAtmos(double AtmHeight);                                 // Includes Atmospheric Production Height;

    //* ....oooOOOOO Standard Functions OOOOOooo.... *//
    //* Get Traverse Distance Functions *//                                     
    double GetTrav_Length(double cosZ);                                         // Atmospheric Baseline(s);
    double GetTrav_Length(double profile[][2], size_t layers);                  // Fixed Baseline

    //* Get Modified Mixing Angles *//     
    double GetEarthMixAngle(int i, int j, int type, double Erg);                // Earth Matter Modified Mixing Angles;
    double GetProfileMixAngle(int i, int j, int type, double Erg);              // Profile Matter Modified Mixing Angles;
    double GetConstMatMixAngle(int i, int j, int type, double dens, double Erg);// Constant Density Matter Modified Mixing Angles;

    //* Get Modified Mass Squared Diffrences *//
    double GetEarthModMass(int i, int j, int type, double Erg);                 // Earth Matter Modified Mass-Sqr Differences;
    double GetProfileModMass(int i, int j, int type, double Erg);               // Profile Matter Modified Mass-Sqr Differences;
    double GetConstMatModMass(int i, int j, int type, double dens, double Erg); // Constant Matter Modified Mass-Sqr Differences;

    //* Oscillation Probability Functions *//
    double EarthMatOscProb(int ij, int jk, int type, double Erg);                             // Earth Matter Oscilation Probabilty;
    double ProfiledMatOscProb(int ij, int jk, int type, double Erg);                          // Profiled Matter Oscilation Probabilty;
    double VaccumOscProb(int ij, int jk, int type, double Erg, double length);                // Vaccum Oscillation Probability;
    double ConstMatOscProb(int ij, int jk, int type, double Erg, double dens, double length); // Const-Dens-Mat Oscillation Probability;

  }; 

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 

    #endif







