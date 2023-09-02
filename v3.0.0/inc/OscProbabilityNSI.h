/**
 *
 *  @file OscProbabilityNSI.h
 *
 *  @brief It deals with Physics associated with Neutral Current Forward Elastic Non-Standard fermionic Interactions.
 *
 *  It provides a unique opportunity to see how the oscillation probabilities get modified with the NCNSI Physics Lists. <br>
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #ifndef _OscProbabilityNSI_H
    #define _OscProbabilityNSI_H

    #include <cmath>
    #include <OscProbability.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Class Definations *//

    class OscProbabilityNSI : virtual OscProbabilitySTD
  {
    private:
    //* ....oooOOOOO NSI Settings OOOOOooo.... *//
    //* NC-NSI Parameters *//
    std::vector < SetBSMIntParams > NCNSI_E;
    std::vector < SetBSMIntParams > NCNSI_U;
    std::vector < SetBSMIntParams > NCNSI_D;

    public:
    virtual ~OscProbabilityNSI() {};
    virtual void SetNCNSI_E(SetBSMIntParams ncE[], size_t _nTe_);
    virtual void SetNCNSI_U(SetBSMIntParams ncU[], size_t _nTu_);
    virtual void SetNCNSI_D(SetBSMIntParams ncD[], size_t _nTd_);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > SetNSIMatrixE(int type, double rho);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > SetNSIMatrixU(int type, double rho);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > SetNSIMatrixD(int type, double rho);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > SetNCNSIMatrix(int typ, double rho);

    //* Oscillation NCNSI Probability Functions *//
    virtual double** EarthMatOscProbNCNSI(int type, double Erg);
    virtual double EarthMatOscProbNCNSI(int ij, int jk, int type, double Erg);                             /// Earth Matter Oscilation Probabilty;
    virtual double** ProfiledMatOscProbNCNSI(int type, double Erg);
    virtual double ProfiledMatOscProbNCNSI(int ij, int jk, int type, double Erg);                          /// Profiled Matter Oscilation Probabilty;
    virtual double** ConstMatOscProbNCNSI(int type, double Erg, double dens, double length);
    virtual double ConstMatOscProbNCNSI(int ij, int jk, int type, double Erg, double dens, double length); /// Const-Dens-Mat Oscillation Probability;

    //* Get Oscillation Parameters with NCNSI *//
    virtual double** GetConstMatNSIOscParms(int type, double dens, double Erg, const char *wh);            /// Constant Matter Modified Matrices;
  };

    #endif

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//





