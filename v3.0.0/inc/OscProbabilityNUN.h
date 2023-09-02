/**
 *
 *  @file OscProbabilityNUN.h
 *
 *  @brief It deals with Physics associated with Indirect Sterile neutrino Searching, or Simply Non Unitarity Mixing.
 *
 *  It provides a unique opportunity to see how the oscillation probabilities get modified with the Non-Unitarity Physics Lists. <br>
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #ifndef _OscProbabilityNUN_H
    #define _OscProbabilityNUN_H

    #include <OscProbability.h>
    #include <OscProbabilityNUN.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Class Definations *//

    class OscProbabilityNUN : virtual OscProbabilitySTD
  {
    private:
    bool _NonNrm_Status_ = true;
    bool _NonUni_Status_ = false;
    double _Formulation_ = +1.00;
    std::vector < SetBSMIntParams > _NON_UNITARY_INFO_;

    public:
    virtual ~OscProbabilityNUN() {};

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Non-Unitarity Settings *//
    virtual void SetNonUnitarityPhys(bool _NsTs_);
    virtual void SetNonUnitarityNorms(bool _NsNs_);
    virtual void SetNonUnitControl(SetBSMIntParams nUn[], size_t _nun_, double _sc_);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > SetNonUnitaryPMNS(int _type_);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic >
            ZeroLengthNorm(const Eigen::Ref < const Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > > &NuN);

    virtual double** EarthMatOscProbNUN(int type, double Erg);
    virtual double EarthMatOscProbNUN(int ij, int jk, int type, double Erg);
    virtual double** ProfiledMatOscProbNUN(int type, double Erg);
    virtual double ProfiledMatOscProbNUN(int ij, int jk, int type, double Erg);
    virtual double** VaccumOscProbNUN(int type, double Erg, double length);
    virtual double VaccumOscProbNUN(int ij, int jk, int type, double Erg, double length);
    virtual double** ConstMatOscProbNUN(int type, double Erg, double dens, double length);
    virtual double ConstMatOscProbNUN(int ij, int jk, int type, double Erg, double dens, double length);
  };

    #endif

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//





