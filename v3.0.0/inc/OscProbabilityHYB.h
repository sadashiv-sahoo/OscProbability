/**
 *
 *  @file OscProbabilityHYB.h
 *
 *  @brief It deals with Hybrid functions of new Physics
 *
 *  It provides a unique opportunity to see how the oscillation probabilities get modified with the Hybrid Physics Lists. <br>
 *  It can handle, LIV + NSI, User's defined new Field, Modification in standard forward elastic scatterings.
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #ifndef _OscProbabilityHYB_H
    #define _OscProbabilityHYB_H

    #include <OscProbability.h>
    #include <OscProbabilityHYB.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Class Definations *//

    class OscProbabilityHYB : virtual OscProbabilitySTD,
                              virtual OscProbabilityLIV,
                              virtual OscProbabilityNSI
  {
    private:
    bool _HybPhy_Status_ = false;
    bool _InElst_Status_ = false;
    bool _UsrInt_Status_ = false;
    int  _MoA_ = 0;
    int  basis = 0;
    double _Ph_CC_ = 0.0;
    double _Ph_NC_ = 0.0;
    std::vector < SetBSMIntParams > _GBLUSER_INTERACTIONS_;

    public:
    virtual ~OscProbabilityHYB() {};

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 
    //* Hybrid Physics Settings *//
    virtual void May_I_Come_In();
    virtual void SetHybridPhyStatus(bool _hsts_);
    virtual void SetFwdInElasticSts(bool _fsts_);
    virtual void SetGBLUsersIntrStatus(bool _usrsts_);
    virtual void SetMatterInElasticPhases(double _Pcc_, double _Pnc_);
    virtual void SetGBLUserInteractions(SetBSMIntParams usr[], size_t _usr_, const char* c);

    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > GetUsrIntMatrix(int _md_);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > GetUsrIntMatrix(int type, double Erg, int _md_);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > GetHbyInteractions(double rho, double PhCC, double PhNC);
    virtual double** GetConstMatHYBOscParms(int type, double dens, double Erg, const char *wh);

    virtual double** EarthMatOscProbHYB(int type, double Erg);
    virtual double EarthMatOscProbHYB(int ij, int jk, int type, double Erg);
    virtual double** ProfiledMatOscProbHYB(int type, double Erg);
    virtual double ProfiledMatOscProbHYB(int ij, int jk, int type, double Erg);
  };

    #endif

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//





