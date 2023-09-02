/**
 *
 *  @file OscProbabilityLIV.h
 *
 *  @brief It deals with Physics associated with Spontaneous Lorentz Invariance Violations (LIV).
 *
 *  It provides a unique opportunity to see how the oscillation probabilities get modified with the LIV Physics Lists. <br>
 *  It can handle both Gauge Conserving and Violating LIV Physics.
 *  This programme is motivivated by a paper <a href="https://arxiv.org/pdf/hep-ph/0309025.pdf"> "Lorentz and CPT violation in neutrinos" </a>
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #ifndef _OscProbabilityLIV_H
    #define _OscProbabilityLIV_H

    #include <cmath>
    #include <OscProbability.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Class Definations *//

    class OscProbabilityLIV : virtual OscProbabilitySTD
  {
    private:
    //* ....oooOOOOO Gauge Invariant LIV Settings OOOOOooo.... *//
    //* GILIV Control Parameter factors *//
    double aL = 0.0;
    double cL = 0.0;

    //* GILIV Parameters *//
    std::vector < SetBSMIntParams > GILIVCPTODD;
    std::vector < SetBSMIntParams > GILIVCPTEVEN;

    //* ....oooOOOOO Gauge Violating LIV Settings OOOOOooo.... *//
    //* GVLIV Control Parameter factors *//
    double _h_ = 0.0;
    double _g_ = 0.0;
    std::complex < double > _eps_;

    //* GVLIV Parameters *//
    std::vector < SetBSMIntParams > GVLIVCPTODD;
    std::vector < SetBSMIntParams > GVLIVCPTEVEN;

    public:
    virtual ~OscProbabilityLIV() {};
    virtual void SetLIVSpTControls(int _mu_, int _nu_);
    virtual void SetGILIVaX(SetBSMIntParams TaX[], size_t _Tav_);
    virtual void SetGILIVcXX(SetBSMIntParams TcXX[], size_t _Tcv_);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > SetGILIVMatrixA(int type, double _units_);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > SetGILIVMatrixC(int type, double _units_);

    virtual void SetGVLIVPolzation(std::complex < double > eps);
    virtual void SetGILIVhXX(SetBSMIntParams ThXX[], size_t _Thv_);
    virtual void SetGILIVgXXX(SetBSMIntParams TgXXX[], size_t _Tgv_);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > SetGVLIVMatrixG();
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > SetGVLIVMatrixH();
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > SetGVLIVMatrixC();
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > SetGVLIVHamiltonian(int _Init_, double Erg);

    //* Oscillation LIV Probability Functions *//
    virtual double** EarthMatOscProbGILIV(int type, double Erg);
    virtual double EarthMatOscProbGILIV(int ij, int jk, int type, double Erg);                             /// Earth Matter Oscilation Probabilty;
    virtual double** ProfiledMatOscProbGILIV(int type, double Erg);
    virtual double ProfiledMatOscProbGILIV(int ij, int jk, int type, double Erg);                          /// Profiled Matter Oscilation Probabilty;
    virtual double** VaccumOscProbGILIV(int type, double Erg, double length);
    virtual double VaccumOscProbGILIV(int ij, int jk, int type, double Erg, double length);                /// Vaccum Oscillation Probability;
    virtual double** ConstMatOscProbGILIV(int type, double Erg, double dens, double length);
    virtual double ConstMatOscProbGILIV(int ij, int jk, int type, double Erg, double dens, double length); /// Const-Dens-Mat Oscillation Probability;
    virtual double** VaccumOscProbGVLIV(int _Init_, double Erg, double length);
    virtual double VaccumOscProbGVLIV(int ij, int jk, int _Init_, double Erg, double length);              /// Vaccum Gauge Violating Oscillation Probability;

    //* Get Oscillation Parameters with LIV *//
    virtual double** GetVacLIVOscParms(int type, double Erg, const char *wh);                              /// Vaccum Modified Matrices;
    virtual double** GetConstMatLIVOscParms(int type, double dens, double Erg, const char *wh);            /// Constant Matter Modified Matrices;
  };

    inline void _LIV_Control_Check_(double _ck_)
  {
      if (_ck_ != 0.) {}
      else
    {
      std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nSpace-Time Indices are not Found!!!\n\033[0m"
                << std::endl;
      std::cout << "\033[0;96mPlease, Make Sure of Enabling \" SetLIVSpTControls(int,int) \".\033[0m\n" 
                << std::endl;
      exit(0);
    }
  }

    #endif

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//






