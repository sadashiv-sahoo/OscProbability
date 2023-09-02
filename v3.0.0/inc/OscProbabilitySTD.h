/**
 *
 *  @file  OscProbabilitySTD.h
 *
 *  @brief OscProbabilitySTD is heart of this whole progamme. Here, "STD" stands for Standard.
 *
 *  For performing the neutrino oscillation, It assigns the values of mass-squared splittings, Mixing Angles, ... etc.
 *  One needs to appropriately assign the oscillation parameters and controlling values via the predefined Struct.
 *  This can return the oscillation probabilities for Vaccum, Constant matter density, User's Profile density and Earth matter density.
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #ifndef _OscProbabilitySTD_H
    #define _OscProbabilitySTD_H

    #include <cmath>
    #include <algorithm>
    #include <OscProbability.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Class Definations *//

    class OscProbabilitySTD
  {
    protected:
    unsigned int dim;
    unsigned int NMix;
    unsigned int NMass;

    //* SM Interaction *//
    double _ne_ = 0.5;
    double _np_ = 0.5;
    double _nn_ = 0.5;

    //* Earth Variables *//                               /// All Scales are in KMs;
    double Depth;
    double AtmosHeight;
    double EarthRadius;

    //* Earth Density Variables *//
    std::vector < std::pair < double, double > > PREM;    /// PREM Global Fixed   Data Store;
    std::vector < std::pair < double, double > > Out_In;  /// PREM Local  Dynamic Data Store;
    std::vector < std::pair < double, double > > Profile; /// User's      Fixed   Data Store;

    //* Oscillation Parameters *//
    std::vector < SetFlavoursPost > FlvPost;
    std::vector < SetMixAngParams > MixAngs;
    std::vector < SetMassSqParams > MSqDiff;

    public:
    virtual ~OscProbabilitySTD() {};
    void SetEarthDensityProfile();                                      /// Loads Earth's PREM Profile for Atmospheric Expt.;
    void SetGivenDensityProfile(double profile[][2], size_t layers);    /// Loads Fixed Baseline Expt.;

    //* Setting Flavour Dimensions *//
    void SetDimension(const int n);
    virtual void _No_of_Mix_Angles_Check(size_t nt);
    virtual void _No_of_MsSq_Diffs_Check(size_t Mt);
    virtual void SetMixAngles(SetMixAngParams MxAg[], size_t At);
    virtual void SetMassSqDiff(SetMassSqParams MSq[], size_t Mt);
    virtual void SetOscFlavours(SetFlavoursPost Fl[], size_t Ft);
    virtual void SetMediumPlrFrac(double Ye, double Yp, double Yn);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > SetPMNS(int type);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > GetInteractions(double rho);
    virtual Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > SetMassHamiltonian(double Erg);

    //* EarthMatOscProb Set Functions *//
    virtual void SetZenith(double CosZ);                                          /// Includes Required Zenith for Analysis;
    virtual void SetDetDepth(double depth);                                       /// Includes Detector Depth;
    virtual void SetEarthAtmos(double AtmHeight);                                 /// Includes Atmospheric Production Height;

    //* Oscillation Probability Functions *//
    virtual double** EarthMatOscProb(int type, double Erg);
    virtual double EarthMatOscProb(int ij, int jk, int type, double Erg);                             /// Earth Matter Oscilation Probabilty;
    virtual double** ProfiledMatOscProb(int type, double Erg);
    virtual double ProfiledMatOscProb(int ij, int jk, int type, double Erg);                          /// Profiled Matter Oscilation Probabilty;
    virtual double** VaccumOscProb(int type, double Erg, double length);
    virtual double VaccumOscProb(int ij, int jk, int type, double Erg, double length);                /// Vaccum Oscillation Probability;
    virtual double** ConstMatOscProb(int type, double Erg, double dens, double length);
    virtual double ConstMatOscProb(int ij, int jk, int type, double Erg, double dens, double length); /// Const-Dens-Mat Oscillation Probability;

    //* Get Standard Functions *//
    virtual double GetTrav_Length(double cosZ);                                                       /// Atmospheric Baseline(s);
    virtual double GetTrav_Length(double profile[][2], size_t layers);                                /// Fixed Baseline

    //* Get Oscillation Parameters *//
    virtual void Shall_I_Proceed_For_Modified_Calc();
    virtual double** GetConstMatOscParms(int type, double dens, double Erg, const char *wh);          /// Constant Matter Modified Matrices;
  };

    //* SM Neutrino (electron, muon, tauon) *//
    inline void _Std_Flavour_Type_Check(const char *ch)
  {
    const char *e = "e";
    const char *m = "mu";
    const char *t = "tau";
    if (int(strcmp(e, ch)) == 0 || int(strcmp(m, ch)) == 0 || int(strcmp(t, ch)) == 0) { }
      else 
    {
      std::cout << "\n\033[0;91m:: Error Warning ::\nFlavour with name (\033[0m" << ch << "\033[0;91m), Can't Proceed !!!\033[0m\n"
                << std::endl;
      std::cout << "\033[0;96mSM Neutrinos must be treated as (\"e\", \"mu\", \"tau\").\033[0m\n" 
                << std::endl;
      exit(0);
    }
  }

    //* Propagation Length in Earth Profile *//
    inline double chordL(double R, double d)
  {
    if(R >= d) { return abs(sqrt(pow(R,2)-pow(d,2))); }
    else       { return 0; }
  }

    //* Concentric Propagations Choice *//
    inline double Travers_Dist(double l1, double l2) { return abs(l1 -l2); }

    #endif

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//






