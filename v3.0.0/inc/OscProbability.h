/**
 *
 *  @file  OscProbability.h
 *
 *  @brief OscProbability assigns the No. of flavors to perform calculations.
 *
 *  For performing the neutrino oscillation, we need atleast 2 neutrinos to act on. It can calulate the oscillation probabilities in a "n" No. of neutrino Paradigsm. <br>
 *  One needs to appropriately assign the oscillation parameters and controlling values via the predefined Struct.
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #ifndef _OscProbability_H
    #define _OscProbability_H

    #include <vector>
    #include <string>
    #include <complex>
    #include <iostream>

    #include <OscConstants.h>
    #include <eigen3/Eigen/Core>
    #include <eigen3/Eigen/Dense>
    #include <OscComplxMatrixExp.h>

    #include <OscProbabilitySTD.h>
    class OscProbabilitySTD;

    #include <OscProbabilityLIV.h>
    class OscProbabilityLIV;

    #include <OscProbabilityNSI.h>
    class OscProbabilityNSI;

    #include <OscProbabilityHYB.h>
    class OscProbabilityHYB;

    #include <OscProbabilityNUN.h>
    class OscProbabilityNUN;

    class OscProbability : virtual public OscProbabilitySTD,
                           virtual public OscProbabilityLIV,
                           virtual public OscProbabilityNSI,
                           virtual public OscProbabilityNUN,
                           virtual public OscProbabilityHYB
   {
      private:
      double X1; // dummy

      protected:
      double X2; // dummy

      public:
      OscProbability(const int n) { OscProbabilitySTD::SetDimension(n); }
      virtual ~OscProbability() {};
   };

    #endif

    #ifndef _MxParameters_H
    #define _MxParameters_H

    //* Flavour Assigning *//
    struct SetFlavoursPost { const char  *Fvs;
                             unsigned int Pos;
                           };

    //* Mixing Angles Variables *//
    struct SetMixAngParams { const double Th;
                             const double Ph;
                             unsigned int P1;
                             unsigned int P2;
                           };

    //* Mass Squared Difference *//
    struct SetMassSqParams { const double Mxx;};

    //* BSM Interaction Parameters *//
    struct SetBSMIntParams {
                             const std::complex < double > Bsm;
                             unsigned int F1;
                             unsigned int F2;
                           };

    #endif

    #ifndef _INLINE_CHECK_H
    #define _INLINE_CHECK_H

    //* Non-Negativity Check *//
    inline void _Non_Zero_Non_Neg_Check(const double X, const char *c)
  {
     const char *l = "Energy";
     int chk   = strcmp(l, c);
     if (X >= 0.)
   {
       if (chk == 0 && X > 0.){}
       else if (chk != 0 && X >= 0.) {}
       else
     {
       std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nUnphysical Value Entered !!!\n\033[0m\n"
                 << c << "\033[0;91m = \033[0m" << X << "\033[0;91m, and \033[0m"
                 << c << "\033[0;91m, Can't be < \"0\"\033[0m\n" 
                 << std::endl;
       std::cout << "\033[0;96mPlease Enter Non-Negative Value.\033[0m\n" 
                 << std::endl;
       exit(-1);
     }
   }
      else
    {
      std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nUnphysical Value Entered !!!\n\033[0m\n"
                << c << "\033[0;91m = \033[0m" << X << "\033[0;91m, and \033[0m"
                << c << "\033[0;91m, Can't be < \"0\"\033[0m\n" 
                << std::endl;
      std::cout << "\033[0;96mPlease Enter Non-Negative Value.\033[0m\n" 
                << std::endl;
      exit(0);
    }
  }

    //* Alignment Check *//
    inline void Alignment(const int ij, const int jk, const int dim, const char* c, const char* z)
  {
    const char *l = "Flavour";
    int chk = strcmp(l, c);
      if(chk == 0)
    {
        if (ij < 0 || ij >= dim || jk < 0 || jk >= dim)
      {
        std::cout << "\n\033[0;91m:: Error Warning ::\nIn " << z << ", " << c << " indices with (" << ij << ", " << jk << "), Can't Proceed !!!\033[0m\n"
                  << std::endl;
        std::cout << "\033[0;96mIndices must be <= \" \033[0m" << dim-1 << "\033[0;96m \".\033[0m\n" 
                  << std::endl;
        exit(-1);
      }
        else {}
    }
      else
    {
        if (ij < 1 || ij > dim || jk < 1 || jk > dim)
      {
        std::cout << "\n\033[0;91m:: Error Warning ::\nIn " << z << ", " << c << " indices with (" << ij << ", " << jk << "), Can't Proceed !!!\033[0m\n"
                  << std::endl;
        std::cout << "\033[0;96mIndices must be >= 1 and <= \" \033[0m" << dim << "\033[0;96m \".\033[0m\n" 
                  << std::endl;
        exit(-1);
      }
        else {}
    }
  }

    //* Number of Elements Check *//
    inline void BSMFlvchk(const size_t _Nb_, const int dim, const char* c)
  {
      const int Nt = dim*dim;
      const int Dl = dim*(dim-1)/2;
      const int Elements = Nt - Dl;
      if (int(_Nb_) <= Elements) {}
      else
    {
      std::cout << "\n\033[0;91m:: Error Warning ::\nNo. of " << c << " Flavour intercation indices with (\033[0m" << _Nb_ << "\033[0;91m), Can't Proceed !!!\033[0m\n"
                << std::endl;
      std::cout << "\033[0;96mTotal Indices must be <= \" \033[0m" << Elements << "\033[0;96m \".\033[0m\n" 
                << std::endl;
      exit(-1);
    }
  }

    //* Self-Adjoint Matrix Check *//
    inline void _BSM_SelfAdjoint_check(const int jk, const char* c)
  {
      if (jk > 0)
    {
      std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed !!!\nEntered element(s) don't make the \"" << c << "\" matrix \"Self-Adjoint\".\033[0m\n"
                << std::endl;
      std::cout << "\033[0;96mPlease follow the conventions used in this version.\033[0m\n" 
                << std::endl;
      exit(0);
    }
      else {}
  }

    //* Non-Vector Check *//
    inline void _SYZ_VECTOR_CHECK_(const size_t jk, const char* c)
  {
      if (jk > 0) {}
      else
    {
      std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed !!!\n"
                << std::endl;
      std::cout << "\033[0;96mPlease Enable \"" << c << "\".\033[0m\n" 
                << std::endl;
      exit(0);
    }
  }

    //* Physics Switch Check *//
    inline void May_I_Proceed(bool sts, const char* c)
  {
      if (sts == false)
    {
      std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed !!!\n"
                << std::endl;
      std::cout << "\033[0;96mPlease Switch ON the \"" << c << "(true)\".\033[0m\n"
                << std::endl;
      exit(0);
    }
      else {}
  }

    #endif

    #ifndef _GLIBCXX_SZ
    #define _GLIBCXX_SZ
    #undef SizeOf
    namespace std 
    _GLIBCXX_VISIBILITY(default)
  {
    _GLIBCXX_BEGIN_NAMESPACE_VERSION
    template <typename _Tp, std::size_t sz>
    std::size_t SizeOf(const _Tp(&)[sz]) { return sz; }
  }
    #endif

    #ifndef _ModMatrix_H
    #define _ModMatrix_H

    Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic >
    GetModMatrix(const Eigen::Ref < const Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > > &A, int ck);

    #endif








