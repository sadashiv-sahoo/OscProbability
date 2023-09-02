/**
 *
 *  @file OscGetSTDFunctions.cxx
 *
 *  @brief It deals with Physics associated with Standard Neutrino Oscillation Physics Functions.
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #include <OscProbability.h>
    #include <OscProbabilitySTD.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Getting Eigen-Functions *//
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic >
      GetModMatrix(const Eigen::Ref < const Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > > &A, int ck)
    {
      int dim = A.rows();
      Eigen::SelfAdjointEigenSolver < Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > > ces(dim);
      ces.compute(A);
      if (ck == 0) { return ces.eigenvectors(); }
      else         { return ces.eigenvalues().asDiagonal(); }
    }

      void OscProbabilitySTD ::
      Shall_I_Proceed_For_Modified_Calc()
    {
       unsigned int jk = 0; 
       for(unsigned int ij = 0; ij < NMix; ij ++)
     {
       if (MixAngs[ij].Ph != 0.0) {jk ++;}
     }
       if (jk == 0) {}
       else
     {
        std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nNon-Zero delta-CP(s) Found !!!\n\033[0m"
                  << std::endl;
        std::cout << "\033[0;96m:: Suggestion ::\nCurrently only Zero delta-CP(s) support(s).\033[0m\n" 
                  << std::endl;
        std::cout << "\033[0;96mPlease try with Mixing-Angles without non-zero delta-CP(s).\033[0m\n" 
                  << std::endl;
        exit(-1);
     }

    }

    //* Get the Oscillation Parameters (Standard) *//
      double** OscProbabilitySTD ::
      GetConstMatOscParms(int type, double dens, double Erg, const char *wh)
    {
        OscProbabilitySTD :: Shall_I_Proceed_For_Modified_Calc();
        if (type != 0 && dens >= 0.0)
      {
        const char *mxang = "MixMatrix";
        const char *msdff = "MassSqMatrix";
        int MxG = strcmp(mxang, wh);
        int MsD = strcmp(msdff, wh);

        double** OscParmRet;
        OscParmRet = new double*[dim];
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > U(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > M(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > N(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > H(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > X(dim,dim);

        int Type = int (type/(abs(type)));
        U  = OscProbabilitySTD :: SetPMNS(Type);
        M  = OscProbabilitySTD :: SetMassHamiltonian(Erg);
        N  = OscProbabilitySTD :: GetInteractions(Type*dens);
        M -= (1.0/dim)*M.trace()*(Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim));
        H  = M + U.inverse()*N*U;

          if (MxG == 0 && MsD != 0)      { X = U*GetModMatrix(H, 0); }
          else if (MxG != 0 && MsD == 0) { X = 2.0*Erg*GeV*GetModMatrix(H, 1); }
          else
        {
          std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nWrong Instruction Found !!!\n\033[0m"
                    << std::endl;
          std::cout << "\033[0;96m:: Suggestion ::\nCurrently only returns full modified Matrix.\033[0m\n"
                    << std::endl;
          std::cout << "\033[0;96mPlease try with either \"\033[0mMixMatrix\033[0;96m\" or \"\033[0mMassSqMatrix\033[0;96m\".\033[0m\n"
                    << std::endl;
          exit(-1);
        }
          for(unsigned int ij = 0; ij < dim; ij++)
        {
          OscParmRet[ij] = new double[dim];
          for(unsigned int jk = 0; jk < dim; jk++) { OscParmRet[ij][jk] = X(ij,jk).real(); }
        }
          return OscParmRet;
          delete OscParmRet;
      }
        else
      {
        std::cout << "\n\033[0;91mError :: Unphysical Value(s) (Negative or Zero) is(are) entered !!!\033[0m\n" 
                  << std::endl;
        exit(0);
      }
    }

    //* Get Neutrino Travel Distance Through Earth *//
      double OscProbabilitySTD ::
      GetTrav_Length(double cosZ)
    {
        if (cosZ >= -1.0 && cosZ <= 1.0)
      {
        double _Er_ = PREM[int(PREM.size())-1].first;
        double R    = _Er_ - Depth;
        return abs(sqrt(pow(R*cosZ,2) + (2.0*_Er_ + AtmosHeight - Depth)*(AtmosHeight + Depth)) - R*cosZ);
      }
        else
      {
        std::cout << "\n\033[0;91m:: Error Warning ::\nZenith Must Be [-1., 1]\033[0m\n"
                  << std::endl;
        exit(0);
      }
    }

    //* Get Neutrino Travel Distance Through Given Profile *//
      double OscProbabilitySTD ::
      GetTrav_Length(double profile[][2], size_t layers)
    {
      double TravL = 0.0;
      for(size_t i = 0; i < layers; i++) { TravL += profile[i][0]; }
      return TravL;
    }

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//





