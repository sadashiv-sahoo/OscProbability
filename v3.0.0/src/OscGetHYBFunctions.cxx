/**
 *
 *  @file OscGetHYBFunctions.cxx
 *
 *  @brief It deals with Physics associated with Hybrid Physics Functions.
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #include <OscProbability.h>
    #include <OscProbabilityHYB.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityHYB ::
      GetUsrIntMatrix(int _md_)
    {
      _SYZ_VECTOR_CHECK_(_GBLUSER_INTERACTIONS_.size(), "SetGBLUserInteractions(SetBSMIntParams , size_t , const char )");
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > H(dim,dim);
      H = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Zero(dim,dim);
        if (_md_ == 999)
      {
          for (size_t ij = 0; ij < _GBLUSER_INTERACTIONS_.size(); ij ++)
        {
            if   (_GBLUSER_INTERACTIONS_[ij].F1 == _GBLUSER_INTERACTIONS_[ij].F2) { H(_GBLUSER_INTERACTIONS_[ij].F1, _GBLUSER_INTERACTIONS_[ij].F2) = _GBLUSER_INTERACTIONS_[ij].Bsm; }
            else
          {
            std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nWrong Valued Enetered !!!\n\033[0m"
                      << std::endl;
            std::cout << "\033[0;96mFor Mass Basis, Please Enter with the \"\033[0;92mDiagonal\033[0;96m\" "
                         "Elements in the \"SetGBLUserInteractions(SetBSMIntParams , size_t , const char )\".\033[0m\n"
                      << std::endl;
            exit(0);
          }
        }
      }
        else if (_md_ == 111) { for (size_t ij = 0; ij < _GBLUSER_INTERACTIONS_.size(); ij ++) { H(_GBLUSER_INTERACTIONS_[ij].F1, _GBLUSER_INTERACTIONS_[ij].F2) = _GBLUSER_INTERACTIONS_[ij].Bsm; } }
        else
      {
        std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nNo USER Valued Enetered !!!\n\033[0m"
                  << std::endl;
        std::cout << "\033[0;96mPlease Enter Elements in the \"SetGBLUserInteractions(SetBSMIntParams , size_t , const char )\".\033[0m\n"
                  << std::endl;
        exit(-1);
      }
        return H;
    }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityHYB ::
      GetUsrIntMatrix(int type, double Erg, int _md_)
    {
      _SYZ_VECTOR_CHECK_(_GBLUSER_INTERACTIONS_.size(), "SetGBLUserInteractions(SetBSMIntParams , size_t , const char )");
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > U(dim,dim);
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > M(dim,dim);
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > H(dim,dim);
      U = OscProbabilitySTD :: SetPMNS(type);
      M = OscProbabilitySTD :: SetMassHamiltonian(Erg);
      H = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Zero(dim,dim);
        if (_md_ == 999)
      {
          for (size_t ij = 0; ij < _GBLUSER_INTERACTIONS_.size(); ij ++)
        { 
            if   (_GBLUSER_INTERACTIONS_[ij].F1 == _GBLUSER_INTERACTIONS_[ij].F2) { H(_GBLUSER_INTERACTIONS_[ij].F1-1, _GBLUSER_INTERACTIONS_[ij].F2-1) = _GBLUSER_INTERACTIONS_[ij].Bsm; }
            else
          {
            std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nWrong Valued Enetered !!!\n\033[0m"
                      << std::endl;
            std::cout << "\033[0;96mFor Mass Basis, Please Enter with the Diagonal Elements in the \"SetGBLUserInteractions(SetBSMIntParams , size_t , const char )\".\033[0m\n"
                      << std::endl;
            exit(0);
          }
        }
        return U*(M+H)*U.inverse();
      }
        else if (_md_ == 111)
      {
        for (size_t ij = 0; ij < _GBLUSER_INTERACTIONS_.size(); ij ++) { H(_GBLUSER_INTERACTIONS_[ij].F1, _GBLUSER_INTERACTIONS_[ij].F2) = _GBLUSER_INTERACTIONS_[ij].Bsm; }
        return U*M*U.inverse() + H;
      }
        else
      {
        std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nNo USER Valued Enetered !!!\n\033[0m"
                  << std::endl;
        std::cout << "\033[0;96mPlease Enter Elements in the \"SetGBLUserInteractions(SetBSMIntParams , size_t , const char )\".\033[0m\n"
                  << std::endl;
        exit(-1);
      }
    }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityHYB ::
      GetHbyInteractions(double rho, double PhCC, double PhNC)
   {
       const char *e = "e";
       const char *m = "mu";
       const char *t = "tau";
       May_I_Proceed(_InElst_Status_, "SetFwdInElasticSts");
       if(FlvPost.size() != 0)
     {
       typedef std::complex < double > r_;
       const std::complex < double > CC = sqrt2GF*rho*_ne_*r_(cos(PhCC), sin(PhCC));
       const std::complex < double > NC = GFbysqrt2*rho*(Cv_by_Ca*(_np_ - _ne_) - _nn_)*r_(cos(PhNC), sin(PhNC));

       Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > Ne(dim, dim);
       Ne = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Zero(dim,dim);

         for(size_t ij = 0; ij < FlvPost.size(); ij ++)
       {
              if(int(strcmp(e, FlvPost[ij].Fvs)) == 0) { Ne(FlvPost[ij].Pos, FlvPost[ij].Pos) = NC + CC; }
         else if(int(strcmp(m, FlvPost[ij].Fvs)) == 0) { Ne(FlvPost[ij].Pos, FlvPost[ij].Pos) = NC;      }
         else if(int(strcmp(t, FlvPost[ij].Fvs)) == 0) { Ne(FlvPost[ij].Pos, FlvPost[ij].Pos) = NC;      }
       }
         return Ne;
     }
       else
     {
       std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed !!!\n\033[0m"
                 << std::endl;
       std::cout << "\033[0;96mElse, Make Sure of Enabling \" SetOscFlavours(#, #) \".\033[0m\n" 
                 << std::endl;
       exit(0);
     }
   }

    //* Get the Oscillation Parameters (NSI) *//
      double** OscProbabilityHYB ::
      GetConstMatHYBOscParms(int type, double dens, double Erg, const char *wh)
    {
        OscProbabilitySTD :: Shall_I_Proceed_For_Modified_Calc();
        if (type != 0 && dens >= 0.0)
      {
        OscProbabilityHYB :: May_I_Come_In();
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

        if (_MoA_ == 1)      /// Hyb
      {
        U  = OscProbabilitySTD :: SetPMNS(Type);
        M  = OscProbabilitySTD :: SetMassHamiltonian(Erg);
        N  = OscProbabilitySTD :: GetInteractions(Type*dens);
        N += OscProbabilityLIV :: SetGILIVMatrixA(Type, GeV);
        N -= OscProbabilityLIV :: SetGILIVMatrixC(Type, Erg*GeV);
        N += OscProbabilityNSI :: SetNCNSIMatrix(Type, dens);
      }
        else if (_MoA_ == 2) /// FwI
      {
        U = OscProbabilitySTD :: SetPMNS(Type);
        M = OscProbabilitySTD :: SetMassHamiltonian(Erg);
        N = OscProbabilityHYB :: GetHbyInteractions(Type*dens, _Ph_CC_, _Ph_NC_);
      }
        else if (_MoA_ == 3) /// Usr
      {
        U = OscProbabilitySTD :: SetPMNS(Type);
        M = OscProbabilitySTD :: SetMassHamiltonian(Erg);
        N = OscProbabilitySTD :: GetInteractions(Type*dens);
          if      (basis == 999) { M += OscProbabilityHYB :: GetUsrIntMatrix(basis);}
          else if (basis == 111) { N += OscProbabilityHYB :: GetUsrIntMatrix(basis);}
          else
        {
          std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nNo Valued Enetered !!!\n\033[0m"
                    << std::endl;
          std::cout << "\033[0;96mPlease Enter Elements in the \"SetGBLUserInteractions(SetBSMIntParams , size_t , const char )\".\033[0m\n"
                    << std::endl;
          exit(-1);
        }
      }
        else
      {
        std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nWrong Instruction Found !!!\n\033[0m"
                  << std::endl;
        std::cout << "\033[0;96m:: Suggestion ::\nCurrently only returns full modified Matrix for a Particular Hybrid mode.\033[0m\n"
                  << std::endl;
        std::cout << "\033[0;96mPlease try Hybrid with either \"\033[0mLIV + NSI\033[0;96m\" or \"\033[0mForward Inelastic Phases\033[0;96m\" or \"\033[0mUser Custom Interactions\033[0;96m\"\033[0m\n"
                  << std::endl;
        exit(-1);
      }

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

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//





