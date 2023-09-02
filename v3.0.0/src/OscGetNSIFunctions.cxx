/**
 *
 *  @file OscGetNSIFunctions.cxx
 *
 *  @brief It deals with Physics associated with NC-NSI related Functions.
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #include <OscProbability.h>
    #include <OscProbabilityNSI.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//

    //* Get the Oscillation Parameters (NSI) *//
      double** OscProbabilityNSI ::
      GetConstMatNSIOscParms(int type, double dens, double Erg, const char *wh)
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
        N += OscProbabilityNSI :: SetNCNSIMatrix(Type, dens);
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






