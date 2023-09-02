/**
 *
 *  @file OscProbabilityNSIs.cxx
 *
 *  @brief It deals with Physics associated with Neutral Current Forward Elastic Non-Standard fermionic Interactions.
 *
 *  It provides a unique opportunity to see how the oscillation probabilities get modified with the NC-NSI Physics Lists. <br>
 *  It returns the oscillation probabilities for Constant matter density, User's Profile density, Earth matter density.
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
    #include <OscComplxMatrixExp.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Oscillation Probability Calculation *//

  //* Const Matter Density Oscillation Function *//
      double OscProbabilityNSI ::
      ConstMatOscProbNCNSI(int ij, int jk, int type, double Erg, double dens, double length)
    {
      Alignment(ij, jk, dim, "Flavour", "ConstMatOscProbNCNSI");
      double** Prob;
      Prob = OscProbabilityNSI :: ConstMatOscProbNCNSI(type, Erg, dens, length);
      return Prob[ij][jk];
      delete Prob;
    }
      double** OscProbabilityNSI ::
      ConstMatOscProbNCNSI(int type, double Erg, double dens, double length)
   {
        if(type != 0 && dens >= 0.00 && length >= 0.00)
     {
        double** P;
        P = new double*[dim];
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > U(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > M(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > HF(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > Ne(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > HMat(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > MatAmp(dim,dim);
        int Type = int (type/(abs(type)));

        U   = OscProbabilitySTD :: SetPMNS(Type);
        M   = OscProbabilitySTD :: SetMassHamiltonian(Erg);
        HF  = U*M*U.inverse();

        Ne     = OscProbabilitySTD :: GetInteractions(Type*dens);
        Ne    += OscProbabilityNSI :: SetNCNSIMatrix(Type, dens);
        HMat   = (HF+Ne)*(length*km2_eV);
        MatAmp = MatExp(img*HMat);
        for(unsigned int ij = 0; ij < dim; ij++)
      {
          P[ij] = new double[dim];
          for(unsigned int jk = 0; jk < dim; jk++)
        {
          P[ij][jk] = pow(abs(MatAmp(ij,jk)),2);
        }
      }
        return P;
        delete P;
     }
        else
      {
        std::cout << "\n\033[0;91mError :: Unphysical Value(s) (Negative or Zero) is(are) entered !!!\033[0m\n" 
                  << std::endl;
        exit(0);
      }
   }

  //* Profiled Matter Density Profile Oscillation Function *//
      double OscProbabilityNSI ::
      ProfiledMatOscProbNCNSI(int ij, int jk, int type, double Erg)
    {
      Alignment(ij, jk, dim, "Flavour", "ProfiledMatOscProbNCNSI");
      double** Prob;
      Prob = OscProbabilityNSI :: ProfiledMatOscProbNCNSI(type, Erg);
      return Prob[ij][jk];
      delete Prob;
    }
      double** OscProbabilityNSI ::
      ProfiledMatOscProbNCNSI(int type, double Erg)
   {
        _SYZ_VECTOR_CHECK_(Profile.size(), "SetGivenDensityProfile(double [][2], size_t )");
        if(type != 0)
     {
        double** P;
        P = new double*[dim];
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > U(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > M(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > B(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > HF(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > Ne(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > HMat(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > MatAmp(dim,dim);

        int Type = int (type/(abs(type)));
        U   = OscProbabilitySTD :: SetPMNS(Type);
        M   = OscProbabilitySTD :: SetMassHamiltonian(Erg);
        HF  = U*M*U.inverse();

        MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
        for(int i = 0; i < int(Profile.size()); i++)
     {
        B       = OscProbabilityNSI :: SetNCNSIMatrix(Type, Profile[i].second);
        Ne      = OscProbabilitySTD :: GetInteractions(Type*Profile[i].second);
        HMat    = (HF+Ne+B)*((Profile[i].first)*km2_eV);
        MatAmp *= MatExp(img*HMat);
     }
        for(unsigned int ij = 0; ij < dim; ij++)
      {
          P[ij] = new double[dim];
          for(unsigned int jk = 0; jk < dim; jk++)
        {
          P[ij][jk] = pow(abs(MatAmp(ij,jk)),2);
        }
      }
        return P;
        delete P;
     }
        else
      {
        std::cout << "\n\033[0;91mError :: Unphysical Value(s) (Negative or Zero) is(are) entered !!!\033[0m\n" 
                  << std::endl;
        exit(0);
      }
   }

  //* Earth Matter Density Profile Oscillation Function *//
      double OscProbabilityNSI ::
      EarthMatOscProbNCNSI(int ij, int jk, int type, double Erg)
    {
      Alignment(ij, jk, dim, "Flavour", "EarthMatOscProbNCNSI");
      double** Prob;
      Prob = OscProbabilityNSI :: EarthMatOscProbNCNSI(type, Erg);
      return Prob[ij][jk];
      delete Prob;
    }
      double** OscProbabilityNSI ::
      EarthMatOscProbNCNSI(int type, double Erg)
   {
        if(type != 0)
     {
        double** P;
        P = new double*[dim];
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > U(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > M(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > B(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > HF(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > Ne(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > HMat(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > MatAmp(dim,dim);

        int Type = int (type/(abs(type)));
        U   = OscProbabilitySTD :: SetPMNS(Type);
        M   = OscProbabilitySTD :: SetMassHamiltonian(Erg);
        HF  = U*M*U.inverse();

        MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
        for(int i = 0; i < int(Out_In.size()); i++)
     {
        B       = OscProbabilityNSI :: SetNCNSIMatrix(Type, Out_In[i].second);
        Ne      = OscProbabilitySTD :: GetInteractions(Type*Out_In[i].second);
        HMat    = (HF+Ne+B)*((Out_In[i].first)*km2_eV);
        MatAmp *= MatExp(img*HMat);
     }
        for(unsigned int ij = 0; ij < dim; ij++)
      {
          P[ij] = new double[dim];
          for(unsigned int jk = 0; jk < dim; jk++)
        {
          P[ij][jk] = pow(abs(MatAmp(ij,jk)),2);
        }
      }
        return P;
        delete P;
     }
        else
      {
        std::cout << "\n\033[0;91mError :: Unphysical Value(s) (Negative or Zero) is(are) entered !!!\033[0m\n" 
                  << std::endl;
        exit(0);
      }
   }

  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//






