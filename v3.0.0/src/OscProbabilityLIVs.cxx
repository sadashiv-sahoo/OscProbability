/**
 *
 *  @file OscProbabilityLIVs.cxx
 *
 *  @brief It deals with Physics associated with Spontaneous Lorentz Invariance Violations (LIV).
 *
 *  It provides a unique opportunity to see how the oscillation probabilities get modified with the LIV Physics Lists. <br>
 *  It returns the oscillation probabilities for Gauge Conserving Vacuum, Constant matter density, User's Profile density, Earth matter density, and a Vacuum Gauge Violations. 
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #include <OscProbability.h>
    #include <OscProbabilityLIV.h>
    #include <OscComplxMatrixExp.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Oscillation Probability Calculation *//

    //* Vaccum Oscillation Function *//
      double OscProbabilityLIV :: 
      VaccumOscProbGILIV(int ij, int jk, int type, double Erg, double length)
    {
      Alignment(ij, jk, dim, "Flavour", "VaccumOscProbGILIV");
      double** Prob;
      Prob = OscProbabilityLIV :: VaccumOscProbGILIV(type, Erg, length);
      return Prob[ij][jk];
      delete Prob;
    }

      double** OscProbabilityLIV ::
      VaccumOscProbGILIV(int type, double Erg, double length)
    {
        if(type != 0 && length >= 0.00)
     {
        double** P;
        P = new double*[dim];
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > U(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > M(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > HF(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > Hvac(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > VacAmp(dim,dim);
        int Type = int (type/(abs(type)));

        U   = OscProbabilitySTD :: SetPMNS(Type);
        M   = OscProbabilitySTD :: SetMassHamiltonian(Erg);
        HF  = U*M*U.inverse();
        HF += OscProbabilityLIV :: SetGILIVMatrixA(Type, GeV);
        HF -= OscProbabilityLIV :: SetGILIVMatrixC(Type, Erg*GeV);

        Hvac   = HF*length*km2_eV;
        VacAmp = MatExp(img*Hvac);

        for(unsigned int ij = 0; ij < dim; ij++)
      {
          P[ij] = new double[dim];
          for(unsigned int jk = 0; jk < dim; jk++)
        {
          P[ij][jk] = pow(abs(VacAmp(ij,jk)),2);
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

    //* Const Matter Density Oscillation Function *//
      double OscProbabilityLIV ::
      ConstMatOscProbGILIV(int ij, int jk, int type, double Erg, double dens, double length)
    {
      Alignment(ij, jk, dim, "Flavour", "ConstMatOscProbGILIV");
      double** Prob;
      Prob = OscProbabilityLIV :: ConstMatOscProbGILIV(type, Erg, dens, length);
      return Prob[ij][jk];
      delete Prob;
    }
      double** OscProbabilityLIV ::
      ConstMatOscProbGILIV(int type, double Erg, double dens, double length)
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
        HF += OscProbabilityLIV :: SetGILIVMatrixA(Type, GeV);
        HF -= OscProbabilityLIV :: SetGILIVMatrixC(Type, Erg*GeV);

        Ne     = OscProbabilitySTD :: GetInteractions(Type*dens);
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
      double OscProbabilityLIV ::
      ProfiledMatOscProbGILIV(int ij, int jk, int type, double Erg)
    {
      Alignment(ij, jk, dim, "Flavour", "ProfiledMatOscProbGILIV");
      double** Prob;
      Prob = OscProbabilityLIV :: ProfiledMatOscProbGILIV(type, Erg);
      return Prob[ij][jk];
      delete Prob;
    }
      double** OscProbabilityLIV ::
      ProfiledMatOscProbGILIV(int type, double Erg)
   {
        _SYZ_VECTOR_CHECK_(Profile.size(), "SetGivenDensityProfile(double [][2], size_t )");
        if(type != 0)
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
        HF += OscProbabilityLIV :: SetGILIVMatrixA(Type, GeV);
        HF -= OscProbabilityLIV :: SetGILIVMatrixC(Type, Erg*GeV);

        MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
        for(int i = 0; i < int(Profile.size()); i++)
     {
        Ne      = OscProbabilitySTD :: GetInteractions(Type*Profile[i].second);
        HMat    = (HF+Ne)*((Profile[i].first)*km2_eV);
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
      double OscProbabilityLIV ::
      EarthMatOscProbGILIV(int ij, int jk, int type, double Erg)
    {
      Alignment(ij, jk, dim, "Flavour", "EarthMatOscProbGILIV");
      double** Prob;
      Prob = OscProbabilityLIV :: EarthMatOscProbGILIV(type, Erg);
      return Prob[ij][jk];
      delete Prob;
    }
      double** OscProbabilityLIV ::
      EarthMatOscProbGILIV(int type, double Erg)
   {
        if(type != 0)
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
        HF += OscProbabilityLIV :: SetGILIVMatrixA(Type, GeV);
        HF -= OscProbabilityLIV :: SetGILIVMatrixC(Type, Erg*GeV);

        MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
        for(int i = 0; i < int(Out_In.size()); i++)
     {
        Ne      = OscProbabilitySTD :: GetInteractions(Type*Out_In[i].second);
        HMat    = (HF+Ne)*((Out_In[i].first)*km2_eV);
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

      double OscProbabilityLIV :: 
      VaccumOscProbGVLIV(int ij, int jk, int _init_, double Erg, double length)
    {
      Alignment(ij, jk, dim, "Flavour", "VaccumOscProbGVLIV");
      double** Prob;
      Prob = OscProbabilityLIV :: VaccumOscProbGVLIV(_init_, Erg, length);
      return Prob[ij][jk];
      delete Prob;
    }

      double** OscProbabilityLIV ::
      VaccumOscProbGVLIV(int _init_, double Erg, double length)
    {
        if(_init_ != 0 && length >= 0.00 && Erg > 0.0)
     {
        double** P;
        P = new double*[dim];
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > HF(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > VacAmp(dim,dim);

        HF = OscProbabilityLIV :: SetGVLIVHamiltonian(_init_, Erg);
        if(_init_ > 0) {-img*sqrt(2.0)*HF;}
        else           {img*sqrt(2.0)*HF.conjugate();}
        VacAmp = MatExp(img*HF*length*km2_eV);

        for(unsigned int ij = 0; ij < dim; ij++)
      {
          P[ij] = new double[dim];
          for(unsigned int jk = 0; jk < dim; jk++)
        {
          P[ij][jk] = pow(abs(VacAmp(ij,jk)),2);
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






