/**
 *
 *  @file OscProbabilityHYBs.cxx
 *
 *  @brief It deals with Physics associated with Hybrid Physics.
 *
 *  It provides a unique opportunity to see how the oscillation probabilities get modified with the Hybrid Physics Lists. <br>
 *  It returns the oscillation probabilities for User's Profile density, Earth matter density.
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
    #include <OscComplxMatrixExp.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Oscillation Probability Calculation *//

    //* Profiled Matter Density Profile Oscillation Function *//
      double OscProbabilityHYB ::
      ProfiledMatOscProbHYB(int ij, int jk, int type, double Erg)
    {
      Alignment(ij, jk, dim, "Flavour", "ProfiledMatOscProbHYB");
      double** Prob;
      Prob = OscProbabilityHYB :: ProfiledMatOscProbHYB(type, Erg);
      return Prob[ij][jk];
      delete Prob;
    }
      double** OscProbabilityHYB ::
      ProfiledMatOscProbHYB(int type, double Erg)
   {
        _SYZ_VECTOR_CHECK_(Profile.size(), "SetGivenDensityProfile(double [][2], size_t )");
        if(type != 0)
     {
        OscProbabilityHYB :: May_I_Come_In();
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

        if (_MoA_ == 1) /// Hyb
      {
          U   = OscProbabilitySTD :: SetPMNS(Type);
          M   = OscProbabilitySTD :: SetMassHamiltonian(Erg);
          HF  = U*M*U.inverse();
          HF += OscProbabilityLIV :: SetGILIVMatrixA(Type, GeV);
          HF -= OscProbabilityLIV :: SetGILIVMatrixC(Type, Erg*GeV);

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Profile.size()); i++)
       {
          B       = OscProbabilityNSI :: SetNCNSIMatrix(Type, Profile[i].second);
          Ne      = OscProbabilitySTD :: GetInteractions(Type*Profile[i].second);
          HMat    = (HF+Ne+B)*((Profile[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
     }
        else if (_MoA_ == 2) /// FwI
      {
          U   = OscProbabilitySTD :: SetPMNS(Type);
          M   = OscProbabilitySTD :: SetMassHamiltonian(Erg);
          HF  = U*M*U.inverse();

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Profile.size()); i++)
       {
          Ne      = OscProbabilityHYB :: GetHbyInteractions(Type*Profile[i].second, _Ph_CC_, _Ph_NC_);
          HMat    = (HF+Ne)*((Profile[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
     }
        else if (_MoA_ == 3) /// Usr
      {
          HF  = OscProbabilityHYB :: GetUsrIntMatrix(Type, Erg, basis);

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Profile.size()); i++)
       {
          Ne      = OscProbabilitySTD :: GetInteractions(Type*Profile[i].second);
          HMat    = (HF+Ne)*((Profile[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
     }
        else if (_MoA_ == 4) /// Hyb, FwI
      {
          U   = OscProbabilitySTD :: SetPMNS(Type);
          M   = OscProbabilitySTD :: SetMassHamiltonian(Erg);
          HF  = U*M*U.inverse();
          HF += OscProbabilityLIV :: SetGILIVMatrixA(Type, GeV);
          HF -= OscProbabilityLIV :: SetGILIVMatrixC(Type, Erg*GeV);

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Profile.size()); i++)
       {
          B       = OscProbabilityNSI :: SetNCNSIMatrix(Type, Profile[i].second);
          Ne      = OscProbabilityHYB :: GetHbyInteractions(Type*Profile[i].second, _Ph_CC_, _Ph_NC_);
          HMat    = (HF+Ne+B)*((Profile[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
     }
        else if (_MoA_ == 5) /// FwI, Usr
      {
          HF  = OscProbabilityHYB :: GetUsrIntMatrix(Type, Erg, basis);

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Profile.size()); i++)
       {
          Ne      = OscProbabilityHYB :: GetHbyInteractions(Type*Profile[i].second, _Ph_CC_, _Ph_NC_);
          HMat    = (HF+Ne)*((Profile[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
     }
        else if (_MoA_ == 6) /// Hyb, Usr
      {
          HF  = OscProbabilityHYB :: GetUsrIntMatrix(Type, Erg, basis);
          HF += OscProbabilityLIV :: SetGILIVMatrixA(Type, GeV);
          HF -= OscProbabilityLIV :: SetGILIVMatrixC(Type, Erg*GeV);

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Profile.size()); i++)
       {
          B       = OscProbabilityNSI :: SetNCNSIMatrix(Type, Profile[i].second);
          Ne      = OscProbabilitySTD :: GetInteractions(Type*Profile[i].second);
          HMat    = (HF+Ne+B)*((Profile[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
     }
        else if (_MoA_ == 7) /// Hyb, FwI, Usr
      {
          HF  = OscProbabilityHYB :: GetUsrIntMatrix(Type, Erg, basis);
          HF += OscProbabilityLIV :: SetGILIVMatrixA(Type, GeV);
          HF -= OscProbabilityLIV :: SetGILIVMatrixC(Type, Erg*GeV);

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Profile.size()); i++)
       {
          B       = OscProbabilityNSI :: SetNCNSIMatrix(Type, Profile[i].second);
          Ne      = OscProbabilityHYB :: GetHbyInteractions(Type*Profile[i].second, _Ph_CC_, _Ph_NC_);
          HMat    = (HF+Ne+B)*((Profile[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
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
      double OscProbabilityHYB ::
      EarthMatOscProbHYB(int ij, int jk, int type, double Erg)
    {
      Alignment(ij, jk, dim, "Flavour", "EarthMatOscProbHYB");
      double** Prob;
      Prob = OscProbabilityHYB :: EarthMatOscProbHYB(type, Erg);
      return Prob[ij][jk];
      delete Prob;
    }
      double** OscProbabilityHYB ::
      EarthMatOscProbHYB(int type, double Erg)
   {
        if(type != 0)
     {
        OscProbabilityHYB :: May_I_Come_In();
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

        if (_MoA_ == 1) /// Hyb
      {
          U   = OscProbabilitySTD :: SetPMNS(Type);
          M   = OscProbabilitySTD :: SetMassHamiltonian(Erg);
          HF  = U*M*U.inverse();
          HF += OscProbabilityLIV :: SetGILIVMatrixA(Type, GeV);
          HF -= OscProbabilityLIV :: SetGILIVMatrixC(Type, Erg*GeV);

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Out_In.size()); i++)
       {
          B       = OscProbabilityNSI :: SetNCNSIMatrix(Type, Out_In[i].second);
          Ne      = OscProbabilitySTD :: GetInteractions(Type*Out_In[i].second);
          HMat    = (HF+Ne+B)*((Out_In[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
     }
        else if (_MoA_ == 2) /// FwI
      {
          U   = OscProbabilitySTD :: SetPMNS(Type);
          M   = OscProbabilitySTD :: SetMassHamiltonian(Erg);
          HF  = U*M*U.inverse();

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Out_In.size()); i++)
       {
          Ne      = OscProbabilityHYB :: GetHbyInteractions(Type*Out_In[i].second, _Ph_CC_, _Ph_NC_);
          HMat    = (HF+Ne)*((Out_In[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
     }
        else if (_MoA_ == 3) /// Usr
      {
          HF  = OscProbabilityHYB :: GetUsrIntMatrix(Type, Erg, basis);

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Out_In.size()); i++)
       {
          Ne      = OscProbabilitySTD :: GetInteractions(Type*Out_In[i].second);
          HMat    = (HF+Ne)*((Out_In[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
     }
        else if (_MoA_ == 4) /// Hyb, FwI
      {
          U   = OscProbabilitySTD :: SetPMNS(Type);
          M   = OscProbabilitySTD :: SetMassHamiltonian(Erg);
          HF  = U*M*U.inverse();
          HF += OscProbabilityLIV :: SetGILIVMatrixA(Type, GeV);
          HF -= OscProbabilityLIV :: SetGILIVMatrixC(Type, Erg*GeV);

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Out_In.size()); i++)
       {
          B       = OscProbabilityNSI :: SetNCNSIMatrix(Type, Out_In[i].second);
          Ne      = OscProbabilityHYB :: GetHbyInteractions(Type*Out_In[i].second, _Ph_CC_, _Ph_NC_);
          HMat    = (HF+Ne+B)*((Out_In[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
     }
        else if (_MoA_ == 5) /// FwI, Usr
      {
          HF  = OscProbabilityHYB :: GetUsrIntMatrix(Type, Erg, basis);

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Out_In.size()); i++)
       {
          Ne      = OscProbabilityHYB :: GetHbyInteractions(Type*Out_In[i].second, _Ph_CC_, _Ph_NC_);
          HMat    = (HF+Ne)*((Out_In[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
     }
        else if (_MoA_ == 6) /// Hyb, Usr
      {
          HF  = OscProbabilityHYB :: GetUsrIntMatrix(Type, Erg, basis);
          HF += OscProbabilityLIV :: SetGILIVMatrixA(Type, GeV);
          HF -= OscProbabilityLIV :: SetGILIVMatrixC(Type, Erg*GeV);

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Out_In.size()); i++)
       {
          B       = OscProbabilityNSI :: SetNCNSIMatrix(Type, Out_In[i].second);
          Ne      = OscProbabilitySTD :: GetInteractions(Type*Out_In[i].second);
          HMat    = (HF+Ne+B)*((Out_In[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
     }
        else if (_MoA_ == 7) /// Hyb, FwI, Usr
      {
          HF  = OscProbabilityHYB :: GetUsrIntMatrix(Type, Erg, basis);
          HF += OscProbabilityLIV :: SetGILIVMatrixA(Type, GeV);
          HF -= OscProbabilityLIV :: SetGILIVMatrixC(Type, Erg*GeV);

          MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
          for(int i = 0; i < int(Out_In.size()); i++)
       {
          B       = OscProbabilityNSI :: SetNCNSIMatrix(Type, Out_In[i].second);
          Ne      = OscProbabilityHYB :: GetHbyInteractions(Type*Out_In[i].second, _Ph_CC_, _Ph_NC_);
          HMat    = (HF+Ne+B)*((Out_In[i].first)*km2_eV);
          MatAmp *= MatExp(img*HMat);
       }
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






