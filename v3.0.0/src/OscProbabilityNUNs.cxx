/**
 *
 *  @file OscProbabilityNUNs.cxx
 *
 *  @brief It deals with Physics associated with Non Unitarity Mixing Physics.
 *
 *  It provides a unique opportunity to see how the oscillation probabilities get modified with the Non-Unitarity Physics Lists. <br>
 *  It returns the oscillation probabilities for Vacuum, Constant matter density, User's Profile density, Earth matter density.
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #include <OscProbability.h>
    #include <OscProbabilityNUN.h>
    #include <OscComplxMatrixExp.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Oscillation Probability Calculation *//

    //* Vaccum Oscillation Function *//
      double OscProbabilityNUN :: 
      VaccumOscProbNUN(int ij, int jk, int type, double Erg, double length)
    {
      Alignment(ij, jk, dim, "Flavour", "VaccumOscProbNUN");
      double** Prob;
      Prob = OscProbabilityNUN :: VaccumOscProbNUN(type, Erg, length);
      return Prob[ij][jk];
      delete Prob;
    }

      double** OscProbabilityNUN ::
      VaccumOscProbNUN(int type, double Erg, double length)
    {
        if(type != 0 && length >= 0.00)
     {
        double** P;
        P = new double*[dim];
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > N(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > M(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > Hvac(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > Norm(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > VacAmp(dim,dim);
        int Type = int (type/(abs(type)));

        N    = OscProbabilityNUN :: SetNonUnitaryPMNS(Type);
        M    = OscProbabilitySTD :: SetMassHamiltonian(Erg);
        Norm = OscProbabilityNUN :: ZeroLengthNorm(N);

        Hvac   = M*length*km2_eV;
        VacAmp = MatExp(img*Hvac);
        VacAmp = N*VacAmp*N.adjoint();

        for(unsigned int ij = 0; ij < dim; ij++)
      {
          P[ij] = new double[dim];
          for(unsigned int jk = 0; jk < dim; jk++)
        {
          P[ij][jk] = pow(abs(VacAmp(ij,jk)),2)/(abs(Norm(ij,ij))*abs(Norm(jk,jk)));
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
      double OscProbabilityNUN :: 
      ConstMatOscProbNUN(int ij, int jk, int type, double Erg, double dens, double length)
    {
      Alignment(ij, jk, dim, "Flavour", "ConstMatOscProbNUN");
      double** Prob;
      Prob = OscProbabilityNUN :: ConstMatOscProbNUN(type, Erg, dens, length);
      return Prob[ij][jk];
      delete Prob;
    }

      double** OscProbabilityNUN ::
      ConstMatOscProbNUN(int type, double Erg, double dens, double length)
    {
        if(type != 0 && dens >= 0.00 && length >= 0.00)
     {
        double** P;
        P = new double*[dim];
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > N(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > M(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > V(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > MatL(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > Norm(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > MatAmp(dim,dim);
        int Type = int (type/(abs(type)));

        N    = OscProbabilityNUN :: SetNonUnitaryPMNS(Type);
        M    = OscProbabilitySTD :: SetMassHamiltonian(Erg);
        V    = OscProbabilitySTD :: GetInteractions(Type*dens);
        Norm = OscProbabilityNUN :: ZeroLengthNorm(N);

        MatL   = (M + N.adjoint()*V*N)*length*km2_eV;
        MatAmp = MatExp(img*MatL);
        MatAmp = N*MatAmp*N.adjoint();

        for(unsigned int ij = 0; ij < dim; ij++)
      {
          P[ij] = new double[dim];
          for(unsigned int jk = 0; jk < dim; jk++)
        {
          P[ij][jk] = pow(abs(MatAmp(ij,jk)),2)/(abs(Norm(ij,ij))*abs(Norm(jk,jk)));
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
      double OscProbabilityNUN ::
      ProfiledMatOscProbNUN(int ij, int jk, int type, double Erg)
    {
      Alignment(ij, jk, dim, "Flavour", "ProfiledMatOscProbNUN");
      double** Prob;
      Prob = OscProbabilityNUN :: ProfiledMatOscProbNUN(type, Erg);
      return Prob[ij][jk];
      delete Prob;
    }
      double** OscProbabilityNUN ::
      ProfiledMatOscProbNUN(int type, double Erg)
   {
        _SYZ_VECTOR_CHECK_(Profile.size(), "SetGivenDensityProfile(double [][2], size_t )");
        if(type != 0)
     {
        double** P;
        P = new double*[dim];
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > N(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > M(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > V(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > MatL(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > Norm(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > MatAmp(dim,dim);
        int Type = int (type/(abs(type)));

        N    = OscProbabilityNUN :: SetNonUnitaryPMNS(Type);
        M    = OscProbabilitySTD :: SetMassHamiltonian(Erg);
        Norm = OscProbabilityNUN :: ZeroLengthNorm(N);

        MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
        for(int i = 0; i < int(Profile.size()); i++)
     {
        V       = OscProbabilitySTD :: GetInteractions(Type*Profile[i].second);
        MatL    = (M + N.adjoint()*V*N)*(Profile[i].first*km2_eV);
        MatAmp *= MatExp(img*MatL);
     }
        MatAmp  = N*MatAmp*N.adjoint();
        for(unsigned int ij = 0; ij < dim; ij++)
      {
          P[ij] = new double[dim];
          for(unsigned int jk = 0; jk < dim; jk++)
        {
          P[ij][jk] = pow(abs(MatAmp(ij,jk)),2)/(abs(Norm(ij,ij))*abs(Norm(jk,jk)));
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
      double OscProbabilityNUN ::
      EarthMatOscProbNUN(int ij, int jk, int type, double Erg)
    {
      Alignment(ij, jk, dim, "Flavour", "EarthMatOscProbNUN");
      double** Prob;
      Prob = OscProbabilityNUN :: EarthMatOscProbNUN(type, Erg);
      return Prob[ij][jk];
      delete Prob;
    }
      double** OscProbabilityNUN ::
      EarthMatOscProbNUN(int type, double Erg)
   {
        _SYZ_VECTOR_CHECK_(Out_In.size(), "SetGivenDensityProfile(double [][2], size_t )");
        if(type != 0)
     {
        double** P;
        P = new double*[dim];
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > N(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > M(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > V(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > MatL(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > Norm(dim,dim);
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > MatAmp(dim,dim);
        int Type = int (type/(abs(type)));

        N    = OscProbabilityNUN :: SetNonUnitaryPMNS(Type);
        M    = OscProbabilitySTD :: SetMassHamiltonian(Erg);
        Norm = OscProbabilityNUN :: ZeroLengthNorm(N);

        MatAmp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim);
        for(int i = 0; i < int(Out_In.size()); i++)
     {
        V       = OscProbabilitySTD :: GetInteractions(Type*Out_In[i].second);
        MatL    = (M + N.adjoint()*V*N)*(Out_In[i].first*km2_eV);
        MatAmp *= MatExp(img*MatL);
     }
        MatAmp  = N*MatAmp*N.adjoint();
        for(unsigned int ij = 0; ij < dim; ij++)
      {
          P[ij] = new double[dim];
          for(unsigned int jk = 0; jk < dim; jk++)
        {
          P[ij][jk] = pow(abs(MatAmp(ij,jk)),2)/(abs(Norm(ij,ij))*abs(Norm(jk,jk)));
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






