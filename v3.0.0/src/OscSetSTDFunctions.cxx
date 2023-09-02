/**
 *
 *  @file  OscSetSTDFunctions.cxx
 *
 *  @brief OscSetSTDFunctions sets all required standard oscillation functions.
 *
 *  For performing the neutrino oscillation, It assigns the values of mass-squared splittings, Mixing Angles, ... etc.
 *  One needs to appropriately assign the oscillation parameters and controlling values via the predefined Struct.
 *  This Loads, User's Profile Density and PREM profile.
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
    //* Setting Oscillation Parameters *//

      void OscProbabilitySTD ::
      SetDimension(const int n)
    {
        if (n > NoMix)
      { dim   = n;
        NMix  = int(0.5*n*(n-1));
        NMass = n-1;
      }
        else if (n < NoMix)
      { 
        std::cout << "\n\033[0;91m:: Error ::\nUnphysical Choice of Flavour = \033[0m" 
                  << n << "\033[0;91m.\033[0m\n"
                  << std::endl;
        std::cout << "\n\033[0;96m:: Suggestion ::\nChoose No. of Flavours >= \033[0m2\n"
                  << std::endl;
        exit(-1);
      }
        else
      { 
        std::cout << "\n\033[0;91m:: Error ::\nWith \033[0m" 
                  << n << " \033[0;91mFlavour, No Mixing Happens !!!\033[0m\n"
                  << std::endl;
        std::cout << "\n\033[0;96m:: Suggestion ::\nAtleast \"Two\" Flavours Require.\033[0m\n"
                  << std::endl;
        exit(0);
      }
    }

      void OscProbabilitySTD ::
      _No_of_Mix_Angles_Check(size_t nt)
    {
      if (int(NMix) == int(nt)){}
        else
      { 
        std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nDimension Ineqaulity Found !!!\n\033[0m"
                  << std::endl;
        std::cout << "\033[0;96m:: Suggestion ::\nFor No. of Flavour = \033[0m" 
                  << dim  << "\033[0;96m,\nThe required No. of Mixing Angles = \033[0m" 
                  << NMix << "\033[0;96m.\033[0m\n" 
                  << std::endl;
        std::cout << "\033[0;96mPlease do check, the No. of Mixing Angles.\033[0m\n" 
                  << std::endl;
        exit(0);
      }
    }

      void OscProbabilitySTD ::
      _No_of_MsSq_Diffs_Check(size_t nt)
    {
      if (int(NMass) == int(nt)){}
        else
      {
        std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nDimension Ineqaulity Found !!!\n\033[0m"
                  << std::endl;
        std::cout << "\033[0;96m:: Suggestion ::\nFor No. of Flavour = \033[0m"
                  << dim   << "\033[0;96m,\nThe required No. of Mass Bases = \033[0m"
                  << NMass << "\033[0;96m.\033[0m\n"
                  << std::endl;
        std::cout << "\033[0;96mPlease do check, the No. of Mass Sq. Differences.\033[0m\n"
                  << std::endl;
        exit(0);
      }
    }

      void OscProbabilitySTD ::
      SetOscFlavours(SetFlavoursPost Flv[], size_t Ft)
    {
        std::vector < int > _flv_;
        if (Ft <= dim && Ft > 1)
      {
          _flv_.clear();
          FlvPost.clear();
          for (unsigned int ij = 0; ij < Ft; ij ++)
        {
          _Std_Flavour_Type_Check(Flv[ij].Fvs);
          Alignment(Flv[ij].Pos, Flv[ij].Pos, dim, "Flavour", "SetOscFlavours");
          _flv_.push_back(Flv[ij].Pos);
          FlvPost.push_back({Flv[ij].Fvs, Flv[ij].Pos});
        }
          int dflv = std::distance(_flv_.begin(), std::unique(_flv_.begin(), _flv_.begin() + _flv_.size()));
          if (Ft - dflv == 0) { }
          else
        {
          std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nIlligal Flavours Settings !!!\n\033[0m"
                    << std::endl;
          std::cout << "\033[0;96mFlavour indices must be Unique.\033[0m\n" 
                    << std::endl;
          std::cout << "\033[0;96mPlease do check, the \" SetFlavoursPost and SetOscFlavours(#, #) \" Settings.\033[0m\n"
                    << std::endl;
          exit(-1);
        }
      }
        else 
      {
        std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nFlavours to Dimension Ineqaulity Found !!!\n\033[0m"
                  << std::endl;
        std::cout << "\033[0;96mIndices must be greater than \" 1 \" and less than or equal to \" \033[0m" << dim << "\033[0;96m \".\033[0m\n" 
                  << std::endl;
        exit(0);
      }
    }

      void OscProbabilitySTD ::
      SetMixAngles(SetMixAngParams n[], size_t nt)
   {
      _No_of_Mix_Angles_Check(nt);
      MixAngs.clear();
      for(unsigned int ij = 0; ij < NMix; ij ++)
    {
      Alignment(n[ij].P1, n[ij].P2, dim, "Mixing-Rotations", "SetMixAngles");
      MixAngs.push_back({n[ij].Th, n[ij].Ph,
                         n[ij].P1, n[ij].P2});
    }
   }

      void OscProbabilitySTD ::
      SetMassSqDiff(SetMassSqParams MSq[], size_t Mt)
   {
      _No_of_MsSq_Diffs_Check(Mt);
      MSqDiff.clear();
      for(unsigned int ij = 0; ij < NMass; ij ++) { MSqDiff.push_back({MSq[ij].Mxx}); }
   }

      void OscProbabilitySTD ::
      SetMediumPlrFrac(double Ye, double Yp, double Yn)
   {
      if(abs(Yn + Yp) == 1.0 && Yn/Yp >= 0.0 && Yn/Yp <= 1.6)
    {
        if(std::copysign(Yp, Ye) == std::copysign(Yp, Yn) && abs(Ye) >= abs(0.9*Yp) && abs(Ye) <= abs(1.1*Yp))
      {
        _ne_ = Ye;
        _np_ = Yp;
        _nn_ = Yn;
      }
        else
      {
        std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed, Unphysical choice !!!\n\033[0m"
                  << std::endl;
        std::cout << "\033[0;96mPlease, Make Sure \"The choice must be either matter or antimatter, and ionic imbalance <= 10%\".\033[0m\n" 
                  << std::endl;
        exit(-1);
      }
    }
      else
    {
      std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed !!!\n\033[0m"
                << std::endl;
      std::cout << "\033[0;96mPlease, Make Sure \"The Sum of Neutron and Proton Molar Fraction = 1, and Ratio <= 1.6\".\033[0m\n" 
                << std::endl;
      exit(0);
    }
   }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilitySTD ::
      GetInteractions(double rho)
   {
       const char *e = "e";
       const char *m = "mu";
       const char *t = "tau";
       if(FlvPost.size() != 0)
     {
       const std::complex < double > CC = sqrt2GF*rho*_ne_;
       const std::complex < double > NC = GFbysqrt2*rho*(Cv_by_Ca*(_np_ - _ne_) - _nn_);

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

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilitySTD ::
      SetPMNS(int type)
   {
      typedef  std::complex < double > r_;
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > R(dim,dim);
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > U(dim,dim);
      U = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim, dim);
      for(unsigned int ij = 0; ij < NMix; ij ++)
    {
      R = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim, dim);

      R(MixAngs[ij].P1-1, MixAngs[ij].P1-1) = r_(cos(rad*MixAngs[ij].Th), 0.);
      R(MixAngs[ij].P2-1, MixAngs[ij].P2-1) = r_(cos(rad*MixAngs[ij].Th), 0.);
      R(MixAngs[ij].P1-1, MixAngs[ij].P2-1) = r_(sin(rad*MixAngs[ij].Th)*cos(rad*MixAngs[ij].Ph), -sin(rad*MixAngs[ij].Th)*sin(rad*MixAngs[ij].Ph));
      R(MixAngs[ij].P2-1, MixAngs[ij].P1-1) = -r_(sin(rad*MixAngs[ij].Th)*cos(rad*MixAngs[ij].Ph), sin(rad*MixAngs[ij].Th)*sin(rad*MixAngs[ij].Ph));
      U *= R;
    }
      if (type == 1) { return U; }
      else           { return U.conjugate();}
   }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilitySTD ::
      SetMassHamiltonian(double Erg)
   {
      _Non_Zero_Non_Neg_Check(Erg, "Energy");
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > M(dim, dim);
      M = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Zero(int(dim), int(dim));
      for(unsigned int ij = 0; ij < NMass; ij ++) { M(int(ij)+1, int(ij)+1) = MSqDiff[ij].Mxx; }
      return M/(2.*Erg*GeV);
   }

    //* Setting Detector Depth from Earth-Surface *//
      void OscProbabilitySTD ::
      SetDetDepth(double depth)
   {
      _Non_Zero_Non_Neg_Check(depth, "Depth");
       if (PREM.size() != 0 && depth <= PREM[int(PREM.size())-1].first) { Depth = depth; }
       else
     {
       std::cout << "\n\033[0;91m:: Error Warning ::\nDepth can't be larger than Earth Radius !!!\n\033[0m\n"
                 << std::endl;
       std::cout << "\033[0;96mElse, Make Sure of Enabling \" SetEarthDensityProfile() \".\033[0m\n" 
                 << std::endl;
       exit(0);
     }
   }

    //* Setting Earth Atomspheric Height Uniformly *//
      void OscProbabilitySTD ::
      SetEarthAtmos(double AtmHeight)
   {
      _Non_Zero_Non_Neg_Check(AtmHeight, "Atmospheric Height");
      AtmosHeight = AtmHeight;
   }

    //* Setting Baseline Matter-Profile *//
      void OscProbabilitySTD ::
      SetGivenDensityProfile(double profile[][2], size_t layers)
   {
      Profile.clear(); 
      for(unsigned int i = 0; i < layers; i++)
    {
       if(profile[i][0] >= 0.00 && profile[i][1] >= 0.00)
     { Profile.push_back(std::make_pair(profile[i][0],profile[i][1])); }
       else
     {
       std::cout << "\n\033[0;91m:: Error Warning ::\nUnphysical Values are entered in Density-Profile\033[0m" 
                 << std::endl;
       exit(0);
     }
    }
   }

    //* Setting Earth Density profile with Zenith *//
      void OscProbabilitySTD ::
      SetZenith(double CosZ) 
   {
        _SYZ_VECTOR_CHECK_(PREM.size(), "SetEarthDensityProfile( )");
        if(CosZ >= -1.0 && CosZ <= 1.0)
     {
        //* oooo000OO000oooo *//
        double d;  // Perpendicular
        double R;  // EarthRadius - Depth;
        double Hr; // EarthRadius + AtmHeight;

        double fL    = 0.0; // Local (Intialize)
        double mL    = 0.0; // Local (Intialize)
        double dens1 = 0.0; // Local (Intialize)
        double dens2 = 0.0; // Local (Intialize)

        double chord1; // Local
        double chord2; // Local

        std::vector < std::pair < double, double > > In;    // Local
        std::vector < std::pair < double, double > > Out;   // Local
        std::vector < std::pair < double, double > > Temp1; // Local
        std::vector < std::pair < double, double > > Temp2; // Local
        //* oooo000OO000oooo *//

        int prem_size = int(PREM.size());
        EarthRadius = PREM[prem_size-1].first;
        R  = EarthRadius - Depth; 
        Hr = EarthRadius + AtmosHeight;
        d  = abs(R*sqrt(1.0-pow(CosZ,2)));

        In.clear();
        Out.clear();
        Temp1.clear();
        Temp2.clear();
        Out_In.clear();

        Temp1.push_back(std::make_pair(Hr,0.0));     // Mandatory(Atmosphere as a "layer" with "0" density);
        if(R != 0.0){                                // If Depth != EarthRadius;
        if(CosZ >= 0.0){
        for(int i = prem_size -1; i >= 0 ; i--){
          if(PREM[i].first > R){
            Temp1.push_back(std::make_pair(PREM[i].first,PREM[i].second));
           }
         }

          for(int i = prem_size -1; i >= 1 ; i--){
             if(R <= PREM[i].first && R > PREM[i-1].first){
                dens1 = PREM[i].second;
               }
            }

          Temp1.push_back(std::make_pair(R, dens1));
          for(int i = 0; i < int(Temp1.size()) -1; i++){
              if(Temp1[i].first >= d){
                  chord1 = chordL(Temp1[i].first, d); 
                  chord2 = chordL(Temp1[i+1].first, d);
                  fL     = Travers_Dist(chord1,chord2);
              if(fL > 0.0){
              Out_In.push_back(std::make_pair(fL, Temp1[i].second));
            }
         }
      }
   }
    else {
      for(int i = prem_size -1; i >= 0 ; i--){
          Temp1.push_back(std::make_pair(PREM[i].first,PREM[i].second));
        }
         for(int i = prem_size -1; i >= 1 ; i--){
           if(R <= PREM[i].first && R > PREM[i-1].first){
              dens1 = PREM[i].second;
             }
          }

        Temp2.push_back(std::make_pair(R, dens1));
        for(int i = prem_size -1; i >= 0 ; i--){
         if(PREM[i].first < R){
            Temp2.push_back(std::make_pair(PREM[i].first,PREM[i].second));
          }
       }
        for(int i = 0; i < int(Temp2.size()) -1; i++){
              if(Temp2[i].first >= d){
                  chord1 = chordL(Temp2[i].first, d);
                  chord2 = chordL(Temp2[i+1].first, d);
                  fL     = Travers_Dist(chord1,chord2);
              if(fL > 0.0){
              Out.push_back(std::make_pair(fL, Temp2[i].second));
            }
         }
      }
        for(int i = 0; i < int(Temp1.size()) -1; i++){
              if(Temp1[i].first >= d){
                  chord1 = chordL(Temp1[i].first, d); 
                  chord2 = chordL(Temp1[i+1].first, d);
                  fL     = Travers_Dist(chord1,chord2);
              if(fL > 0.0){
              In.push_back(std::make_pair(fL, Temp1[i].second));
            }
         }
      }
        if(In[int(In.size())-1].second == Out[int(Out.size())-1].second){
            mL    = In[int(In.size())-1].first + Out[int(Out.size())-1].first;
            dens2 = In[int(In.size())-1].second;
          }
           for(int j = 0; j < int(In.size())-1; j++){
               Out_In.push_back(std::make_pair(In[j].first, In[j].second));
             }
               Out_In.push_back(std::make_pair(mL,dens2));
               for(int j = int(Out.size())-2; j >= 0; j--){
                   Out_In.push_back(std::make_pair(Out[j].first, Out[j].second));
                }
           }
      }
        else {
          for(int i = prem_size -1; i >= 0 ; i--){
              Temp1.push_back(std::make_pair(PREM[i].first,PREM[i].second));
            }
              for(int j = 0; j < int(Temp1.size())-1; j++){
              if(abs(Temp1[j].first - Temp1[j+1].first) > 0.0){
                  Out_In.push_back(std::make_pair(abs(Temp1[j].first - Temp1[j+1].first), Temp1[j].second));
          }
       }
    }
      }
        else
      {
        std::cout << "\n\033[0;91m:: Error Warning ::\nZenith Must Be [-1., 1]\033[0m\n"
                  << std::endl;
        exit(0);
      }
    }

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//









