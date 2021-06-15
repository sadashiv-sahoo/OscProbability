
  //////////////////////////////////////////////////////////////////
  //   .........oooOOO 18th-Aug-2020 Last Update OOOooo........   //
  //                        SADASHIV SAHOO                        //
  //               India-based Neutrino Observatory               //
  //        Homi Bhabha National Institute, Mumbai, INDIA         // 
  //          @email: sadashiv.sahoo@tifr.res.in                  //
  //                : sadashiv.sahoo@iopb.res.in                  //
  //////////////////////////////////////////////////////////////////

  #include <OscProbability.h>
  #include <NeutrinoOscProbSTD.h>

  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 
  //* Setting Oscillation Parameters *//
   
  //* Setting MassSqrDiff *//
      void NeutrinoOscProbSTD :: 
                   SetMassSqDiff(double Sol, double Atm)
   {
           m21 = Sol;
           m31 = Atm;
   }

  //* Setting Mixing Angles *//
      void NeutrinoOscProbSTD :: 
                   SetMixAngles(double Th12, double Th23, double Th13, double CP)
   {
           theta12  = Th12*rad;
           theta23  = Th23*rad;
           theta13  = Th13*rad;
           Delta_cp = CP*rad;

           complex <double> deltacp(0., Delta_cp);
           cpterm = exp(deltacp); 
   }

  //* PMNS Matrix Construction *//
      Matrix<complex<double>,dim,dim> NeutrinoOscProbSTD :: SetPMNS()
   {
           Matrix<complex<double>,dim,dim> U12;
           typedef complex <double> u12;
           U12 << u12(cos(theta12),0.),  u12(sin(theta12),0.), u12(0.,0.),
                  u12(-sin(theta12),0.), u12(cos(theta12),0.), u12(0.,0.),
                  u12(0.,0.),            u12(0.,0.),           u12(1.,0.);

           Matrix<complex<double>,dim,dim> U23;
           typedef complex <double> u23;
           U23 << u23(1.,0.), u23(0.,0.),            u23(0.,0.),
                  u23(0.,0.), u23(cos(theta23),0.),  u23(sin(theta23),0.),
                  u23(0.,0.), u23(-sin(theta23),0.), u23(cos(theta23),0.);

           Matrix<complex<double>,dim,dim> U13;
           typedef complex <double> u13;
           U13 << u13(cos(theta13),0.),             u13(0.,0.),   u13(sin(theta13)*cpterm.real(),
                                                                    sin(theta13)*-cpterm.imag()),
                  u13(0.,0.),                       u13(1.,0.),                       u13(0.,0.),
                  u13(-sin(theta13)*cpterm.real(),
                     -sin(theta13)*cpterm.imag()),  u13(0.,0.),             u13(cos(theta13),0.);

           Matrix<complex<double>,dim,dim> P;           // Currently, Dirac Phase is considered
           typedef complex <double> p;
           P  <<  p(1.,0.), p(0.,0.), p(0.,0.),
                  p(0.,0.), p(1.,0.), p(0.,0.),
                  p(0.,0.), p(0.,0.), p(1.,0.); 

           Matrix<complex<double>,dim,dim> U(dim, dim); // PMNS Matrix;
           U = U23*U13*U12*P;

           return U;
   }
  
  //* Hamiltonian at Mass Basis *//
      Matrix<complex<double>,dim,dim> NeutrinoOscProbSTD :: SetMassHamiltonian(double Erg)
   {
        if(Erg > 0.0)
     {
        Matrix<complex<double>,dim,dim> M;
        typedef complex <double> m;
        double Energy = abs(Erg)*1E9;

        M  <<  m(0.,0.),  m(0.,0.),                           m(0.,0.),
               m(0.,0.),  m(m21/(2.*Energy),0.),              m(0.,0.),
               m(0.,0.),  m(0.,0.),              m(m31/(2.*Energy),0.);

        return M;
     }
        else 
      { cout << "\033[0;91mError :: Unphysical Value, Energy shoud have +ve\033[0m" 
             << endl;
             exit(0);
      }
   }

  //* Setting Baseline Matter-Profile *//
      void NeutrinoOscProbSTD :: 
                   SetGivenDensityProfile(double profile[][2], size_t layers)
       {
             Profile.clear(); 
             for(unsigned int i = 0; i < layers; i++)
           {
               if(profile[i][0] >= 0.00 && profile[i][1] >= 0.00)
             {
                Profile.push_back(make_pair(profile[i][0],profile[i][1]));
             }
                else
              {
                cout << "\033[0;91mError :: Unphysical Values are entered in Density-Profile\033[0m" 
                     << endl;
                     exit(0);
              }
           } 
       }

  //* Setting Detector Depth from Earth-Surface *//
      void NeutrinoOscProbSTD :: 
                   SetDetDepth(double depth)
   {
        if(depth >= 0.0)
     {
           if(abs(depth) <= PREM[int(PREM.size())-1].first)
        {
            Depth = abs(depth);
        }   
            else 
          { cout << "\033[0;91mError :: Depth can't be larger than Earth Radius\033[0m" 
                 << endl;
                 exit(-1);
          } 
     }   
          else 
        { cout << "\033[0;91mError :: Unphysical Value, Depth should have +ve\033[0m" 
                   << endl;
                   exit(0);
        }
   }

  //* Setting Earth Atomspheric Height Uniformly *//
      void NeutrinoOscProbSTD :: 
                   SetEarthAtmos(double AtmHeight)
   {
        if(AtmHeight >= 0.0)
     {
        AtmosHeight = abs(AtmHeight);
     }
        else 
      { cout << "\033[0;91mError :: Unphysical Value, AtmosHeight should have +ve\033[0m" 
             << endl;
             exit(0);
      }
   }

  //* Propagation Length in Earth Profile *//
      double chordL(double R, double d)
   {
      if(R >= d) { return abs(sqrt(pow(R,2)-pow(d,2))); }
      else       { return 0; }
   }

  //* Concentric Propagations Choice *//
      double Travers_Dist(double l1, double l2) { return abs(l1 -l2); }

  //* Setting Earth Density profile with Zenith *//
      void NeutrinoOscProbSTD :: SetZenith(double CosZ) 
   {
        if(CosZ >= -1. && CosZ <= 1.)
     {
        double zenith = CosZ;

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

        vector < pair < double, double > > In;    // Local
        vector < pair < double, double > > Out;   // Local
        vector < pair < double, double > > Temp1; // Local
        vector < pair < double, double > > Temp2; // Local
        //* oooo000OO000oooo *//

        int prem_size = int(PREM.size());
        EarthRadius = PREM[prem_size-1].first;
        R  = EarthRadius - Depth; 
        Hr = EarthRadius + AtmosHeight;
        d  = abs(R*sqrt(1.0-pow(zenith,2)));

        In.clear();
        Out.clear();
        Temp1.clear();
        Temp2.clear();
        Out_In.clear();

        Temp1.push_back(make_pair(Hr,0.0));          // Mandatory(Atmosphere as a "layer" with "0" density);
        if(R != 0.0){                                // If Depth != EarthRadius;
        if(zenith >= 0.0){
        for(int i = prem_size -1; i >= 0 ; i--){
          if(PREM[i].first > R){
            Temp1.push_back(make_pair(PREM[i].first,PREM[i].second));
           }
         }

          for(int i = prem_size -1; i >= 1 ; i--){
             if(R <= PREM[i].first && R > PREM[i-1].first){
                dens1 = PREM[i].second;
               }            
            }

          Temp1.push_back(make_pair(R, dens1));
          for(int i = 0; i < int(Temp1.size()) -1; i++){
              if(Temp1[i].first >= d){
                  chord1 = chordL(Temp1[i].first, d); 
                  chord2 = chordL(Temp1[i+1].first, d);
                  fL     = Travers_Dist(chord1,chord2);
              if(fL > 0.0){
              Out_In.push_back(make_pair(fL, Temp1[i].second));
            }
         }
      }
   }
    else {
      for(int i = prem_size -1; i >= 0 ; i--){     
          Temp1.push_back(make_pair(PREM[i].first,PREM[i].second));
        }
         for(int i = prem_size -1; i >= 1 ; i--){
           if(R <= PREM[i].first && R > PREM[i-1].first){
              dens1 = PREM[i].second;
             }            
          }

        Temp2.push_back(make_pair(R, dens1));     
        for(int i = prem_size -1; i >= 0 ; i--){
         if(PREM[i].first < R){
            Temp2.push_back(make_pair(PREM[i].first,PREM[i].second));
          }
       }
        for(int i = 0; i < int(Temp2.size()) -1; i++){
              if(Temp2[i].first >= d){
                  chord1 = chordL(Temp2[i].first, d); 
                  chord2 = chordL(Temp2[i+1].first, d);
                  fL     = Travers_Dist(chord1,chord2);
              if(fL > 0.0){
              Out.push_back(make_pair(fL, Temp2[i].second));
            }
         }
      }
        for(int i = 0; i < int(Temp1.size()) -1; i++){
              if(Temp1[i].first >= d){
                  chord1 = chordL(Temp1[i].first, d); 
                  chord2 = chordL(Temp1[i+1].first, d);
                  fL     = Travers_Dist(chord1,chord2);
              if(fL > 0.0){
              In.push_back(make_pair(fL, Temp1[i].second));
            }
         }
      }
        if(In[int(In.size())-1].second == Out[int(Out.size())-1].second){
            mL    = In[int(In.size())-1].first + Out[int(Out.size())-1].first;
            dens2 = In[int(In.size())-1].second;
          }
           for(int j = 0; j < int(In.size())-1; j++){
               Out_In.push_back(make_pair(In[j].first, In[j].second));
             }
               Out_In.push_back(make_pair(mL,dens2));
               for(int j = int(Out.size())-2; j >= 0; j--){
                   Out_In.push_back(make_pair(Out[j].first, Out[j].second));
                }
           }
      }
        else {
          for(int i = prem_size -1; i >= 0 ; i--){
              Temp1.push_back(make_pair(PREM[i].first,PREM[i].second));
            }
              for(int j = 0; j < int(Temp1.size())-1; j++){
              if(abs(Temp1[j].first - Temp1[j+1].first) > 0.0){
                  Out_In.push_back(make_pair(abs(Temp1[j].first - Temp1[j+1].first), Temp1[j].second));
          }
       }
    }
//      for(int i = 0; i < int(Out_In.size()); i++){
//         cout << i+1 << "\t" << Out_In[i].first << "\t" << Out_In[i].second << endl;
//      }
          } else { cout << "\n\033[0;91mError :: Zenith Must Be [-1., 1]\033[0m\n" << endl;
                   exit(0);
       }

    }

  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//  





