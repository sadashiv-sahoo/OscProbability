
  //////////////////////////////////////////////////////////////////
  //   .........oooOOO 18th-Aug-2020 Last Update OOOooo........   //
  //                        SADASHIV SAHOO                        //
  //               India-based Neutrino Observatory               //
  //        Homi Bhabha National Institute, Mumbai, INDIA         // 
  //          @email: sadashiv.sahoo@tifr.res.in                  //
  //                : sadashiv.sahoo@iopb.res.in                  //
  //////////////////////////////////////////////////////////////////

  #include <OscProbability.h>
  #include <NeutrinoOscProbSME.h>

  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 
  //* Currently this decay phenomena is Not-Supported by Author *//

  //* Setting Non-Standard Neutrino Decay Parameters *//
  //* NS-NuDecay Controls Parameters *//
      void NeutrinoOscProbDEC ::
                   SetNuDecay( complex <double> a_11, complex <double> a_12, complex <double> a_13, 
                                  complex <double> a_21, complex <double> a_22, complex <double> a_23,
                                    complex <double> a_31, complex <double> a_32, complex <double> a_33 )
   {
          _a_11_        = a_11;  // _a_11_ -Mass Basis;
          _a_12_        = a_12;  // _a_12_ -Mass Basis;
          _a_13_        = a_13;  // _a_13_ -Mass Basis;
          _a_21_        = a_21;  // _a_21_ -Mass Basis;
          _a_22_        = a_22;  // _a_22_ -Mass Basis;
          _a_23_        = a_23;  // _a_23_ -Mass Basis;
          _a_31_        = a_31;  // _a_31_ -Mass Basis;
          _a_32_        = a_32;  // _a_32_ -Mass Basis;
          _a_33_        = a_33;  // _a_33_ -Mass Basis;
   }

  //* NS NuDecay Matrix *//
      Matrix<complex<double>,dim,dim> NeutrinoOscProbDEC :: SetNuDecayMatrix(double Erg)
   {
       if(Erg > 0.0)
     {
         Matrix<complex<double>,dim,dim> NuDecay;
         NuDecay << _a_11_, _a_12_, _a_13_,
                    _a_21_, _a_22_, _a_23_,
                    _a_31_, _a_32_, _a_33_;
         NuDecay *= (1.0/(2.0*(abs(Erg))*1E9));
         return NuDecay;
     }
         else 
       { cout << "\033[0;91mError :: Unphysical Value, Energy shoud have +ve\033[0m" 
              << endl;
              exit(0);
       }
   }

  //* Vaccum Oscillation Function *//
      double NeutrinoOscProbDEC ::
                     VaccumOscProbNuDecay(int ij, int jk, int type,
                                                  double Erg, double length)
   {
        if(type != 0 && length >= 0.00)
     {
        Matrix<complex<double>,dim,dim> U;
        Matrix<complex<double>,dim,dim> M;
        Matrix<complex<double>,dim,dim> HF;
        Matrix<complex<double>,dim,dim> Hvac;
        Matrix<complex<double>,dim,dim> VacAmp;
        int Type = int (type/(abs(type)));

        U  = NeutrinoOscProbSTD :: SetPMNS();
        M  = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);
        M += NeutrinoOscProbDEC :: SetNuDecayMatrix(Erg);
        HF = (U*M*U.inverse());

        Hvac   = (-Type)*HF*length*km2_eV;
        VacAmp = MatExp(-img*Hvac);

        return  pow(abs(VacAmp(ij,jk)),2);
     }
        else 
      { cout << "\n\033[0;91mError :: Unphysical Value(s) (Negative or Zero) is(are) entered !!!\033[0m\n" 
             << endl;
             exit(0);
      }

   }

  //* Const Matter Density Oscillation Function *//
      double NeutrinoOscProbDEC :: 
                     ConstMatOscProbNuDecay(int ij, int jk, int type,
                                                 double Erg, double dens, double length)
   {
        if(type != 0 && dens >= 0.00 && length >= 0.00)
     {
        Matrix<complex<double>,dim,dim> U;
        Matrix<complex<double>,dim,dim> M;
        Matrix<complex<double>,dim,dim> HF;
        Matrix<complex<double>,dim,dim> HMat;
        Matrix<complex<double>,dim,dim> MatAmp;
        Matrix<complex<double>,dim,dim> Ne;
        typedef complex <double> e_;
        int Type = int (type/(abs(type)));

        U  = NeutrinoOscProbSTD :: SetPMNS();
        M  = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);
        M += NeutrinoOscProbDEC :: SetNuDecayMatrix(Erg);
        HF = (U*M*U.inverse());

        Ne << e_(Type*sqrt2GF*dens*eS,0.), e_(0.,0.), e_(0.,0.),
                                e_(0.,0.), e_(0.,0.), e_(0.,0.),
                                e_(0.,0.), e_(0.,0.), e_(0.,0.);
        HMat   = (-Type)*(HF+Ne)*(length*km2_eV);
        MatAmp = MatExp(-img*HMat);

        return  pow(abs(MatAmp(ij,jk)),2);
     }
        else
      { cout << "\n\033[0;91mError :: Unphysical Value(s) (Negative or Zero) is(are) entered !!!\033[0m\n" 
             << endl;
             exit(0);
      }
   }

  //* Profiled Matter Density Oscillation Function *//
      double NeutrinoOscProbDEC :: 
                     ProfiledMatOscProbNuDecay(int ij, int jk, int type, double Erg)
   {
        if(type != 0)
     {
        Matrix<complex<double>,dim,dim> U;
        Matrix<complex<double>,dim,dim> M;
        Matrix<complex<double>,dim,dim> HF;
        Matrix<complex<double>,dim,dim> HmatL;
        Matrix<complex<double>,dim,dim> MatAmp;
        Matrix<complex<double>,dim,dim> Ne;
        typedef complex <double> e_;
        int Type = int (type/(abs(type)));

        U        = NeutrinoOscProbSTD :: SetPMNS();
        M        = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);
        M       += NeutrinoOscProbDEC :: SetNuDecayMatrix(Erg);
        HF       = (U*M*U.inverse());

          MatAmp   = Matrix<complex<double>,dim,dim> :: Identity();
          for(int i = 0; i < int(Profile.size()); i++)
       {
          Ne << e_(Type*sqrt2GF*(Profile[i].second)*eS,0.), e_(0.,0.), e_(0.,0.),  
                                                 e_(0.,0.), e_(0.,0.), e_(0.,0.),
                                                 e_(0.,0.), e_(0.,0.), e_(0.,0.);
          HmatL   = (-Type)*(HF+Ne)*((Profile[i].first)*km2_eV);
          MatAmp *= MatExp(-img*HmatL); 
       }
          return pow(abs(MatAmp(ij,jk)),2);
     }
          else
       {  cout << "\n\033[0;91mError :: Type should be integer(-1 or +1) and not equal to 0\033[0m\n"
               << endl;
               exit(0);
       }
   }

  //* Earth Matter Density Profile Oscillation Function *//
      double NeutrinoOscProbDEC :: 
                     EarthMatOscProbNuDecay(int ij, int jk, int type, double Erg)
   {
        if(type != 0)
     {
        Matrix<complex<double>,dim,dim> U;
        Matrix<complex<double>,dim,dim> M;
        Matrix<complex<double>,dim,dim> HF;
        Matrix<complex<double>,dim,dim> HmatL;
        Matrix<complex<double>,dim,dim> MatAmp;
        Matrix<complex<double>,dim,dim> Ne;
        typedef complex <double> e_;
        int Type = int (type/(abs(type)));

        U        = NeutrinoOscProbSTD :: SetPMNS();
        M        = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);
        M       += NeutrinoOscProbDEC :: SetNuDecayMatrix(Erg);
        HF       = (U*M*U.inverse());

          MatAmp   = Matrix<complex<double>,dim,dim> :: Identity();
          for(int i = 0; i < int(Out_In.size()); i++)
       {
          Ne << e_(Type*sqrt2GF*(Out_In[i].second)*eS,0.), e_(0.,0.), e_(0.,0.),
                                                e_(0.,0.), e_(0.,0.), e_(0.,0.),
                                                e_(0.,0.), e_(0.,0.), e_(0.,0.);
          HmatL   = (-Type)*(HF+Ne)*((Out_In[i].first)*km2_eV);
          MatAmp *= MatExp(-img*HmatL);
       }
          return pow(abs(MatAmp(ij,jk)),2);
    }
          else
       {  cout << "\n\033[0;91mError :: Type should be integer(-1 or +1) and not equal to 0\033[0m\n"
               << endl;
               exit(0);
       }
   }


  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//






