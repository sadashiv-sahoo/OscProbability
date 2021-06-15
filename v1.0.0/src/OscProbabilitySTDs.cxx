
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
  #include <OscComplxMatrixExp.h>

  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
  //* Oscillation Probability Calculation *//
 
  //* Vaccum Oscillation Function *//
      double NeutrinoOscProbSTD ::
                 VaccumOscProb(int ij, int jk, int type,
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
      double NeutrinoOscProbSTD :: 
                     ConstMatOscProb(int ij, int jk, int type,
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
      double NeutrinoOscProbSTD :: 
                     ProfiledMatOscProb(int ij, int jk, int type, double Erg)
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
      double NeutrinoOscProbSTD :: 
                     EarthMatOscProb(int ij, int jk, int type, double Erg)
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






