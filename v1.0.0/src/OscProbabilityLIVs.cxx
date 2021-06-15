
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
  #include <OscComplxMatrixExp.h>

  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
  //* Oscillation Probability Calculation *//
 
  //* Vaccum Oscillation Function *//
      double NeutrinoOscProbSME ::
                     VaccumOscProbGILIV(int ij, int jk, int type,
                                                double Erg, double length)
   {
        if(type != 0 && length >= 0.00)
     {
        Matrix<complex<double>,dim,dim> U;
        Matrix<complex<double>,dim,dim> M;
        Matrix<complex<double>,dim,dim> A;
        Matrix<complex<double>,dim,dim> C;
        Matrix<complex<double>,dim,dim> HF;
        Matrix<complex<double>,dim,dim> Hvac;
        Matrix<complex<double>,dim,dim> VacAmp;
        int Type = int (type/(abs(type)));

        U    = NeutrinoOscProbSTD :: SetPMNS();
        A    = NeutrinoOscProbSME :: SetGILIVMatrixA(Type);
        C    = NeutrinoOscProbSME :: SetGILIVMatrixC(Type);
        M    = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);

        HF   = U*M*(U.inverse());
        HF  += (Type)*A;
        HF  -= (Erg*1E9)*C;

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
      double NeutrinoOscProbSME :: 
                     ConstMatOscProbGILIV(int ij, int jk, int type,
                                                 double Erg, double dens, double length)
   {
        if(type != 0 && dens >= 0.00 && length >= 0.00)
     {
        Matrix<complex<double>,dim,dim> U;
        Matrix<complex<double>,dim,dim> M;
        Matrix<complex<double>,dim,dim> A;
        Matrix<complex<double>,dim,dim> C;
        Matrix<complex<double>,dim,dim> HF;
        Matrix<complex<double>,dim,dim> HMat;
        Matrix<complex<double>,dim,dim> MatAmp;
        int Type = int (type/(abs(type)));

        U    = NeutrinoOscProbSTD :: SetPMNS();
        A    = NeutrinoOscProbSME :: SetGILIVMatrixA(Type);
        C    = NeutrinoOscProbSME :: SetGILIVMatrixC(Type);
        M    = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);

        HF   = U*M*(U.inverse());
        HF  += (Type)*A;
        HF  -= (Erg*1E9)*C;

        Matrix<complex<double>,dim,dim> Ne;
        typedef complex <double> e_;
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
      double NeutrinoOscProbSME :: 
                     ProfiledMatOscProbGILIV(int ij, int jk, int type, double Erg)
   {
        if(type != 0)
     {
        Matrix<complex<double>,dim,dim> U;
        Matrix<complex<double>,dim,dim> M;
        Matrix<complex<double>,dim,dim> A;
        Matrix<complex<double>,dim,dim> C;
        Matrix<complex<double>,dim,dim> HF;
        Matrix<complex<double>,dim,dim> HmatL;
        Matrix<complex<double>,dim,dim> MatAmp;
        Matrix<complex<double>,dim,dim> Ne;
        typedef complex <double> e_;
        int Type = int (type/(abs(type))); 

        U    = NeutrinoOscProbSTD :: SetPMNS();
        A    = NeutrinoOscProbSME :: SetGILIVMatrixA(Type);
        C    = NeutrinoOscProbSME :: SetGILIVMatrixC(Type);
        M    = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);

        HF   = U*M*(U.inverse());
        HF  += (Type)*A;
        HF  -= (Erg*1E9)*C;

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
      double NeutrinoOscProbSME :: 
                     EarthMatOscProbGILIV(int ij, int jk, int type, double Erg)
   {
        if(type != 0)
     {
        Matrix<complex<double>,dim,dim> U;
        Matrix<complex<double>,dim,dim> M;
        Matrix<complex<double>,dim,dim> A;
        Matrix<complex<double>,dim,dim> C;
        Matrix<complex<double>,dim,dim> HF;
        Matrix<complex<double>,dim,dim> HmatL;
        Matrix<complex<double>,dim,dim> MatAmp;
        Matrix<complex<double>,dim,dim> Ne;
        typedef complex <double> e_;
        int Type = int (type/(abs(type))); 

        U    = NeutrinoOscProbSTD :: SetPMNS();
        A    = NeutrinoOscProbSME :: SetGILIVMatrixA(Type);
        C    = NeutrinoOscProbSME :: SetGILIVMatrixC(Type);
        M    = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);

        HF   = U*M*(U.inverse());
        HF  += (Type)*A;
        HF  -= (Erg*1E9)*C;

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

  //* GVLIV Oscillation Probability *//
      double NeutrinoOscProbSME :: 
                     VaccumOscProbGVLIV(int ij, int jk, double Erg, double length)
   {
         if(ij != jk && abs(ij) > 0 && abs(ij) < 4 && abs(jk) > 0 && abs(jk) < 4)
      {
            if(Erg > 0.0 && length >= 0.0)
         {
              complex <double> PhaseEvolScale;
              Matrix<complex<double>,dim,dim> GV_H_;         // GVLIV Vacuum Proper weighted Hamiltonian;
              Matrix<complex<double>,dim,dim> GVAmp;         // GVLIV Vacuum Oscillation Amplitude;
              Matrix<complex<double>,dim,dim> GVHamiltonian; // GVLIV Vacuum Partial-Hamiltonian; 

              GVHamiltonian = NeutrinoOscProbSME :: GetGVLIVHamiltonian(ij, jk, Erg);
              if(ij > 0 && jk < 0){GV_H_ = (-img*sqrt(2.0))*(GVHamiltonian);}
              else  {GV_H_ = (img*sqrt(2.0))*(GVHamiltonian.conjugate());}             

              PhaseEvolScale = -img*length*km2_eV;
              GVAmp = MatExp(PhaseEvolScale*GV_H_);
                 
              return  pow(abs(GVAmp(abs(ij)-1,abs(jk)-1)),2);                 
         } 
              else 
           {  cout << "\033[0;91mUnphysical Values Entered; Energy & Length must > 0 [Units]\033[0m" 
                   << endl;
                  exit(-1);
           }
      } 
          else 
        { cout << "\033[0;91mOnly Inter-Phase-Transitions are allowed, not Intra-Phase-Transitions...([1, 2, 3] <=> [-1,-2,-3])\033[0m" 
               << endl;
               exit(0);
        }
   }


  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//





