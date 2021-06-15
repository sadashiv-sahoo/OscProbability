
  //////////////////////////////////////////////////////////////////
  //   .........oooOOO 18th-Aug-2020 Last Update OOOooo........   //
  //                        SADASHIV SAHOO                        //
  //               India-based Neutrino Observatory               //
  //        Homi Bhabha National Institute, Mumbai, INDIA         // 
  //          @email: sadashiv.sahoo@tifr.res.in                  //
  //                : sadashiv.sahoo@iopb.res.in                  //
  //////////////////////////////////////////////////////////////////

  #include <OscProbability.h>
  #include <NeutrinoOscProbNSI.h>
  #include <OscComplxMatrixExp.h>

  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
  //* Oscillation Probability Calculation *// 

  //* Const Matter Density Oscillation Function *//
      double NeutrinoOscProbNSI :: 
                     ConstMatOscProbNSI(int ij, int jk, int type,
                                                 double Erg, double dens, double length)
   {
        if(type != 0 && dens >= 0.00 && length >= 0.00)
     {
        double _CCe_ = dens*eS;
        double _CCu_ = rS*_CCe_;
        double _CCd_ = _CCu_;

        Matrix<complex<double>,dim,dim> U;
        Matrix<complex<double>,dim,dim> H;
        Matrix<complex<double>,dim,dim> HF;
        Matrix<complex<double>,dim,dim> HMat;
        Matrix<complex<double>,dim,dim> MatAmp;
        Matrix<complex<double>,dim,dim> NSI_E;
        Matrix<complex<double>,dim,dim> NSI_U;
        Matrix<complex<double>,dim,dim> NSI_D;
        Matrix<complex<double>,dim,dim> NSI_TOTAL;
        Matrix<complex<double>,dim,dim> Ne;
        typedef complex <double> e_;
        int Type = int (type/(abs(type)));

        U     = NeutrinoOscProbSTD :: SetPMNS();
        NSI_E = NeutrinoOscProbNSI :: SetNSIMatrixE();
        NSI_U = NeutrinoOscProbNSI :: SetNSIMatrixU();
        NSI_D = NeutrinoOscProbNSI :: SetNSIMatrixD();
        H     = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);
        HF    = U*H*(U.inverse());

        Ne << e_(Type*sqrt2GF*dens*eS,0.), e_(0.,0.), e_(0.,0.),
                                e_(0.,0.), e_(0.,0.), e_(0.,0.),
                                e_(0.,0.), e_(0.,0.), e_(0.,0.);

        NSI_TOTAL  = ((Type*sqrt2GF*_CCe_))*NSI_E;
        NSI_TOTAL += ((Type*sqrt2GF*_CCu_))*NSI_U;
        NSI_TOTAL += ((Type*sqrt2GF*_CCd_))*NSI_D;

        HMat   = (-Type)*(HF+Ne+NSI_TOTAL)*(length*km2_eV);
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
      double NeutrinoOscProbNSI :: 
                     ProfiledMatOscProbNSI(int ij, int jk, int type, double Erg)
   {
        if(type != 0)
     {
        double _CCe_;
        double _CCu_;
        double _CCd_;

        Matrix<complex<double>,dim,dim> U;
        Matrix<complex<double>,dim,dim> H;
        Matrix<complex<double>,dim,dim> HF;
        Matrix<complex<double>,dim,dim> HmatL;
        Matrix<complex<double>,dim,dim> MatAmp;
        Matrix<complex<double>,dim,dim> NSI_E;
        Matrix<complex<double>,dim,dim> NSI_U;
        Matrix<complex<double>,dim,dim> NSI_D;
        Matrix<complex<double>,dim,dim> NSI_TOTAL;
        Matrix<complex<double>,dim,dim> Ne;
        typedef complex <double> e_;
        int Type = int (type/(abs(type)));

        U     = NeutrinoOscProbSTD :: SetPMNS();
        NSI_E = NeutrinoOscProbNSI :: SetNSIMatrixE();
        NSI_U = NeutrinoOscProbNSI :: SetNSIMatrixU();
        NSI_D = NeutrinoOscProbNSI :: SetNSIMatrixD();
        H     = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);
        HF    = U*H*(U.inverse());

          MatAmp   = Matrix<complex<double>,dim,dim> :: Identity();
          for(int i = 0; i < int(Profile.size()); i++)
       {
          Ne << e_(Type*sqrt2GF*(Profile[i].second)*eS,0.), e_(0.,0.), e_(0.,0.),
                                                 e_(0.,0.), e_(0.,0.), e_(0.,0.),
                                                 e_(0.,0.), e_(0.,0.), e_(0.,0.);
          _CCe_ = (Profile[i].second)*eS;
          _CCu_ = rS*_CCe_;
          _CCd_ = _CCu_;

          NSI_TOTAL  = ((Type*sqrt2GF*_CCe_))*NSI_E;
          NSI_TOTAL += ((Type*sqrt2GF*_CCu_))*NSI_U;
          NSI_TOTAL += ((Type*sqrt2GF*_CCd_))*NSI_D;

          HmatL   = (-Type)*(HF+Ne+NSI_TOTAL)*((Profile[i].first)*km2_eV);
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
      double NeutrinoOscProbNSI :: 
                     EarthMatOscProbNSI(int ij, int jk, int type, double Erg)
   {
        if(type != 0)
     {
        double _CCe_;
        double _CCu_;
        double _CCd_;

        Matrix<complex<double>,dim,dim> U;
        Matrix<complex<double>,dim,dim> H;
        Matrix<complex<double>,dim,dim> HF;
        Matrix<complex<double>,dim,dim> HmatL;
        Matrix<complex<double>,dim,dim> MatAmp;
        Matrix<complex<double>,dim,dim> NSI_E;
        Matrix<complex<double>,dim,dim> NSI_U;
        Matrix<complex<double>,dim,dim> NSI_D;
        Matrix<complex<double>,dim,dim> NSI_TOTAL;
        Matrix<complex<double>,dim,dim> Ne;
        typedef complex <double> e_;
        int Type = int (type/(abs(type)));

        U     = NeutrinoOscProbSTD :: SetPMNS();
        NSI_E = NeutrinoOscProbNSI :: SetNSIMatrixE();
        NSI_U = NeutrinoOscProbNSI :: SetNSIMatrixU();
        NSI_D = NeutrinoOscProbNSI :: SetNSIMatrixD();
        H     = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);
        HF    = U*H*(U.inverse());

          MatAmp   = Matrix<complex<double>,dim,dim> :: Identity();
          for(int i = 0; i < int(Out_In.size()); i++)
       {
          Ne << e_(Type*sqrt2GF*(Out_In[i].second)*eS,0.), e_(0.,0.), e_(0.,0.),
                                                e_(0.,0.), e_(0.,0.), e_(0.,0.),
                                                e_(0.,0.), e_(0.,0.), e_(0.,0.);
          _CCe_ = (Out_In[i].second)*eS;
          _CCu_ = rS*_CCe_;
          _CCd_ = _CCu_;

          NSI_TOTAL  = ((Type*sqrt2GF*_CCe_))*NSI_E;
          NSI_TOTAL += ((Type*sqrt2GF*_CCu_))*NSI_U;
          NSI_TOTAL += ((Type*sqrt2GF*_CCd_))*NSI_D;

          HmatL   = (-Type)*(HF+Ne+NSI_TOTAL)*((Out_In[i].first)*km2_eV);
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





