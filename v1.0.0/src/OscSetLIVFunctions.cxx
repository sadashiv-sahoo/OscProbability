
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
  //* Setting SME-GILIV Parameters *//
  
  //* SME-GILIV Controls over Parameters *//
      void NeutrinoOscProbSME :: SetGILIVControls(int _mu_, int _nu_)
   {
        if(_mu_ >= 0 && _mu_ < 4 && _nu_ >= 0 && _nu_ < 4)
     {
          int i = 1; // (initiation);
          int j = 1; // (initiation);

          if(_nu_ == 0){j = 1;} else {j = -1;}
          if(_mu_ == 0){i = 1; aL = 1.0;} else {i = -1; aL = -1.0;}
          if(_mu_ == 0 && _nu_ == 0){cL = (4.0/3.0);} else {cL = 1.0;}
          cL *= i*j; // cout << "aL = " << aL << "\t" << "cL = " << cL << endl;
     }   
         else 
       { cout << "\n\033[0;91mError :: Four-Momenta Index (0,1,2,3) \033[0m\n" 
              << endl;
              exit(0); 
       }
   }

  //* GILIV ODD Parameters *//
      void NeutrinoOscProbSME :: 
                   SetGILIVaX( complex <double> _ee, complex <double> _mumu,
                                complex <double> _tautau, complex <double> _emu, 
                                   complex <double> _etau, complex <double> _mutau)
   {
          aX_ee         = _ee*1E9;      // in GeV;
          aX_mumu       = _mumu*1E9;    // in GeV;
          aX_tautau     = _tautau*1E9;  // in GeV;
          aX_emu        = _emu*1E9;     // in GeV;
          aX_etau       = _etau*1E9;    // in GeV;
          aX_mutau      = _mutau*1E9;   // in GeV;
   }

  //* GILIV EVEN Parameters *//
      void NeutrinoOscProbSME :: 
                   SetGILIVcXX( complex <double> _ee_, complex <double> _mumu_,
                                 complex <double> _tautau_, complex <double> _emu_, 
                                    complex <double> _etau_, complex <double> _mutau_)
   {
          cXX_ee        = _ee_;         // Without Unit;
          cXX_mumu      = _mumu_;       // Without Unit;
          cXX_tautau    = _tautau_;     // Without Unit;
          cXX_emu       = _emu_;        // Without Unit;
          cXX_etau      = _etau_;       // Without Unit;
          cXX_mutau     = _mutau_;      // Without Unit;
   }

  //* GILIV ODD Matrix *//
      Matrix<complex<double>,dim,dim> NeutrinoOscProbSME :: SetGILIVMatrixA()
   {
          Matrix<complex<double>,dim,dim> LIVMatrixA;
          LIVMatrixA << aX_ee,         aX_emu,         aX_etau,
                        conj(aX_emu),  aX_mumu,        aX_mutau,
                        conj(aX_etau), conj(aX_mutau), aX_tautau;
          LIVMatrixA *= aL;
          return LIVMatrixA;
   }

  //* GILIV EVEN Matrix *//
      Matrix<complex<double>,dim,dim> NeutrinoOscProbSME :: SetGILIVMatrixC()
   {
          Matrix<complex<double>,dim,dim> LIVMatrixC;
          LIVMatrixC << cXX_ee,         cXX_emu,         cXX_etau,
                        conj(cXX_emu),  cXX_mumu,        cXX_mutau,
                        conj(cXX_etau), conj(cXX_mutau), cXX_tautau;
          LIVMatrixC *= cL;
          return LIVMatrixC;
   }

  //* GILIV ODD Matrix* *//
      Matrix<complex<double>,dim,dim> NeutrinoOscProbSME :: SetGILIVMatrixA(int type)
   {
          Matrix<complex<double>,dim,dim> LIVMatrixAType;
          if(type == 1) { LIVMatrixAType = NeutrinoOscProbSME :: SetGILIVMatrixA(); } 
          else { LIVMatrixAType = NeutrinoOscProbSME :: SetGILIVMatrixA().conjugate(); }
          return LIVMatrixAType;
   }

  //* GILIV EVEN Matrix* *//
      Matrix<complex<double>,dim,dim> NeutrinoOscProbSME :: SetGILIVMatrixC(int type)
   {
          Matrix<complex<double>,dim,dim> LIVMatrixCType;
          if(type == 1){ LIVMatrixCType = NeutrinoOscProbSME :: SetGILIVMatrixC(); } 
          else { LIVMatrixCType = NeutrinoOscProbSME :: SetGILIVMatrixC().conjugate(); }
          return LIVMatrixCType;
   }


  //* Setting SME-GVLIV Parameters *//
  //* SME-GVLIV Controls over Parameters *//
      void NeutrinoOscProbSME :: SetGVLIVControls(int _mu_, int _sigma_)
   {
        if(_mu_ >= 0 && _mu_ < 4 && _sigma_ >= 0 && _sigma_ < 4)
     {
         int i = 1; // (initiation);
         int j = 1; // (initiation);

         if(_sigma_ == 0){j = 1;} else {j = -1;}
         if(_mu_ == 0){i = 1; _h_ = 1.0;} else {i = -1; _h_ = -1.0;}
         _g_  = i*j;
     }  
         else 
      {  cout << "\n\033[0;91mError :: Four-Momenta Index (0,1,2,3) \033[0m\n" 
              << endl;
              exit(0); 
      }

   }

  //* GVLIV Vacuum Direction Polarization *//
  //* WithOut Spatial Polarization, Gauge-Violation Mixing is Highly Suprressed at Low Energy !!! *//
      void NeutrinoOscProbSME :: SetGVLIVPolzation(complex <double> eps) 
   {                                                                     
        if(abs(eps) != 0.0){ epsilon = eps/(abs(eps)); }
        else { cout << "\033[0;91mVacuum Direction Polarization is necessary for Gauge Violating Oscillations\033[0m" 
                    << endl;
                   exit(-1);
             }
  }

  //* GVLIV ODD Parameters *//
      void NeutrinoOscProbSME :: 
                   SetGVLIVgXXX(complex <double> _eebar, complex <double> _mumubar,
                                   complex <double> _tautaubar, complex <double> _emubar, 
                                      complex <double> _etaubar, complex <double> _mutaubar)
   {
        gXXX_eebar     = _eebar*1.0E9;
        gXXX_mumubar   = _mumubar*1.0E9;         
        gXXX_tautaubar = _tautaubar*1.0E9;
        gXXX_emubar    = _emubar*1.0E9;
        gXXX_etaubar   = _etaubar*1.0E9;
        gXXX_mutaubar  = _mutaubar*1.0E9;
   }

  //* GVLIV EVEN Parameters *//
      void NeutrinoOscProbSME ::
                   SetGVLIVhXX(complex <double> _eebar_, complex <double> _mumubar_,
                                  complex <double> _tautaubar_, complex <double> _emubar_, 
                                     complex <double> _etaubar_, complex <double> _mutaubar_)
   {
        hXX_eebar     = _eebar_;
        hXX_mumubar   = _mumubar_;         
        hXX_tautaubar = _tautaubar_;
        hXX_emubar    = _emubar_;
        hXX_etaubar   = _etaubar_;
        hXX_mutaubar  = _mutaubar_;
   }

  //* GVLIV ODD Matrix *//
      Matrix<complex<double>,dim,dim>  NeutrinoOscProbSME :: SetGVLIVMatrixG()
   {
        Matrix<complex<double>,dim,dim> LIVMatrixG;
LIVMatrixG << gXXX_eebar,         gXXX_emubar,         gXXX_etaubar,
                      conj(gXXX_emubar),  gXXX_mumubar,        gXXX_mutaubar,
                      conj(gXXX_etaubar), conj(gXXX_mutaubar), gXXX_tautaubar;
        LIVMatrixG *= _g_;
        return LIVMatrixG;
   }    

  //* GVLIV EVEN Matrix *//
      Matrix<complex<double>,dim,dim> NeutrinoOscProbSME :: SetGVLIVMatrixH()
   {
        Matrix<complex<double>,dim,dim> LIVMatrixH;
        LIVMatrixH << hXX_eebar,         hXX_emubar,         hXX_etaubar,
                      conj(hXX_emubar),  hXX_mumubar,        hXX_mutaubar,
                      conj(hXX_etaubar), conj(hXX_mutaubar), hXX_tautaubar;
        LIVMatrixH *= _h_;
        return LIVMatrixH;
   }

  //* GVLIV Charge-Conguation Matrix *//    
      Matrix<complex<double>,dim,dim> NeutrinoOscProbSME :: SetGVLIVMatrixC()
   {
        Matrix<complex<double>,dim,dim> I;
        I = Matrix<complex<double>,dim,dim> :: Identity(); // Till now It acts as Idenitity Matrix;
        return I;                                          // One may modify for further Research Advances;
   }

  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//






 
