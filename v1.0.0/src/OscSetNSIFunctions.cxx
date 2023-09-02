
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

  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *// 
  //* Setting NSI Parameters *//

  //* NSI Parameters for electron *//
      void NeutrinoOscProbNSI :: 
                   SetNSI_E( complex <double> E_ee_, complex <double> E_mumu_,
                                complex <double> E_tautau_, complex <double> E_emu_, 
                                   complex <double> E_etau_, complex <double> E_mutau_)
   {
            nsiE_ee     = E_ee_;
            nsiE_mumu   = E_mumu_;
            nsiE_tautau = E_tautau_;
            nsiE_emu    = E_emu_;
            nsiE_etau   = E_etau_;
            nsiE_mutau  = E_mutau_;
   }

  //* NSI Parameters for u-quark *//
      void NeutrinoOscProbNSI :: 
                   SetNSI_U( complex <double> U_ee_, complex <double> U_mumu_,
                                complex <double> U_tautau_, complex <double> U_emu_, 
                                   complex <double> U_etau_, complex <double> U_mutau_)
   {
            nsiU_ee     = U_ee_;
            nsiU_mumu   = U_mumu_;
            nsiU_tautau = U_tautau_;
            nsiU_emu    = U_emu_;
            nsiU_etau   = U_etau_;
            nsiU_mutau  = U_mutau_;
   }

  //* NSI Parameters for d-quark *//
      void NeutrinoOscProbNSI :: 
                   SetNSI_D( complex <double> D_ee_, complex <double> D_mumu_,
                                complex <double> D_tautau_, complex <double> D_emu_, 
                                   complex <double> D_etau_, complex <double> D_mutau_)
   {
            nsiD_ee     = D_ee_;
            nsiD_mumu   = D_mumu_;
            nsiD_tautau = D_tautau_;
            nsiD_emu    = D_emu_;
            nsiD_etau   = D_etau_;
            nsiD_mutau  = D_mutau_;
   }

  //* NSI Matrix for electron *//
      Matrix<complex<double>,dim,dim> NeutrinoOscProbNSI :: SetNSIMatrixE()
   {
            Matrix<complex<double>,dim,dim> NSIMatrixE;
            NSIMatrixE << nsiE_ee,         nsiE_emu,         nsiE_etau,
                          conj(nsiE_emu),  nsiE_mumu,        nsiE_mutau,
                          conj(nsiE_etau), conj(nsiE_mutau), nsiE_tautau;
            return NSIMatrixE;
   }

  //* NSI Matrix for u-quark *//
      Matrix<complex<double>,dim,dim> NeutrinoOscProbNSI :: SetNSIMatrixU()
   {
            Matrix<complex<double>,dim,dim> NSIMatrixU;
            NSIMatrixU << nsiU_ee,         nsiU_emu,         nsiU_etau,
                          conj(nsiU_emu),  nsiU_mumu,        nsiU_mutau,
                          conj(nsiU_etau), conj(nsiU_mutau), nsiU_tautau;
            return NSIMatrixU;
   }

  //* NSI Matrix for d-quark *//
      Matrix<complex<double>,dim,dim> NeutrinoOscProbNSI :: SetNSIMatrixD()
   {
            Matrix<complex<double>,dim,dim> NSIMatrixD;
            NSIMatrixD << nsiD_ee,         nsiD_emu,         nsiD_etau,
                          conj(nsiD_emu),  nsiD_mumu,        nsiD_mutau,
                          conj(nsiD_etau), conj(nsiD_mutau), nsiD_tautau;
            return NSIMatrixD;
   }

  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//





