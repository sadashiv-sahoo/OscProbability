/**
 *
 *  @file OscSetNSIFunctions.cxx
 *
 *  @brief It deals with Physics associated with NC-NSI Physics Functions.
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #include <OscProbability.h>
    #include <OscProbabilityNSI.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Setting SME-GILIV Parameters *//

      void OscProbabilityNSI ::
      SetNCNSI_E(SetBSMIntParams eNN[], size_t _eTnn_)
    {
        int ii = 0;
        BSMFlvchk(_eTnn_, dim, "NC-NSI-Electron");
        NCNSI_E.clear();
        for (size_t ij = 0; ij < _eTnn_; ij++)
      {
        Alignment(eNN[ij].F1, eNN[ij].F2, dim, "Flavour", "NC-NSI-Electron");
        if (eNN[ij].F1 > eNN[ij].F2) {ii++;}
        NCNSI_E.push_back({eNN[ij].Bsm, eNN[ij].F1, eNN[ij].F2});
      }
        _BSM_SelfAdjoint_check(ii, "NC-NSI-Electron");
    }

      void OscProbabilityNSI ::
      SetNCNSI_U(SetBSMIntParams uNN[], size_t _uTnn_)
    {
        int ii = 0;
        BSMFlvchk(_uTnn_, dim, "NC-NSI-UpQuark");
        NCNSI_U.clear();
        for (size_t ij = 0; ij < _uTnn_; ij++)
      {
        Alignment(uNN[ij].F1, uNN[ij].F2, dim, "Flavour", "NC-NSI-UpQuark");
        if (uNN[ij].F1 > uNN[ij].F2) {ii++;}
        NCNSI_U.push_back({uNN[ij].Bsm, uNN[ij].F1, uNN[ij].F2});
      }
        _BSM_SelfAdjoint_check(ii, "NC-NSI-UpQuark");
    }

      void OscProbabilityNSI ::
      SetNCNSI_D(SetBSMIntParams dNN[], size_t _dTnn_)
    {
        int ii = 0;
        BSMFlvchk(_dTnn_, dim, "NC-NSI-DownQuark");
        NCNSI_D.clear();
        for (size_t ij = 0; ij < _dTnn_; ij++)
      {
        Alignment(dNN[ij].F1, dNN[ij].F2, dim, "Flavour", "NC-NSI-DownQuark");
        if (dNN[ij].F1 > dNN[ij].F2) {ii++;}
        NCNSI_D.push_back({dNN[ij].Bsm, dNN[ij].F1, dNN[ij].F2});
      }
        _BSM_SelfAdjoint_check(ii, "NC-NSI-DownQuark");
    }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityNSI ::
      SetNSIMatrixE(int type, double _rho_)
    {
      const double _NSIInt_ = sqrt2GF*_ne_*type*_rho_;
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > NSIE(dim, dim);
      NSIE = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Zero(dim,dim);
      for (size_t ij = 0; ij < NCNSI_E.size(); ij ++) { NSIE(NCNSI_E[ij].F1, NCNSI_E[ij].F2) = NCNSI_E[ij].Bsm; }
      for (unsigned int i = 0; i < dim; i++) { for (unsigned int j = 0; j < dim; j++) { if (i > j) { NSIE(i,j) = conj(NSIE(j,i)); } } }
      if (type == 1) { return _NSIInt_*NSIE; }
      else           { return _NSIInt_*NSIE.conjugate(); }
    }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityNSI ::
      SetNSIMatrixU(int type, double _rho_)
    {
      const double _NSIInt_ = sqrt2GF*(2.*_np_ + _nn_)*type*_rho_;
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > NSIU(dim, dim);
      NSIU = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Zero(dim,dim);
      for (size_t ij = 0; ij < NCNSI_U.size(); ij ++) { NSIU(NCNSI_U[ij].F1, NCNSI_U[ij].F2) = NCNSI_U[ij].Bsm; }
      for (unsigned int i = 0; i < dim; i++) { for (unsigned int j = 0; j < dim; j++) { if (i > j) { NSIU(i,j) = conj(NSIU(j,i)); } } }
      if (type == 1) { return _NSIInt_*NSIU; }
      else           { return _NSIInt_*NSIU.conjugate(); }
    }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityNSI ::
      SetNSIMatrixD(int type, double _rho_)
    {
      const double _NSIInt_ = sqrt2GF*(_np_ + 2.0*_nn_)*type*_rho_;
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > NSID(dim, dim);
      NSID = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Zero(dim,dim);
      for (size_t ij = 0; ij < NCNSI_D.size(); ij ++) { NSID(NCNSI_D[ij].F1, NCNSI_D[ij].F2) = NCNSI_D[ij].Bsm; }
      for (unsigned int i = 0; i < dim; i++) { for (unsigned int j = 0; j < dim; j++) { if (i > j) { NSID(i,j) = conj(NSID(j,i)); } } }
      if (type == 1) { return _NSIInt_*NSID; }
      else           { return _NSIInt_*NSID.conjugate(); }
    }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityNSI ::
      SetNCNSIMatrix(int type, double _rho_) { return OscProbabilityNSI :: SetNSIMatrixE(type,_rho_) + OscProbabilityNSI :: SetNSIMatrixU(type,_rho_) + OscProbabilityNSI :: SetNSIMatrixD(type,_rho_); }

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//





