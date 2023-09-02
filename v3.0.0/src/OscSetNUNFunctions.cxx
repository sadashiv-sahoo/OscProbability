/**
 *
 *  @file OscSetNUNFunctions.cxx
 *
 *  @brief It deals with Physics associated with Non-Unitarity Physics Functions.
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

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Setting Non-Unitarity Parameters *//

      void OscProbabilityNUN ::
      SetNonUnitarityPhys(bool _NsTs_) { _NonUni_Status_ = _NsTs_; }

      void OscProbabilityNUN ::
      SetNonUnitarityNorms(bool _NsNs_) { _NonNrm_Status_ = _NsNs_; }

      void OscProbabilityNUN ::
      SetNonUnitControl(SetBSMIntParams nUn[], size_t _nun_, double _sc_)
    {
        May_I_Proceed(_NonUni_Status_, "SetNonUnitarityPhys");
        int ii = 0;
        _Formulation_ = _sc_;
        BSMFlvchk(_nun_, dim, "elements in SetNonUnitControl");
        _NON_UNITARY_INFO_.clear();
        for (size_t ij = 0; ij < _nun_; ij++)
      {
        Alignment(nUn[ij].F1, nUn[ij].F2, dim, "Non-Unitarity", "SetNonUnitControl");
        if (nUn[ij].F1 < nUn[ij].F2) {ii++;}    /// Lower Triangular Matrix;
        _NON_UNITARY_INFO_.push_back({nUn[ij].Bsm, nUn[ij].F1, nUn[ij].F2});
      }
        _BSM_SelfAdjoint_check(ii, "SetNonUnitControl");
    }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityNUN ::
      SetNonUnitaryPMNS(int _type_)
    {
      _SYZ_VECTOR_CHECK_(_NON_UNITARY_INFO_.size(), "SetNonUnitControl(SetBSMIntParams , size_t , double )");
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > Alp(dim,dim);
      Alp = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Zero(dim,dim);
      for (size_t ij = 0; ij < _NON_UNITARY_INFO_.size(); ij ++) { Alp(_NON_UNITARY_INFO_[ij].F1-1, _NON_UNITARY_INFO_[ij].F2-1) = _NON_UNITARY_INFO_[ij].Bsm; }
      return (Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim) + _Formulation_*Alp)*OscProbabilitySTD :: SetPMNS(_type_);
    }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityNUN ::
      ZeroLengthNorm(const Eigen::Ref < const Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > > &NuN)
    {
      if   (_NonNrm_Status_ == true) { return NuN*NuN.adjoint(); }
      else                           { return Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim); }
    }

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//





