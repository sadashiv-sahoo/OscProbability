/**
 *
 *  @file OscSetLIVFunctions.cxx
 *
 *  @brief It deals with Physics associated with LIV Physics Functions.
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #include <OscProbability.h>
    #include <OscProbabilityLIV.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Setting SME-GILIV Parameters *//

    //* GILIV Controls over Parameters *//
      void OscProbabilityLIV ::
      SetLIVSpTControls(int _mu_, int _nu_)
    {
        if (_mu_ >= 0 && _mu_ < 4 && _nu_ >= 0 && _nu_ < 4)
      {
        int i = 1; // (initiation);
        int j = 1; // (initiation);

        if (_nu_ == 0){j = 1;} else {j = -1;}
        if (_mu_ == 0){i = 1; aL = 1.0; _h_ = 1.0;} else {i = -1; aL = -1.0; _h_ = 1.0;}
        if (_mu_ == 0 && _nu_ == 0){cL = 4.0/3.0;} else {cL = 1.0;}
        cL *= i*j;
        _g_ = i*j;
      }
        else
      {
        std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nUnphysical Value Entered !!!\n\033[0m\n"
                  << "\033[0;96mFour-Momenta Index :: \033[0;91m(\033[0m" << _mu_ << "," << _nu_ << "\033[0;91m)\033[0m"
                  << std::endl;
        std::cout << "\033[0;96mPlease Enter :: Four-Momenta Index (0,1,2,3).\033[0m\n" 
                  << std::endl;
        exit(0);
      }
    }

      void OscProbabilityLIV ::
      SetGILIVaX(SetBSMIntParams aX[], size_t _av_)
    {
        int ii = 0;
        BSMFlvchk(_av_, dim, "CPTV-GILIV");
        GILIVCPTODD.clear();
        for (size_t ij = 0; ij < _av_; ij++)
      {
        Alignment(aX[ij].F1, aX[ij].F2, dim, "Flavour", "CPTV-GILIV");
        if (aX[ij].F1 > aX[ij].F2) {ii++;}
        GILIVCPTODD.push_back({aX[ij].Bsm, aX[ij].F1, aX[ij].F2});
      }
        _BSM_SelfAdjoint_check(ii, "CPTV-GILIV");
    }

      void OscProbabilityLIV ::
      SetGILIVcXX(SetBSMIntParams cXX[], size_t _cv_)
    {
        int ii = 0;
        BSMFlvchk(_cv_, dim, "CPTC-GILIV");
        GILIVCPTEVEN.clear();
        for (size_t ij = 0; ij < _cv_; ij++)
      {
        Alignment(cXX[ij].F1, cXX[ij].F2, dim, "Flavour", "CPTC-GILIV");
        if (cXX[ij].F1 > cXX[ij].F2) {ii++;}
        GILIVCPTEVEN.push_back({cXX[ij].Bsm, cXX[ij].F1, cXX[ij].F2});
      }
        _BSM_SelfAdjoint_check(ii, "CPTC-GILIV");
    }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityLIV ::
      SetGILIVMatrixA(int type, double _units_)
    {
      _LIV_Control_Check_(aL);
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > CPTV(dim, dim);
      CPTV = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Zero(dim,dim);
      for (size_t ij = 0; ij < GILIVCPTODD.size(); ij ++) { CPTV(GILIVCPTODD[ij].F1, GILIVCPTODD[ij].F2) = GILIVCPTODD[ij].Bsm; }
      for (unsigned int i = 0; i < dim; i++) { for (unsigned int j = 0; j < dim; j++) { if (i > j) { CPTV(i,j) = conj(CPTV(j,i)); } } }
      if (type == 1) { return aL*type*_units_*CPTV; }
      else           { return aL*type*_units_*CPTV.conjugate(); }
    }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityLIV ::
      SetGILIVMatrixC(int type, double _units_)
    {
      _LIV_Control_Check_(cL);
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > CPTC(dim, dim);
      CPTC = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Zero(dim,dim);
      for (size_t ij = 0; ij < GILIVCPTEVEN.size(); ij ++) { CPTC(GILIVCPTEVEN[ij].F1, GILIVCPTEVEN[ij].F2) = GILIVCPTEVEN[ij].Bsm; }
      for (unsigned int i = 0; i < dim; i++) { for (unsigned int j = 0; j < dim; j++) { if (i > j) { CPTC(i,j) = conj(CPTC(j,i)); } } }
      if (type == 1) { return cL*_units_*CPTC; }
      else           { return cL*_units_*CPTC.conjugate(); }
    }

    //*  Setting SME-GVLIV Parameters  *//
    //* GVLIV Controls over Parameters *//
      void OscProbabilityLIV ::
      SetGVLIVPolzation(std::complex < double > eps)
    {
        if (abs(eps) != 0.0) {_eps_ = eps/abs(eps);}
        else
      {
        std::cout << "\033[0;91mVacuum Polarization is necessary for Gauge Violating Oscillations\033[0m"
                  << std::endl;
        exit(0);
      }
    }

      void OscProbabilityLIV ::
      SetGILIVhXX(SetBSMIntParams hXX[], size_t _hv_)
    {
        int ii = 0;
        BSMFlvchk(_hv_, dim, "CPTV-GVLIV");
        GVLIVCPTODD.clear();
        for (size_t ij = 0; ij < _hv_; ij++)
      {
        Alignment(hXX[ij].F1, hXX[ij].F2, dim, "Flavour", "CPTV-GVLIV");
        if (hXX[ij].F1 > hXX[ij].F2) {ii++;}
        GVLIVCPTODD.push_back({hXX[ij].Bsm, hXX[ij].F1, hXX[ij].F2});
      }
        _BSM_SelfAdjoint_check(ii, "CPTV-GVLIV");
    }

      void OscProbabilityLIV ::
      SetGILIVgXXX(SetBSMIntParams gXXX[], size_t _gv_)
    {
        int ii = 0;
        BSMFlvchk(_gv_, dim, "CPTC-GVLIV");
        GVLIVCPTEVEN.clear();
        for (size_t ij = 0; ij < _gv_; ij++)
      {
        Alignment(gXXX[ij].F1, gXXX[ij].F2, dim, "Flavour", "CPTC-GVLIV");
        if (gXXX[ij].F1 > gXXX[ij].F2) {ii++;}
        GVLIVCPTEVEN.push_back({gXXX[ij].Bsm, gXXX[ij].F1, gXXX[ij].F2});
      }
        _BSM_SelfAdjoint_check(ii, "CPTC-GVLIV");
    }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityLIV ::
      SetGVLIVMatrixG()
    {
      _LIV_Control_Check_(_g_);
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > CPTV(dim, dim);
      CPTV = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Zero(dim,dim);
      for (size_t ij = 0; ij < GVLIVCPTODD.size(); ij ++) { CPTV(GVLIVCPTODD[ij].F1, GVLIVCPTODD[ij].F2) = GVLIVCPTODD[ij].Bsm; }
      for (unsigned int i = 0; i < dim; i++) { for (unsigned int j = 0; j < dim; j++) { if (i > j) { CPTV(i,j) = conj(CPTV(j,i)); } } }
      return _g_*CPTV;
    }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityLIV ::
      SetGVLIVMatrixH()
    {
      _LIV_Control_Check_(_h_);
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > CPTC(dim, dim);
      CPTC = Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Zero(dim,dim);
      for (size_t ij = 0; ij < GVLIVCPTEVEN.size(); ij ++) { CPTC(GVLIVCPTEVEN[ij].F1, GVLIVCPTEVEN[ij].F2) = GVLIVCPTEVEN[ij].Bsm; }
      for (unsigned int i = 0; i < dim; i++) { for (unsigned int j = 0; j < dim; j++) { if (i > j) { CPTC(i,j) = conj(CPTC(j,i)); } } }
      return _h_*CPTC;
    }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityLIV ::
      SetGVLIVMatrixC() { return Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > :: Identity(dim,dim); }

      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > OscProbabilityLIV ::
      SetGVLIVHamiltonian(int _Init_, double Erg)
    {
        Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > H(dim, dim);
        H = OscProbabilityLIV :: SetGVLIVMatrixG()*Erg*GeV;
        if      (_Init_ == -1) { H += OscProbabilityLIV :: SetGVLIVMatrixH(); } /// Hamiltonian = G + H; [nubar -> nu];
        else if (_Init_ == +1) { H -= OscProbabilityLIV :: SetGVLIVMatrixH(); } /// Hamiltonian = G - H; [nu -> nubar];
        else
      {
        std::cout << "\033[0;91mInitital Flavour must be either +1 or -1 only!!!\033[0m"
                  << std::endl;
        exit(0);
      }
        return H *= OscProbabilityLIV :: SetGVLIVMatrixC()*_eps_;
    }

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//






