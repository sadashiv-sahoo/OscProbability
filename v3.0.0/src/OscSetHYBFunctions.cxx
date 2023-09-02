/**
 *
 *  @file OscSetHYBFunctions.cxx
 *
 *  @brief It deals with Physics associated with Hybrid Physics Functions.
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #include <OscProbability.h>
    #include <OscProbabilityHYB.h>

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//
    //* Setting Hybrid Parameters *//

      void OscProbabilityHYB ::
      May_I_Come_In()
    {
        if (_HybPhy_Status_ == false && _InElst_Status_ == false && _UsrInt_Status_ == false)
      {
        std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nNo Hybrid-physics Found !!!\n\033[0m\n" 
                  << std::endl;
        std::cout << "\033[0;96mPlease Enter any Hybrid Physics Conditions.\033[0m\n"
                  << "\033[0;96mTurn ON, either \"SetHybridPhyStatus(true)\" or \"SetFwdInElasticSts(true)\" or \"SetGBLUsersIntrStatus(true)\".\033[0m\n"
                  << std::endl;
        exit(-1);
      }
        else if (_HybPhy_Status_ == true  && _InElst_Status_ == false && _UsrInt_Status_ == false) {_MoA_ = 1;}
        else if (_HybPhy_Status_ == false && _InElst_Status_ == true  && _UsrInt_Status_ == false) {_MoA_ = 2;}
        else if (_HybPhy_Status_ == false && _InElst_Status_ == false && _UsrInt_Status_ == true ) {_MoA_ = 3;}
        else if (_HybPhy_Status_ == true  && _InElst_Status_ == true  && _UsrInt_Status_ == false) {_MoA_ = 4;}
        else if (_HybPhy_Status_ == false && _InElst_Status_ == true  && _UsrInt_Status_ == true ) {_MoA_ = 5;}
        else if (_HybPhy_Status_ == true  && _InElst_Status_ == false && _UsrInt_Status_ == true ) {_MoA_ = 6;}
        else {_MoA_ = 7;}
    }

      void OscProbabilityHYB ::
      SetHybridPhyStatus(bool _hsts_) { _HybPhy_Status_ = _hsts_; }

      void OscProbabilityHYB ::
      SetFwdInElasticSts(bool _fsts_) { _InElst_Status_ = _fsts_; }

      void OscProbabilityHYB ::
      SetGBLUsersIntrStatus(bool _usts_) { _UsrInt_Status_ = _usts_; }

      void OscProbabilityHYB ::
      SetMatterInElasticPhases(double _Pcc_, double _Pnc_) { _Ph_CC_ = _Pcc_*rad; _Ph_NC_ = _Pnc_*rad; }

      void OscProbabilityHYB ::
      SetGBLUserInteractions(SetBSMIntParams usr[], size_t _usr_, const char* c)
    {
         May_I_Proceed(_UsrInt_Status_, "SetGBLUsersIntrStatus");
         const size_t Zd = dim*dim;
         if (_usr_ <= Zd && _usr_ > 0)
      {
         const char* l1 = "Flavour";
         const char* l2 = "Mass";
         int chk1 = strcmp(l1, c);
         int chk2 = strcmp(l2, c);
         _GBLUSER_INTERACTIONS_.clear();
          if      (chk1 == 0 && chk2 != 0) { basis = 111; } 
          else if (chk1 != 0 && chk2 == 0) { basis = 999; } 
          else
        {
          std::cout << "\n\033[0;91m:: Error Warning ::\nCan't Proceed,\nWrong Instruction Found !!!\n\033[0m"
                    << std::endl;
          std::cout << "\033[0;96mPlease try with either \"\033[0mMass\033[0;96m\" or \"\033[0mFlavour\033[0;96m\".\033[0m\n"
                    << std::endl;
          exit(0);
        }
          for (size_t ij = 0; ij < _usr_; ij++)
        {
          Alignment(usr[ij].F1, usr[ij].F2, dim, c, "SetGBLUserInteractions");
          _GBLUSER_INTERACTIONS_.push_back({usr[ij].Bsm, usr[ij].F1, usr[ij].F2});
        }
      }
        else 
      {
        std::cout << "\n\033[0;91m:: Error Warning ::\nIn SetGBLUserInteractions, with entered No. of element(s) " << _usr_ << ", Can't Proceed !!!\033[0m\n"
                  << std::endl;
        std::cout << "\033[0;96mNo. of Elements must be > 0 and <= \" \033[0m" << Zd << "\033[0;96m \".\033[0m\n" 
                  << std::endl;
        exit(-1);
      }
    }

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//





