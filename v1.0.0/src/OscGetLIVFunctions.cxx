
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
  //* GetLIVFunctions *//

  //* Vacuum Modified Mixing Angles *//
      double NeutrinoOscProbSME :: 
                     GetVacGILIVMixAngle(int i, int j, int type, double Erg)
   {
        if(type != 0)
     {
          if((i == 1 || i == 2) && (j == 2 || j == 3))
       {
            if(i != j)
         {
            int ij = i;
            int jk = j;

            double Ue2; double Ue3; double UMu3;
            double ModTh12; double ModTh13; double ModTh23;
            double Sin2Th12; double Sin2Th13; double Sin2Th23;

            Matrix<complex<double>,dim,dim> U;
            Matrix<complex<double>,dim,dim> V;
            Matrix<complex<double>,dim,dim> H;
            Matrix<complex<double>,dim,dim> A;
            Matrix<complex<double>,dim,dim> C;
            Matrix<complex<double>,dim,dim> HF;
            int Type = int (type/(abs(type)));

            U    = NeutrinoOscProbSTD :: SetPMNS();
            A    = NeutrinoOscProbSME :: SetGILIVMatrixA(Type);
            C    = NeutrinoOscProbSME :: SetGILIVMatrixC(Type);
            H    = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);

            HF   = U*H*(U.inverse());
            HF  += (Type)*A;
            HF  -= (Erg*1E9)*C;

            ComplexEigenSolver <MatrixXcd> ces;
            ces.compute(HF);
            U    = ces.eigenvectors();

            Ue2  = abs(U(0,1));  // For Theta12
            Ue3  = abs(U(0,2));  // For Theta13
            UMu3 = abs(U(1,2));  // For Theta23

               if(ij == 2 && jk ==3){                                 //                       |UMu3|^2
                   Sin2Th23 = ((pow(UMu3,2)))/(1.0 - (pow(Ue3,2)));   //    Sin^2(Th23) =  ---------------
                   ModTh23  = asin(abs(sqrt(abs(Sin2Th23))));         //                    1.0 - |Ue3|^2
                 return (ModTh23)*deg;
             }
               if(ij == 1 && jk == 3){
                   Sin2Th13 = pow(Ue3,2);                             //    Sin^2(Th13) =      |Ue3|^2  
                   ModTh13  = asin(abs(sqrt(abs(Sin2Th13))));         //        Th13    = Sin^-1(Sqrt(Sin^2(Th13)))
                 return (ModTh13)*deg;
             }
               else {                                                 //                       |Ue2|^2
                   Sin2Th12 = ((pow(Ue2,2)))/(1.0 - (pow(Ue3,2)));    //    Sin^2(Th12) =   --------------
                   ModTh12  = asin(abs(sqrt(abs(Sin2Th12))));         //                     1.0 - |Ue3|^2
                 return (ModTh12)*deg;
             }
         }    
               else 
             { cout << "\n\033[0;91mError :: [i] can't be equal to [j] :: (2,2) not allowed \033[0m\n" 
                    << endl;
                    exit (-2);
             }
          }    else 
             { cout << "\n\033[0;91mError :: Inputs Must be (1,2) or (1,3) or (2,3) as per MixAngle indices \033[0m\n" 
                    << endl;
                   exit(-1);
             }
         }    else
            { cout << "\n\033[0;91mError :: Type should be integer and not equal to 0\033[0m\n" 
                   << endl;
                   exit(0);
            }
    }

  //* ConstMatter Modified Mixing Angles *//
      double NeutrinoOscProbSME :: 
                     GetConstMatGILIVMixAngle(int i, int j, int type, double dens, double Erg)
   {
        if(type != 0)
     {
          if((i == 1 || i == 2) && (j == 2 || j == 3))
       {
            if(i != j && dens >= 0.0)
         {
            int ij = i;
            int jk = j;

            double Ue2; double Ue3; double UMu3;                        
            double ModTh12; double ModTh13; double ModTh23;                   
            double Sin2Th12; double Sin2Th13; double Sin2Th23;                    

            Matrix<complex<double>,dim,dim> U;  
            Matrix<complex<double>,dim,dim> V;  
            Matrix<complex<double>,dim,dim> H;
            Matrix<complex<double>,dim,dim> F;
            Matrix<complex<double>,dim,dim> A;
            Matrix<complex<double>,dim,dim> C;
            Matrix<complex<double>,dim,dim> HF;
            int Type = int (type/(abs(type)));

            U    = NeutrinoOscProbSTD :: SetPMNS();
            A    = NeutrinoOscProbSME :: SetGILIVMatrixA(Type);
            C    = NeutrinoOscProbSME :: SetGILIVMatrixC(Type);
            H    = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);

            HF   = U*H*(U.inverse());
            HF  += (Type)*A;
            HF  -= (Erg*1E9)*C;

            ComplexEigenSolver <MatrixXcd> ces;
            typedef complex <double> v_;
            V  << v_(Type*sqrt2GF*dens*eS,0.), v_(0.,0.), v_(0.,0.),  
                                    v_(0.,0.), v_(0.,0.), v_(0.,0.),
                                    v_(0.,0.), v_(0.,0.), v_(0.,0.);
            F    = HF + V;
            ces.compute(F);
            U    = ces.eigenvectors();

            Ue2  = abs(U(0,1));  // For Theta12
            Ue3  = abs(U(0,2));  // For Theta13
            UMu3 = abs(U(1,2));  // For Theta23

               if(ij == 2 && jk ==3){                                 //                       |UMu3|^2
                   Sin2Th23 = ((pow(UMu3,2)))/(1.0 - (pow(Ue3,2)));   //    Sin^2(Th23) =  ---------------
                   ModTh23  = asin(abs(sqrt(abs(Sin2Th23))));         //                    1.0 - |Ue3|^2
                 return (ModTh23)*deg;
             }
               if(ij == 1 && jk == 3){
                   Sin2Th13 = pow(Ue3,2);                             //    Sin^2(Th13) =      |Ue3|^2  
                   ModTh13  = asin(abs(sqrt(abs(Sin2Th13))));         //        Th13    = Sin^-1(Sqrt(Sin^2(Th13)))
                 return (ModTh13)*deg;
             }
               else {                                                 //                       |Ue2|^2
                   Sin2Th12 = ((pow(Ue2,2)))/(1.0 - (pow(Ue3,2)));    //    Sin^2(Th12) =   --------------
                   ModTh12  = asin(abs(sqrt(abs(Sin2Th12))));         //                     1.0 - |Ue3|^2
                 return (ModTh12)*deg;
             }
         }    
               else 
             { cout << "\n\033[0;91mError :: [i] can't be equal to [j] :: (2,2) not allowed \033[0m\n" 
                    << endl;
                    exit (-2);
             }
          }    else 
             { cout << "\n\033[0;91mError :: Inputs Must be (1,2) or (1,3) or (2,3) as per MixAngle indices \033[0m\n" 
                    << endl;
                   exit(-1);
             }
         }    else
            { cout << "\n\033[0;91mError :: Type should be integer and not equal to 0\033[0m\n" 
                   << endl;
                   exit(0);
            }
    }

  //* Profile Modified Mixing Angles *//
      double NeutrinoOscProbSME ::
                     GetProfileGILIVMixAngle(int i, int j, int type, double Erg)
   {
        if(type != 0)
     {
          if((i == 1 || i == 2) && (j == 2 || j == 3))
       {
            if(i != j)
         {
            int ij = i;
            int jk = j;

            double Ue2; double Ue3; double UMu3;
            double ModTh12; double ModTh13; double ModTh23;
            double Sin2Th12; double Sin2Th13; double Sin2Th23;

            Matrix<complex<double>,dim,dim> U;
            Matrix<complex<double>,dim,dim> V;
            Matrix<complex<double>,dim,dim> H;
            Matrix<complex<double>,dim,dim> R;
            Matrix<complex<double>,dim,dim> F;
            Matrix<complex<double>,dim,dim> A;
            Matrix<complex<double>,dim,dim> C;
            Matrix<complex<double>,dim,dim> Um;
            Matrix<complex<double>,dim,dim> HF;
            vector < Matrix<complex<double>,dim,dim> > T;
            typedef complex <double> v_;
            int Type = int (type/(abs(type)));

            U    = NeutrinoOscProbSTD :: SetPMNS();
            A    = NeutrinoOscProbSME :: SetGILIVMatrixA(Type);
            C    = NeutrinoOscProbSME :: SetGILIVMatrixC(Type);
            H    = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);

            HF   = U*H*(U.inverse());
            HF  += (Type)*A;
            HF  -= (Erg*1E9)*C;

            T.push_back(U);
            ComplexEigenSolver <MatrixXcd> ces;
            for(int ii = 0; ii < int(Profile.size()); ii++)
         {
            V  << v_(Type*sqrt2GF*(Profile[ii].second)*eS,0.), v_(0.,0.), v_(0.,0.),
                                                    v_(0.,0.), v_(0.,0.), v_(0.,0.),
                                                    v_(0.,0.), v_(0.,0.), v_(0.,0.);
            F  = HF + V;
            ces.compute(F);
            R  = ces.eigenvectors();
            T.push_back(R);
         }
                R    = Matrix<complex<double>,dim,dim> :: Identity();
                for(int ij = int(T.size())-1; ij > 0; ij--) { R *= T[ij]*(T[ij-1].inverse()); }

                Um   = R*U;                                           // PMNS (Modified)
                Ue2  = abs(Um(0,1));
                Ue3  = abs(Um(0,2));
                UMu3 = abs(Um(1,2));

               if(ij == 2 && jk ==3){                                 //                       |UMu3|^2
                   Sin2Th23 = ((pow(UMu3,2)))/(1.0 - (pow(Ue3,2)));   //    Sin^2(Th23) =  ---------------
                   ModTh23  = asin(abs(sqrt(abs(Sin2Th23))));         //                    1.0 - |Ue3|^2
                 return (ModTh23)*deg;
             }
               if(ij == 1 && jk == 3){
                   Sin2Th13 = pow(Ue3,2);                             //    Sin^2(Th13) =      |Ue3|^2  
                   ModTh13  = asin(abs(sqrt(abs(Sin2Th13))));         //        Th13    = Sin^-1(Sqrt(Sin^2(Th13)))
                 return (ModTh13)*deg;
             }
               else {                                                 //                       |Ue2|^2
                   Sin2Th12 = ((pow(Ue2,2)))/(1.0 - (pow(Ue3,2)));    //    Sin^2(Th12) =   --------------
                   ModTh12  = asin(abs(sqrt(abs(Sin2Th12))));         //                     1.0 - |Ue3|^2
                 return (ModTh12)*deg;
             }
         }    
               else 
             { cout << "\n\033[0;91mError :: [i] can't be equal to [j] :: (2,2) not allowed \033[0m\n" 
                    << endl;
                    exit (-2);
             }
          }    else 
             { cout << "\n\033[0;91mError :: Inputs Must be (1,2) or (1,3) or (2,3) as per MixAngle indices \033[0m\n" 
                    << endl;
                   exit(-1);
             }
         }    else
            { cout << "\n\033[0;91mError :: Type should be integer and not equal to 0\033[0m\n" 
                   << endl;
                   exit(0);
            }
    }

  //* Earth Modified Mixing Angles *//
      double NeutrinoOscProbSME ::
                     GetEarthGILIVMixAngle(int i, int j, int type, double Erg)
   {
        if(type != 0)
     {
         if((i == 1 || i == 2) && (j == 2 || j == 3))
       {
            if(i != j)
         {
            int ij = i;
            int jk = j;

            double Ue2; double Ue3; double UMu3;
            double ModTh12; double ModTh13; double ModTh23;
            double Sin2Th12; double Sin2Th13; double Sin2Th23;

            Matrix<complex<double>,dim,dim> U;  
            Matrix<complex<double>,dim,dim> V;  
            Matrix<complex<double>,dim,dim> H;  
            Matrix<complex<double>,dim,dim> R;  
            Matrix<complex<double>,dim,dim> F;
            Matrix<complex<double>,dim,dim> A;
            Matrix<complex<double>,dim,dim> C;
            Matrix<complex<double>,dim,dim> Um;
            Matrix<complex<double>,dim,dim> HF;
            vector < Matrix<complex<double>,dim,dim> > T;
            typedef complex <double> v_;
            int Type = int (type/(abs(type)));

            U    = NeutrinoOscProbSTD :: SetPMNS();
            A    = NeutrinoOscProbSME :: SetGILIVMatrixA(Type);
            C    = NeutrinoOscProbSME :: SetGILIVMatrixC(Type);
            H    = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);

            HF   = U*H*(U.inverse());
            HF  += (Type)*A;
            HF  -= (Erg*1E9)*C;

            T.push_back(U);
            ComplexEigenSolver <MatrixXcd> ces;
            for(int ii = 0; ii < int(Out_In.size()); ii++)
         {
            V  << v_(Type*sqrt2GF*(Out_In[ii].second)*eS,0.), v_(0.,0.), v_(0.,0.),  
                                                   v_(0.,0.), v_(0.,0.), v_(0.,0.),
                                                   v_(0.,0.), v_(0.,0.), v_(0.,0.);
            F  = HF + V;
            ces.compute(F);
            R  = ces.eigenvectors();
            T.push_back(R);
         }
                R  = Matrix<complex<double>,dim,dim> :: Identity();
                for(int ij = int(T.size())-1; ij > 0; ij--) { R *= T[ij]*(T[ij-1].inverse()); }

                Um   = R*U;                                           // PMNS (Modified)
                Ue2  = abs(Um(0,1));
                Ue3  = abs(Um(0,2));
                UMu3 = abs(Um(1,2));

               if(ij == 2 && jk ==3){                                 //                       |UMu3|^2
                   Sin2Th23 = ((pow(UMu3,2)))/(1.0 - (pow(Ue3,2)));   //    Sin^2(Th23) =  ---------------
                   ModTh23  = asin(abs(sqrt(abs(Sin2Th23))));         //                    1.0 - |Ue3|^2
                 return (ModTh23)*deg;
             }
               if(ij == 1 && jk == 3){
                   Sin2Th13 = pow(Ue3,2);                             //    Sin^2(Th13) =      |Ue3|^2  
                   ModTh13  = asin(abs(sqrt(abs(Sin2Th13))));         //        Th13    = Sin^-1(Sqrt(Sin^2(Th13)))
                 return (ModTh13)*deg;
             }
               else {                                                 //                       |Ue2|^2
                   Sin2Th12 = ((pow(Ue2,2)))/(1.0 - (pow(Ue3,2)));    //    Sin^2(Th12) =   --------------
                   ModTh12  = asin(abs(sqrt(abs(Sin2Th12))));         //                     1.0 - |Ue3|^2
                 return (ModTh12)*deg;
             }
         }    
               else 
             { cout << "\n\033[0;91mError :: [i] can't be equal to [j] :: (2,2) not allowed \033[0m\n" 
                    << endl;
                    exit (-2);
             }
          }    else 
             { cout << "\n\033[0;91mError :: Inputs Must be (1,2) or (1,3) or (2,3) as per MixAngle indices \033[0m\n" 
                    << endl;
                   exit(-1);
             }
         }    else
            { cout << "\n\033[0;91mError :: Type should be integer and not equal to 0\033[0m\n" 
                   << endl;
                   exit(0);
            }
    }

  //* Vacuum Modified Mass *//
      double NeutrinoOscProbSME :: 
                     GetVacGILIVModMass(int i, int j, int type, double Erg)
   {
        if(i > 0 && i < 4 && j > 0 && j < 4)
     {
          if(type != 0)
       {
           double MMM;
           Matrix<complex<double>,dim,dim> U;  
           Matrix<complex<double>,dim,dim> V;  
           Matrix<complex<double>,dim,dim> H;
           Matrix<complex<double>,dim,dim> A;
           Matrix<complex<double>,dim,dim> C;  
           Matrix<complex<double>,dim,dim> M;  
           Matrix<complex<double>,dim,dim> HF;
           int Type = int (type/(abs(type)));

           U   = NeutrinoOscProbSTD :: SetPMNS();
           A   = NeutrinoOscProbSME :: SetGILIVMatrixA(Type);
           C   = NeutrinoOscProbSME :: SetGILIVMatrixC(Type);
           H   = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);

           HF  = U*H*(U.inverse());
           HF += (Type)*A;
           HF -= (Erg*1E9)*C;

           ComplexEigenSolver <MatrixXcd> ces;
           ces.compute(HF);
           M.diagonal() << (ces.eigenvalues())(0,0), //
                           (ces.eigenvalues())(1,0), //
                           (ces.eigenvalues())(2,0); //
           M  -= ((1.0/dim)*(M.trace()))*(Matrix<complex<double>,dim,dim> :: Identity());
           MMM = 2.0*(Erg*1E9)*((M(i-1,i-1).real())-(M(j-1,j-1).real()));
           return MMM;
       }
           else
         { cout << "\n\033[0;91mError :: Type should be integer and not equal to 0\033[0m\n"
                << endl;
               exit(-1);
         }
      }
          else 
        { cout << "\n\033[0;91mError :: Mass Basis ranges [M21; M31; M32] \033[0m\n"
               << endl;
               exit(0);
        }
   }

  //* ConstMatter Modified Mass *//
      double NeutrinoOscProbSME :: 
                     GetConstMatGILIVModMass(int i, int j, int type, double dens, double Erg)
   {
        if(i > 0 && i < 4 && j > 0 && j < 4)
     {
          if(type != 0 && dens >= 0.0)
       {
           double MMM;
           Matrix<complex<double>,dim,dim> U;
           Matrix<complex<double>,dim,dim> V;
           Matrix<complex<double>,dim,dim> H;
           Matrix<complex<double>,dim,dim> F;
           Matrix<complex<double>,dim,dim> A;
           Matrix<complex<double>,dim,dim> C;
           Matrix<complex<double>,dim,dim> M;
           Matrix<complex<double>,dim,dim> HF;
           typedef complex <double> v_;
           int Type = int (type/(abs(type)));

           U   = NeutrinoOscProbSTD :: SetPMNS();
           A   = NeutrinoOscProbSME :: SetGILIVMatrixA(Type);
           C   = NeutrinoOscProbSME :: SetGILIVMatrixC(Type);
           H   = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);

           HF  = U*H*(U.inverse());
           HF += (Type)*A;
           HF -= (Erg*1E9)*C;

           ComplexEigenSolver <MatrixXcd> ces;
           V  << v_(Type*sqrt2GF*dens*eS,0.), v_(0.,0.), v_(0.,0.),
                                   v_(0.,0.), v_(0.,0.), v_(0.,0.),
                                   v_(0.,0.), v_(0.,0.), v_(0.,0.);
           F   = HF + V;
           ces.compute(F);
           M.diagonal() << (ces.eigenvalues())(0,0), //
                           (ces.eigenvalues())(1,0), //
                           (ces.eigenvalues())(2,0); //
           M  -= ((1.0/dim)*(M.trace()))*(Matrix<complex<double>,dim,dim> :: Identity());
           MMM = 2.0*(Erg*1E9)*((M(i-1,i-1).real())-(M(j-1,j-1).real()));
           return MMM;
       }
           else
         { cout << "\n\033[0;91mError :: Type should be integer and not equal to 0\033[0m\n"
                << endl;
               exit(-1);
         }
      }
          else
        { cout << "\n\033[0;91mError :: Mass Basis ranges [M21; M31; M32] \033[0m\n"
               << endl;
               exit(0);
        }
   }

  //* Profiled Modified Mass *//
      double NeutrinoOscProbSME ::
                     GetProfileGILIVModMass(int i, int j, int type, double Erg)
   {
        if(i > 0 && i < 4 && j > 0 && j < 4)
     {
          if(type != 0)
       {
           double MMM;
           Matrix<complex<double>,dim,dim> U;
           Matrix<complex<double>,dim,dim> V;
           Matrix<complex<double>,dim,dim> H;
           Matrix<complex<double>,dim,dim> F;
           Matrix<complex<double>,dim,dim> A;
           Matrix<complex<double>,dim,dim> C;
           Matrix<complex<double>,dim,dim> M;
           Matrix<complex<double>,dim,dim> HF;
           vector < Matrix<complex<double>,dim,dim> > T;
           typedef complex <double> v_;
           int Type = int (type/(abs(type)));

            U    = NeutrinoOscProbSTD :: SetPMNS();
            A    = NeutrinoOscProbSME :: SetGILIVMatrixA(Type);
            C    = NeutrinoOscProbSME :: SetGILIVMatrixC(Type);
            H    = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);

            HF   = U*H*(U.inverse());
            HF  += (Type)*A;
            HF  -= (Erg*1E9)*C;

            T.push_back(H);
            ComplexEigenSolver <MatrixXcd> ces;
            for(int ii = 0; ii < int(Profile.size()); ii++)
         {
             V  << v_(Type*sqrt2GF*(Profile[ii].second)*eS,0.), v_(0.,0.), v_(0.,0.),
                                                     v_(0.,0.), v_(0.,0.), v_(0.,0.),
                                                     v_(0.,0.), v_(0.,0.), v_(0.,0.);
             F   = U*H*(U.inverse()) + V;
             ces.compute(F);
             M.diagonal() << (ces.eigenvalues())(0,0), // 
                             (ces.eigenvalues())(1,0), // 
                             (ces.eigenvalues())(2,0); //
             T.push_back(M);
         }
           for(int ij = int(T.size())-1; ij > 0; ij--) { H += T[ij]-T[ij-1]; }
           H  -= ((1.0/dim)*(H.trace()))*(Matrix<complex<double>,dim,dim> :: Identity());
           MMM = 2.0*(Erg*1E9)*((H(i-1,i-1).real())-(H(j-1,j-1).real()));
           return MMM;
       }
           else
         { cout << "\n\033[0;91mError :: Type should be integer and not equal to 0\033[0m\n"
                << endl;
               exit(-1);
         }
      }
          else
        { cout << "\n\033[0;91mError :: Mass Basis ranges [M21; M31; M32] \033[0m\n"
               << endl;
               exit(0);
        }
   }

  //* Earth Modified Mass *//
      double NeutrinoOscProbSME :: 
                     GetEarthGILIVModMass(int i, int j, int type, double Erg)
   {
        if(i > 0 && i < 4 && j > 0 && j < 4)
     {
          if(type != 0)
       {
           double MMM;
           Matrix<complex<double>,dim,dim> U;  
           Matrix<complex<double>,dim,dim> V;  
           Matrix<complex<double>,dim,dim> H;
           Matrix<complex<double>,dim,dim> F;
           Matrix<complex<double>,dim,dim> A;
           Matrix<complex<double>,dim,dim> C;  
           Matrix<complex<double>,dim,dim> M;  
           Matrix<complex<double>,dim,dim> HF;
           vector < Matrix<complex<double>,dim,dim> > T;
           typedef complex <double> v_;
           int Type = int (type/(abs(type)));

            U    = NeutrinoOscProbSTD :: SetPMNS();
            A    = NeutrinoOscProbSME :: SetGILIVMatrixA(Type);
            C    = NeutrinoOscProbSME :: SetGILIVMatrixC(Type);
            H    = NeutrinoOscProbSTD :: SetMassHamiltonian(Erg);

            HF   = U*H*(U.inverse());
            HF  += (Type)*A;
            HF  -= (Erg*1E9)*C;

            T.push_back(H);
            ComplexEigenSolver <MatrixXcd> ces;
            for(int ii = 0; ii < int(Out_In.size()); ii++)
         {
             V  << v_(Type*sqrt2GF*(Out_In[ii].second)*eS,0.), v_(0.,0.), v_(0.,0.),  
                                                    v_(0.,0.), v_(0.,0.), v_(0.,0.),
                                                    v_(0.,0.), v_(0.,0.), v_(0.,0.);
             HF = U*H*(U.inverse()) + V;
             ces.compute(HF);
             M.diagonal() << (ces.eigenvalues())(0,0), // 
                             (ces.eigenvalues())(1,0), // 
                             (ces.eigenvalues())(2,0); //
             T.push_back(M);
         }
           for(int ij = int(T.size())-1; ij > 0; ij--) { H += T[ij]-T[ij-1]; }
           H  -= ((1.0/dim)*(H.trace()))*(Matrix<complex<double>,dim,dim> :: Identity());
           MMM = 2.0*(Erg*1E9)*((H(i-1,i-1).real())-(H(j-1,j-1).real()));
           return MMM;
       }
           else
         { cout << "\n\033[0;91mError :: Type should be integer and not equal to 0\033[0m\n"
                << endl;
               exit(-1);
         }
      }   
          else
        { cout << "\n\033[0;91mError :: Mass Basis ranges [M21; M31; M32] \033[0m\n"
               << endl;
               exit(0);
        }
   }

  //* GVLIV GetHamiltonian *//
      Matrix<complex<double>,dim,dim> NeutrinoOscProbSME :: 
                     GetGVLIVHamiltonian(int ij, int jk, double Erg)
   {
        Matrix<complex<double>,dim,dim> H;
        H  = (NeutrinoOscProbSME :: SetGVLIVMatrixG())*Erg;                   // Hamiltonian = G;
        if(ij > 0 && jk < 0){ H -= NeutrinoOscProbSME :: SetGVLIVMatrixH(); } // Hamiltonian = G - H; [nu -> nubar]
        else                { H += NeutrinoOscProbSME :: SetGVLIVMatrixH(); } // Hamiltonian = G + H; [nubar -> nu]
        H *= (NeutrinoOscProbSME :: SetGVLIVMatrixC())*epsilon;               // Hamiltonian = [G - H]*[C] || [G + H]*[C];
        return H;
   }

  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//





