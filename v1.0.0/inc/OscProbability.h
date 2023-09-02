
    //////////////////////////////////////////////////////////////////
    //   .........oooOOO 18th-Aug-2020 Last Update OOOooo........   //
    //                        SADASHIV SAHOO                        //
    //               India-based Neutrino Observatory               //
    //        Homi Bhabha National Institute, Mumbai, INDIA         // 
    //          @email: sadashiv.sahoo@tifr.res.in                  //
    //                : sadashiv.sahoo@iopb.res.in                  //
    //////////////////////////////////////////////////////////////////

    #ifndef _OscProbability_H
    #define _OscProbability_H


    #include <vector>
    #include <complex>
    #include <iostream>

    #include <OscConstants.h>
    #include <eigen3/Eigen/Core>
    #include <eigen3/Eigen/Dense>
    #include <OscComplxMatrixExp.h>

    using namespace std;
    using namespace Eigen;

    #include <NeutrinoOscProbSTD.h>
    class NeutrinoOscProbSTD;

    #include <NeutrinoOscProbSME.h>
    class NeutrinoOscProbSME;

    #include <NeutrinoOscProbNSI.h>
    class NeutrinoOscProbNSI;

    #include <NeutrinoOscProbDEC.h>
    class NeutrinoOscProbDEC;

    class OscProbability : virtual public NeutrinoOscProbSTD,
                           virtual public NeutrinoOscProbSME,
                           virtual public NeutrinoOscProbNSI,
                           virtual public NeutrinoOscProbDEC
   { 
      private:
      double X1; // dummy

      protected:
      double X2; // dummy

      public:
      virtual ~OscProbability() {}; 
 
   };

   #endif







