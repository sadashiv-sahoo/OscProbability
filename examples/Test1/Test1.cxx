
//* v1.0.0                              *//
//* Test1.cxx                           *//
//* cd examples                         *//
//* cmake -DSET_TARGET=Test1 ../v1.0.0/ *//
//* CMake 3.20.3 version is required    *//

//* Start *//

  #include <new>
  #include <fstream>
  #include <iomanip>
  #include <progress.h>
  #include <OscProbability.h>

  using namespace std;

    int main(void){

    setprecision(64);
    OscProbability *Y1 = new OscProbability();

    double E;
    double Emin = 1.E-3;  // (0.001 GeV)
    double Emax = 25.00;  // (25.00 GeV)
    double dE   = 1.E-3;  // (0.001 GeV)
    int Esteps  = int((Emax-Emin)/dE);

    double L    = 4500.;  // (4500.0 Km)

    //* Mandatory Input Parameters *//
    double M21  = 7.55E-5;  // Solar(in eV)
    double M31  = 2.50E-3;  // Atmosphere(in eV)

    double Th12 = 34.50;    // Solar      (in degree)
    double Th23 = 47.70;    // Atmosphere (in degree)
    double Th13 = 8.45;     // Reactor    (in degree)
    double Dcp  = 0.00;     // Leptonic   (in degree)

    //* Default Settings of Parameters *//
    Y1 -> SetMassSqDiff(M21, M31);
    Y1 -> SetMixAngles(Th12, Th23, Th13, Dcp);

     for(int j = 0; j <= Esteps; j++)
  { 
     E = Emin + dE*j;
     cout << E 
          << "\t" << Y1 -> VaccumOscProb(1, 0, 1, E, L)
          << endl;
  }
     delete Y1;
     return 0;
  }

// * End *//







