
//* v3.0.0                              *//
//* Test2.cxx                           *//
//* cd examples                         *//
//* cmake -DSET_TARGET=Test2 ../v3.0.0/ *//
//* CMake 3.20.3 version is required    *//

//* Start *//

  #include <new>
  #include <fstream>
  #include <iomanip>
  #include <OscProbability.h>

  using namespace std;

    int main(void){

    setprecision(64);
    const int n = 3;
    OscProbability *Y1  = new OscProbability(n);

    double E;
    double Emin = 1.E-3;   // (0.001 GeV)
    double Emax = 25.00;   // (25.00 GeV)
    double dE   = 1.E-3;   // (0.001 GeV)
    int Esteps  = int((Emax-Emin)/dE);

    double L    = 4500.;  // (4500.0 Km)
    double rho  = 2.500;  // (2.50 g/cc)

    //* Mandatory Input Parameters *//
    double M21  = 7.55E-5; // Solar       (in eV^2)
    double M31  = 2.50E-3; // Atmosphere  (in eV^2)

    double Th12 = 34.00;   // Solar       (in degree)
    double Th23 = 45.00;   // Atmosphere  (in degree)
    double Th13 = 08.50;   // Reactor     (in degree)
    double Dcp  = 45.00;   // Leptonic CP (in degree)

    SetFlavoursPost F[] = {{"e", 0}, {"mu", 1}, {"tau", 2}};
    SetMixAngParams R[] = {{Th23, 0.0, 2, 3},
                           {Th13, Dcp, 1, 3},
                           {Th12, 0.0, 1, 2}};
    SetMassSqParams M[] = {{M21}, {M31}};

    //* Default Settings of Parameters *//
    Y1 -> SetMixAngles  (R, SizeOf(R));
    Y1 -> SetMassSqDiff (M, SizeOf(M));
    Y1 -> SetOscFlavours(F, SizeOf(F));

     for(int j = 0; j <= Esteps; j++)
  { 
     E = Emin + dE*j;
     cout << E 
          << "\t" << Y1 -> ConstMatOscProb(1, 0, 1, E, rho, L)
          << endl;
  }
     delete Y1;
     return 0;
  }

// * End *//







