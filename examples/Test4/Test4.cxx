
//* v3.0.0                              *//
//* Test4.cxx                           *//
//* cd examples                         *//
//* cmake -DSET_TARGET=Test4 ../v3.0.0/ *//
//* CMake 3.20.3 version is required    *//

//* Start *//

  #include <new>
  #include <fstream>
  #include <iomanip>
  #include <bits/stdc++.h>
  #include <OscProbability.h>

  using namespace std;

    int main(void)
  {
    setprecision(64);
    const int n = 3;
    OscProbability *Y1  = new OscProbability(n);

    double AtmH = 15.0;
    double DptH = 0.00;

    double E;
    double Emin = 0.001; // 1  MeV (Minimum Energy)
    double Emax = 25.00; // 25 GeV (Maximum Energy)
    double dE   = 0.001; // 1  MeV (Resolution Energy)
    int Esteps  = int((Emax-Emin)/dE);

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

    //* Mandatory Settings for EarthMatOscProb() *//
    Y1 -> SetEarthDensityProfile();
    Y1 -> SetDetDepth(DptH);
    Y1 -> SetEarthAtmos(AtmH);
    Y1 -> SetZenith(-1.0);

    int j;
    double MuE;
    vector < pair < double, double > > P;
    #pragma omp parallel for private(E, MuE)  // Parallel Open-MP
         for(j = 0; j <= Esteps; j++)
      {         
         E   = Emin + dE*j;
         MuE = Y1 -> EarthMatOscProb(1, 0, 1, E);
         #pragma omp critical                 // Parallel Open-MP
       {
//       cout << E 
//            << "\t" << MuE
//            << endl;
         P.push_back(make_pair(E, MuE));
       }
     }

     sort(P.begin(), P.end());
     for(int i = 0; i < int(P.size()); i++)
   {
              cout << P[i].first 
                   << "\t" << P[i].second
                   << endl;
      }

   delete Y1;
     return 0;

  }

//* End *//



