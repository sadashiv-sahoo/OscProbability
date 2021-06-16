
//* v1.0.0                              *//
//* Test4.cxx                           *//
//* cd examples                         *//
//* cmake -DSET_TARGET=Test4 ../v1.0.0/ *//
//* CMake 3.20.3 version is required    *//

//* Start *//

  #include <new>
  #include <fstream>
  #include <iomanip>
  #include <bits/stdc++.h>
  #include <OscProbability.h>

  using namespace std;

  int main(void){

      setprecision(64);
      OscProbability *Y1 = new OscProbability();

      double AtmH = 15.0;
      double DptH = 0.00;

      double E;
      double Emin = 0.001; // 1  MeV (Minimum Energy)
      double Emax = 25.00; // 25 GeV (Maximum Energy)
      double dE   = 0.001; // 1  MeV (Resolution Energy)
      int Esteps  = int((Emax-Emin)/dE);

      //* Default Input Parameters *//
      double M21  = 7.40E-5; // Solar       (in eV)
      double M31  = 2.53E-3; // Atmospheric (in eV)

      double Theta12 = 33.808729;   // Solar      (in degree)
      double Theta23 = 45.000000;   // Atmosphere (in degree)
      double Theta13 = 8.6028772;   // Reactor    (in degree)
      double Dcp     = 0.0000000;   // Leptonic   (in degree)

      //* Default Settings of Parameters *//
      Y1 -> SetMassSqDiff(M21, M31);
      Y1 -> SetMixAngles(Theta12, Theta23, Theta13, Dcp);

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
//          cout << E 
//               << "\t" << MuE
//               << endl;
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



