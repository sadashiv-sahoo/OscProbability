
  //////////////////////////////////////////////////////////////////
  //   .........oooOOO 18th-Aug-2020 Last Update OOOooo........   //
  //                        SADASHIV SAHOO                        //
  //               India-based Neutrino Observatory               //
  //        Homi Bhabha National Institute, Mumbai, INDIA         // 
  //          @email: sadashiv.sahoo@tifr.res.in                  //
  //                : sadashiv.sahoo@iopb.res.in                  //
  //////////////////////////////////////////////////////////////////

  //* Please Don't Change Written Format *//
  //* All Settings at Editor Tab Width:8 *//

  #include <ctime>
  #include <vector>
  #include <string>
  #include <fstream>
  #include <iostream>
  #include <algorithm>

  using namespace std;

     string PrintTime(char* a, int size)
  { 
     string s = ""; 
     for (int i = 0; i < size; i++) { s = s + a[i]; }
     return s; 
  }

   void getData(ifstream &f1, vector < pair < double, double > > &T)
 {
     double r; double rho;
     while(f1 >> r >> rho) 
   {
     if(r >= 0.00 && rho >= 0.00) { T.push_back(make_pair(r, rho)); } 
     else {cout << "\033[0;91mError :: Unphysical Values are entered in PREM-Profile\033[0m" << endl; exit (-1); }
   }
 }
     int main(void)
 {
     cout << "\n\033[0;93mThe PREM Profile Should Contain Increasing Order Of Layer-Radius!!!\n"
          << "Layer Radii Must be in KM Unit,"
          << " & Densities Must be in g/cc Unit\033[0m\n" << endl;

       string file;
       cout << "\033[0;96mEnter PREM source file : \033[0m";
       cin >> file;
       ifstream PREM_FILE(file);
       vector < pair < double, double > > PREM;

       if(PREM_FILE.good())
     {
       cout << "\033[0;96mAnalysing Entered PREM source File : \033[0m" << file << endl;
       PREM.clear();
       getData(PREM_FILE, PREM);
       cout << "\033[0;93m\n";
       for(int ij = 0; ij < int(PREM.size()); ij ++) { cout << ij + 1 << "\t" << PREM[ij].first << "\t" << PREM[ij].second << endl; }
       cout << "\033[0m\n";
     } else { cout << "\033[0;91mError Warning :: Can't Open Source File (\033[0m" 
                   << file << "\033[0;91m)\nPlease Check!!!\033[0m\n" 
                   << endl;
              exit(0);
            }

         //* Necessary Condition *//
         if(is_sorted(PREM.begin(), PREM.end()))
      {
         cout << "\n\033[0;92mEarth Layer Radii Are In Required Order.\nInput Data Is Accepted\033[0m" << endl;

         //* Don't Change File Directory *//
         ofstream outfile0("../src/PREM.info");
         ofstream outfile1("../src/OscEarthPREMProfle.cxx");

         const static int N = int(PREM.size());
         cout << "\033[0;96mNo. of PREM Layers :\t" << N << "\033[0m\n"<< endl;

         //* Day Month Date Hour:Minutes:Seconds Year *//
         int  tL     = 24;
         time_t timetoday;
         time (&timetoday);
         string gloT = PrintTime(asctime(gmtime(&timetoday)), tL);
         string locT = PrintTime(asctime(localtime(&timetoday)), tL);

         //* Don't Change Writing patterns *// 
         outfile0 << "\t///////////////////////////////////////////////////" << endl;
         outfile0 << "\t//    ...oooOOO Earth PREM Profile OOOooo...     //" << endl;
         outfile0 << "\t// Source File        : " << file << "\t         //" << endl;
         outfile0 << "\t// Generated on LOCAL : " << locT << " //" << endl;
         outfile0 << "\t//              GMT   : " << gloT << " //" << endl;
         outfile0 << "\t// No. of layers in Earth PREM Profile: " << N << "\t //"<< endl; 
         outfile0 << "\t///////////////////////////////////////////////////\n" << endl;
         
         //* Don't Change Writing patterns *//
         outfile1 << "  ///////////////////////////////////////////////////" << endl;
         outfile1 << "  //    ...oooOOO Earth PREM Profile OOOooo...     //" << endl;
         outfile1 << "  // Source File        : " << file << "\t           //" << endl;
         outfile1 << "  // Generated on LOCAL : " << locT << " //" << endl;
         outfile1 << "  //              GMT   : " << gloT << " //" << endl;
         outfile1 << "  // No. of layers in Earth PREM Profile: " << N << "\t   //"<< endl;
         outfile1 << "  ///////////////////////////////////////////////////\n" << endl;

         //* Don't Change Writing patterns *//
         outfile1 << "  //* Radius in Km & Density in g/cc *//\n" << endl;
         outfile1 << "  #include <OscProbability.h>\n" << endl;
         outfile1 << "    const double LOADED_PREM["<< N <<"][2] = "<< "\n" << "   {" << endl;

         //* Don't Change Writing patterns *//
         for(int i = 0; i < N-1; i++) { outfile1 << "     {" << PREM[i].first << "," << PREM[i].second << "}," << endl; }

         //* Don't Change Writing pattern *//
         outfile1 << "     {" << PREM[N-1].first << "," << PREM[N-1].second << "}" << endl;
         outfile1 << "   };\n\n" << endl;
         outfile1 << "    //* Loading PREM Profile data *//" << endl; 
         outfile1 << "    void NeutrinoOscProbSTD :: SetEarthDensityProfile()\n  {" << endl;
         outfile1 << "    const int rows =  sizeof(LOADED_PREM)/sizeof(LOADED_PREM[0]);" << endl;
         outfile1 << "    PREM.clear();" << endl;
         outfile1 << "    for(int i = 0; i < rows; i++) { PREM.push_back(make_pair(LOADED_PREM[i][0], LOADED_PREM[i][1])); }" << endl;
         outfile1 << "  }\n\n\n\n\n\n" << endl;

         cout << "\033[0;92mOscEarthPREMProfle.cxx & PREM.info Files are Successfully Generated\033[0m\n\033[0;93mDone!!!\033[0m\n" 
              << endl;
      }  else { cout << "\033[0;91mError Warning :: Earth Layer Radii Are Not In Required Order,\nCheck Source File (\033[0m" 
                     << file << "\033[0;91m) !!!\033[0m\n" 
                     << endl;
                exit (1);
              }

    PREM_FILE.close();
    return 0;
  }

  //* Please Don't Change Written Format *// 






