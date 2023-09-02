
  //////////////////////////////////////////////////////////////////
  //   .........oooOOO 18th-Aug-2020 Last Update OOOooo........   //
  //                        SADASHIV SAHOO                        //
  //               India-based Neutrino Observatory               //
  //        Homi Bhabha National Institute, Mumbai, INDIA         // 
  //          @email: sadashiv.sahoo@tifr.res.in                  //
  //                : sadashiv.sahoo@iopb.res.in                  //
  //////////////////////////////////////////////////////////////////

  #include <OscProbability.h>
  #include <OscComplxMatrixExp.h>

  //*....ooooOOOO GSL Complex Matrix Exponetial OOOOoooo....*// 
  void gsl_complex_matrix_exponential
           (gsl_matrix_complex *M, gsl_matrix_complex *eM)
    {
       gsl_complex s;
       gsl_matrix *Mat    = gsl_matrix_alloc(2*dim,2*dim);   // | N | N | => 2N | Real | -Img | <==> [ Here,User Choice is ] => | Real | Img  |
       gsl_matrix *expMat = gsl_matrix_alloc(2*dim,2*dim);   // | N | N | => 2N | Img  | Real | <==> [ Conjugate Transpose ] => | -Img | Real |

  //*...ooo00 Matrix ---> [Mat.real(), Mat.imag(); -Mat.imag(), Mat.real()] 00ooo...*//
        for(int ij = 0; ij < dim; ij++)
     {
          for(int jk = 0; jk < dim; jk++)
        {
           s = gsl_matrix_complex_get(M, ij, jk);            // https://en.wikipedia.org/wiki/Conjugate_transpose @Motivation -- (section)
           gsl_matrix_set(Mat, ij, jk, GSL_REAL(s));         // Here, the Conj.Trsp of Matrix is taken as Origin [User Convenience] 
           gsl_matrix_set(Mat, ij, dim+jk, GSL_IMAG(s));     // Std. form of the whole Matrix in [Unit] Column & [2xUnits] Rows (By Choice)
           gsl_matrix_set(Mat, dim+ij, jk, -GSL_IMAG(s));    // The Matrix is in [Unit] Row & [2xUnits] Columns; Its upto User.  
           gsl_matrix_set(Mat, dim+ij, dim+jk, GSL_REAL(s)); // GSL Method is faster than 3.3.9-Eigen's Unsupported Exponetial.
        }
     }

           gsl_linalg_exponential_ss(Mat,expMat,.01);        // Precesion setting with 0.01 @ real Matrix of (2N X 2N) dim;

           double Re;
           double Im;

            for(int kl = 0; kl < dim; kl++)
         {
              for(int lm = 0; lm < dim; lm++)
            {
               Re = gsl_matrix_get(expMat, kl, lm);
               Im = gsl_matrix_get(expMat, kl, dim+lm);
               gsl_matrix_complex_set(eM, kl, lm, gsl_complex_rect(Re,Im));
            }
         }

               gsl_matrix_free(Mat);
               gsl_matrix_free(expMat);
    }

  //*....ooooOOOO Matrix Exponetial OOOOoooo....*//
  Matrix <complex<double>,dim,dim> 
         MatExp(Matrix<complex<double>,dim,dim> A)
    {
         complex <double> UV[dim][dim];  // Inter-Mediate Matrix From Eigen to GSL;
         complex <double>  B[dim][dim];  // Mapping Matrix From GSL to Eigen;

         gsl_matrix_complex *T  = gsl_matrix_complex_alloc(dim,dim); // Complex-Matrix
         gsl_matrix_complex *eT = gsl_matrix_complex_alloc(dim,dim); // Exponential of Complex-Matrix

           for(int i = 0; i < dim; i++)
        {
             for(int j = 0; j < dim; j++)
          {
              UV[i][j] = A(i,j);
              gsl_matrix_complex_set(T, i, j, gsl_complex_rect(UV[i][j].real(), UV[i][j].imag()));
//            cout << "(" << UV[i][j].real() << "," <<UV[i][j].imag() << ")" << "\t";

          }
//            cout << endl;
        }

              gsl_complex_matrix_exponential(T, eT); // Exponential Function of Complex Matrix;
              gsl_complex gslNum;

                 for(int i = 0; i < dim; i++)
              {
                   for(int j = 0; j < dim; j++)
                {
                    gslNum  = gsl_matrix_complex_get(eT,i,j);
                    B[i][j] = complex <double> (GSL_REAL(gslNum), GSL_IMAG(gslNum));
//                  cout << "(" << B[i][j].real() << "," << B[i][j].imag() << ")" << "\t";
                }
//                  cout << endl;
              }

                    gsl_matrix_complex_free(T);
                    gsl_matrix_complex_free(eT);

  //*...oooOOO C should change with dim OOOooo...*//
  Matrix <complex<double>,dim,dim> C;  // Write for 3 X 3 Matrix, Increase/Decrease with Dimensions

         C << B[0][0], B[0][1], B[0][2], // Please Modify according to dimensions
              B[1][0], B[1][1], B[1][2], // Increase or decrease elements accordindly
              B[2][0], B[2][1], B[2][2]; // This process is efficient than auto mapping;

         return C;
  }

  //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//





