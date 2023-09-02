/**
 *
 *  @file OscComplxMatrixExp.cxx
 *
 *  @brief It calculates the exponential function of complex Matrices.
 *
 *  User can set the Precesion of calculations as given in 
 *  <a href="https://www.gnu.org/software/gsl/doc/html/specfunc.html?highlight=gsl_prec_double#c.gsl_mode_t.GSL_PREC_DOUBLE"> GNU GSL manual </a>
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #include <OscComplxMatrixExp.h>

    //*....ooooOOOO GSL Complex Matrix Exponetial OOOOoooo....*// 
      inline void gsl_complex_matrix_exponential
      (gsl_matrix_complex *M, gsl_matrix_complex *eM, int dim)
    {
      gsl_complex s;
      gsl_matrix *Mat    = gsl_matrix_alloc(2*dim,2*dim);      /// | N | N | => 2N | Real | -Img | <==> [ Here,User Choice is ] => | Real | Img  |
      gsl_matrix *expMat = gsl_matrix_alloc(2*dim,2*dim);      /// | N | N | => 2N | Img  | Real | <==> [ Conjugate Transpose ] => | -Img | Real |

    //*...ooo00 Matrix ---> [Mat.real(), Mat.imag(); -Mat.imag(), Mat.real()] 00ooo...*//
      for(int ij = 0; ij < dim; ij++)
   {
       for(int jk = 0; jk < dim; jk++)
     {
       s = gsl_matrix_complex_get(M, ij, jk);                  /// https://en.wikipedia.org/wiki/Conjugate_transpose at motivation -- (section)
       gsl_matrix_set(Mat, ij, jk, GSL_REAL(s));               /// Here, the Conj.Trsp of Matrix is taken as Origin [User Convenience]
       gsl_matrix_set(Mat, ij, dim+jk, GSL_IMAG(s));           /// Std. form of the whole Matrix in [Unit] Column & [2xUnits] Rows (By Choice)
       gsl_matrix_set(Mat, dim+ij, jk, -GSL_IMAG(s));          /// The Matrix is in [Unit] Row & [2xUnits] Columns; Its upto User.
       gsl_matrix_set(Mat, dim+ij, dim+jk, GSL_REAL(s));       /// Using the GSL Method is faster than 3.3.9-Eigen's Unsupported Exponetial.
     }
   }
      double Re; double Im;
      gsl_linalg_exponential_ss(Mat, expMat, GSL_PREC_SINGLE); /// Precesion GSL_PREC_DOUBLE : 0, GSL_PREC_SINGLE : 1, GSL_PREC_APPROX : 2.
      for(int kl = 0; kl < dim; kl++)
    {
        for(int lm = 0; lm < dim; lm++)
      {
        Re = gsl_matrix_get(expMat, kl, lm);
        Im = gsl_matrix_get(expMat, kl, dim+lm);
        gsl_matrix_complex_set(eM, kl, lm, gsl_complex_rect(Re, Im));
      }
    }
      gsl_matrix_free(Mat);
      gsl_matrix_free(expMat);
    }

    //*....ooooOOOO Matrix Exponetial OOOOoooo....*//
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic >
      MatExp(const Eigen::Ref < const Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > > &A)
   {
      int dim = A.rows();
      gsl_matrix_complex *T  = gsl_matrix_complex_alloc(dim,dim);
      gsl_matrix_complex *eT = gsl_matrix_complex_alloc(dim,dim);
      for(int i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
      {
        gsl_matrix_complex_set(T, i, j, gsl_complex_rect(A(i, j).real(), A(i, j).imag()));
      }
    }
      gsl_complex gslNum;
      gsl_complex_matrix_exponential(T, eT, dim);
      Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > C(dim,dim);
      for(int i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
      {
        gslNum = gsl_matrix_complex_get(eT,i,j);
        C(i,j) = std::complex < double > (GSL_REAL(gslNum), GSL_IMAG(gslNum));
      }
    }
      gsl_matrix_complex_free(T);
      gsl_matrix_complex_free(eT);
      return C;
   }

    //* ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... *//






