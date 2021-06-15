
    //////////////////////////////////////////////////////////////////
    //   .........oooOOO 18th-Aug-2020 Last Update OOOooo........   //
    //                        SADASHIV SAHOO                        //
    //               India-based Neutrino Observatory               //
    //        Homi Bhabha National Institute, Mumbai, INDIA         // 
    //          @email: sadashiv.sahoo@tifr.res.in                  //
    //                : sadashiv.sahoo@iopb.res.in                  //
    //////////////////////////////////////////////////////////////////

    #ifndef _OscComplxMatrixExp_H
    #define _OscComplxMatrixExp_H

    #include <complex>
    #include <iostream>
    #include <OscConstants.h>
    #include <gsl/gsl_blas.h>
    #include <gsl/gsl_linalg.h>
    #include <gsl/gsl_matrix.h>
    #include <eigen3/Eigen/Dense>
    #include <gsl/gsl_complex_math.h>

    using namespace std;
    using namespace Eigen;

    //* Complex Matrix Exponential *//
    //* Essential Internal Function to calculate exp(Matrix) *//
    Matrix <complex<double>,dim,dim> MatExp(Matrix<complex<double>,dim,dim> A);

    //* https://en.wikipedia.org/wiki/Conjugate_transpose  *// 
    //* Numerical Recipes in C: The Art of Scientific Computing *//
    //* https://stackoverflow.com/questions/10049160/complex-matrix-exponential-in-c [S. Joseph] *//
    void gsl_complex_matrix_exponential                
                               (gsl_matrix_complex *M, gsl_matrix_complex *eM);


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // This method is faster than Eigen's Unsupported Exponetial Function at large "PREM Profile" //
    // One may choose Eigen's Unsupported Expoential Function, In order to aviod GSL-Installation //
    ////////////////////////////////////////////////////////////////////////////////////////////////

    #endif






