/**
 *
 *  @file  OscComplxMatrixExp.h
 *
 *  @brief It calculates the exponential function of complex Matrices.
 *
 *  Reference: <br>
 *  The mxn matrix of complex numbers could be well represented by a 2mx2n matrix of real numbers 
 *  <a href="https://en.wikipedia.org/wiki/Conjugate_transpose"> (See Motivation Section) </a>. <br>
 *  Numerical Recipes in C: The Art of Scientific Computing [ISBN : 978-0521431088]. <br>
 *  Using <a href="https://www.gnu.org/software/gsl/"> GNU GSL </a>, 
 *  I modified the publicly available code at <a href="https://stackoverflow.com/questions/10049160/complex-matrix-exponential-in-c"> stackoverflow.com</a>. <br>
 *
 *  GSL_PREC_DOUBLE : Double-precision, a relative accuracy of approximately 2.0E-16 <br>
 *  GSL_PREC_SINGLE : Single-precision, a relative accuracy of approximately 1.0E-7  <br>
 *  GSL_PREC_APPROX : Approximate values, a relative accuracy of approximately 5.0E-4 <br>
 *  User can set the precision of calculations as given in 
 *  <a href="https://www.gnu.org/software/gsl/doc/html/specfunc.html?highlight=gsl_prec_double#c.gsl_mode_t.GSL_PREC_DOUBLE"> GNU GSL manual</a>. <br>
 *
 *  @author
 *
 *  SADASHIV SAHOO <br>
 *  India-based Neutrino Observatory <br>
 *  Homi Bhabha National Institute, Mumbai, INDIA
 *
 **/

    #ifndef _OscComplxMatrixExp_H
    #define _OscComplxMatrixExp_H

    #include <complex>
    #include <iostream>
    #include <OscConstants.h>
    #include <gsl/gsl_blas.h>
    #include <gsl/gsl_linalg.h>
    #include <gsl/gsl_matrix.h>
    #include <eigen3/Eigen/Core>
    #include <eigen3/Eigen/Dense>
    #include <gsl/gsl_complex_math.h>

    //* Complex Matrix Exponential *//
    //* Essential Internal Function to calculate exp(Matrix) *//
    Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic >
    MatExp(const Eigen::Ref < const Eigen::Matrix < std::complex < double >, Eigen::Dynamic, Eigen::Dynamic > > &A);

    inline void gsl_complex_matrix_exponential (gsl_matrix_complex *M, gsl_matrix_complex *eM, int dim);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // This method is faster than Eigen's Unsupported Exponetial Function at large "PREM Profile" //
    // One may choose Eigen's Unsupported Expoential Function, In order to aviod GSL-Installation //
    ////////////////////////////////////////////////////////////////////////////////////////////////

    #endif





