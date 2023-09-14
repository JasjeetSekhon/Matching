/* 
   DASUM - BLAS level one, sums the absolute values of the elements of
   a double precision vector
   
   * cblas_dasum.c
   *
   * The program is a C interface to dasum
   * It calls the fortran wrapper before calling dasum.
   *
   * Written by Keita Teranishi.  2/11/1998
   
   CBLAS provides a C interface to the BLAS routines, which were
   originally written in FORTRAN.  CBLAS wrappers are already provided
   on Windows and OS X, but not on other UNIX-like operating systems
   (such as Linux). For most platforms (particularly AMD chips), I
   recommend Kazushige Goto's High-Performance BLAS Library:

   http://www.cs.utexas.edu/users/flame/goto/
   
   For more information on BLAS (including function definitions) see:
   http://www.netlib.org/blas/
   
   Note that I have only included wrappers for the BLAS functions which
   the Matching package actually uses.

   Jas Sekhon
   <sekhon@berkeley.edu>
   https://www.jsekhon.com
   August 1, 2007
*/


#include <R.h>
#include <R_ext/Applic.h> /* No longer contains R blas declarations */
#include <R_ext/BLAS.h> /* New R blas declarations */

//#include "cblas.h"
double cblas_dasum( const int N, const double *X, const int incX) 
{
   double asum;
#ifdef F77_INT
   F77_INT F77_N=N, F77_incX=incX;
#else 
   #define F77_N N
   #define F77_incX incX
#endif
   asum = F77_CALL(dasum)( &F77_N, X, &F77_incX);
   return asum;
}
