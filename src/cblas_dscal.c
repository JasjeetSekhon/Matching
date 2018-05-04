/* 
   DSCAL - BLAS level one, scales a double precision vector
   
   * cblas_dscal.c
   *
   * The program is a C interface to dscal.
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
   http://sekhon.berkeley.edu
   August 1, 2007
*/


#include <R.h>
#include <R_ext/Applic.h> /* R blas declarations */

#include "cblas.h"
void cblas_dscal( const int N, const double alpha, double *X, 
                       const int incX)
{
#ifdef F77_INT
   F77_INT F77_N=N, F77_incX=incX;
#else 
   #define F77_N N
   #define F77_incX incX
#endif
   F77_CALL(dscal)( &F77_N, &alpha, X, &F77_incX);
}
