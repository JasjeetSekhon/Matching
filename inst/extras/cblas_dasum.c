/*
 * cblas_dasum.c
 *
 * The program is a C interface to dasum.
 * It calls the fortran wrapper before calling dasum.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#include <R.h>
#include <R_ext/Applic.h> /* R blas declarations */

#include "cblas.h"
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
