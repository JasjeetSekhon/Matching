/* 
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


#ifndef CBLAS_H
#define CBLAS_H
#endif
#include <stddef.h>

/*
 * Enumerated and derived types
 */
#define CBLAS_INDEX size_t  /* this may vary between platforms */

enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};

