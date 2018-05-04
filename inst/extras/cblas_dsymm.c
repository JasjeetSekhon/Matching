/*
 *
 * cblas_dsymm.c
 * This program is a C interface to dsymm.
 * Written by Keita Teranishi
 * 4/8/1998
 *
 */

#include <R.h>
#include <R_ext/Applic.h> /* R blas declarations */

#include "cblas.h"
void cblas_dsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const double alpha, const double  *A, const int lda,
                 const double  *B, const int ldb, const double beta,
                 double  *C, const int ldc)
{
   char SD, UL;   
#ifdef F77_CHAR
   F77_CHAR F77_SD, F77_UL;
#else
   #define F77_SD &SD  
   #define F77_UL &UL  
#endif

#ifdef F77_INT
   F77_INT F77_M=M, F77_N=N, F77_lda=lda, F77_ldb=ldb;
   F77_INT F77_ldc=ldc;
#else
   #define F77_M M
   #define F77_N N
   #define F77_lda lda
   #define F77_ldb ldb
   #define F77_ldc ldc
#endif

   int CBLAS_CallFromC;
   int RowMajorStrg;
   RowMajorStrg = 0;
   CBLAS_CallFromC = 1;

   if( Order == CblasColMajor )
   {
      if( Side == CblasRight) SD='R';
      else if ( Side == CblasLeft ) SD='L';
      else 
      {
	error("cblas_dsymm: Illegal Side setting, %d\n", Side);
	/* cblas_xerbla(2, "cblas_dsymm","Illegal Side setting, %d\n", Side);*/
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      if( Uplo == CblasUpper) UL='U';
      else if ( Uplo == CblasLower ) UL='L';
      else 
      {
	error("cblas_dsymm: Illegal Uplo setting, %d\n", Uplo);
	/* cblas_xerbla(3, "cblas_dsymm","Illegal Uplo setting, %d\n", Uplo); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      #ifdef F77_CHAR
         F77_UL = C2F_CHAR(&UL);
         F77_SD = C2F_CHAR(&SD);
      #endif

	 F77_CALL(dsymm)(F77_SD, F77_UL, &F77_M, &F77_N, &alpha, A, &F77_lda,
                      B, &F77_ldb, &beta, C, &F77_ldc);
   } else if (Order == CblasRowMajor)
   {
      RowMajorStrg = 1;
      if( Side == CblasRight) SD='L';
      else if ( Side == CblasLeft ) SD='R';
      else 
      {
	error("cblas_dsymm: Illegal Side setting, %d\n", Side);
	/* cblas_xerbla(2, "cblas_dsymm","Illegal Side setting, %d\n", Side); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      if( Uplo == CblasUpper) UL='L';
      else if ( Uplo == CblasLower ) UL='U';
      else 
      {
	error("cblas_dsymm: Illegal Uplo setting, %d\n", Uplo);
	/* cblas_xerbla(3, "cblas_dsymm","Illegal Uplo setting, %d\n", Uplo); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      #ifdef F77_CHAR
         F77_UL = C2F_CHAR(&UL);
         F77_SD = C2F_CHAR(&SD);
      #endif

	 F77_CALL(dsymm)(F77_SD, F77_UL, &F77_N, &F77_M, &alpha, A, &F77_lda, B,
                 &F77_ldb, &beta, C, &F77_ldc);
   } 
   else 
     {
       error("cblas_dsymm: Illegal Order setting, %d\n", Order);
       /* cblas_xerbla(1, "cblas_dsymm","Illegal Order setting, %d\n", Order);  */
     }
   CBLAS_CallFromC = 0;
   RowMajorStrg = 0;
   return;
} 
