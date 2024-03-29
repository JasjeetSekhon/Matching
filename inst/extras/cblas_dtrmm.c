/*
 *
 * cblas_dtrmm.c
 * This program is a C interface to dtrmm.
 * Written by Keita Teranishi
 * 4/6/1998
 *
 */

#include <R.h>
#include <R_ext/Applic.h> /* No longer contains R blas declarations */
#include <R_ext/BLAS.h> /* New R blas declarations */

#include "cblas.h"
void cblas_dtrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const  enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double  *A, const int lda,
                 double  *B, const int ldb)
{
   char UL, TA, SD, DI;   
#ifdef F77_CHAR
   F77_CHAR F77_TA, F77_UL, F77_SD, F77_DI;
#else
   #define F77_TA &TA  
   #define F77_UL &UL  
   #define F77_SD &SD
   #define F77_DI &DI
#endif

#ifdef F77_INT
   F77_INT F77_M=M, F77_N=N, F77_lda=lda, F77_ldb=ldb;
#else
   #define F77_M M
   #define F77_N N
   #define F77_lda lda
   #define F77_ldb ldb
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
	error("cblas_dtrmm: Illegal Side setting, %d\n", Side);
	/* cblas_xerbla(2, "cblas_dtrmm","Illegal Side setting, %d\n", Side); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }
      if( Uplo == CblasUpper) UL='U';
      else if ( Uplo == CblasLower ) UL='L';
      else 
      {
	error("cblas_dtrmm: Illegal Uplo setting, %d\n", Uplo);
	/* cblas_xerbla(3, "cblas_dtrmm","Illegal Uplo setting, %d\n", Uplo); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      if( TransA == CblasTrans) TA ='T';
      else if ( TransA == CblasConjTrans ) TA='C';
      else if ( TransA == CblasNoTrans )   TA='N';
      else 
      {
	error("cblas_dtrmm: Illegal Trans setting, %d\n", TransA);
	/* cblas_xerbla(4, "cblas_dtrmm","Illegal Trans setting, %d\n", TransA); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      if( Diag == CblasUnit ) DI='U';
      else if ( Diag == CblasNonUnit ) DI='N';
      else 
      {
	error("cblas_dtrmm: Illegal Diag setting, %d\n", Diag);
	/* cblas_xerbla(5, "cblas_dtrmm","Illegal Diag setting, %d\n", Diag); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      #ifdef F77_CHAR
         F77_UL = C2F_CHAR(&UL);
         F77_TA = C2F_CHAR(&TA);
         F77_SD = C2F_CHAR(&SD);
         F77_DI = C2F_CHAR(&DI);
      #endif

	 F77_CALL(dtrmm)(F77_SD, F77_UL, F77_TA, F77_DI, &F77_M, &F77_N, &alpha, A, &F77_lda, B, &F77_ldb);
   } else if (Order == CblasRowMajor)
   {
      RowMajorStrg = 1;
      if( Side == CblasRight) SD='L';
      else if ( Side == CblasLeft ) SD='R';
      else 
      {
	error("cblas_dtrmm: Illegal Side setting, %d\n", Side);
	/* cblas_xerbla(2, "cblas_dtrmm","Illegal Side setting, %d\n", Side); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      if( Uplo == CblasUpper) UL='L';
      else if ( Uplo == CblasLower ) UL='U';
      else 
      {
	error("cblas_dtrmm: Illegal Uplo setting, %d\n", Uplo);
	/* cblas_xerbla(3, "cblas_dtrmm","Illegal Uplo setting, %d\n", Uplo); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      if( TransA == CblasTrans) TA ='T';
      else if ( TransA == CblasConjTrans ) TA='C';
      else if ( TransA == CblasNoTrans )   TA='N';
      else 
      {
	error("cblas_dtrmm: Illegal Trans setting, %d\n", TransA);
	/* cblas_xerbla(4, "cblas_dtrmm","Illegal Trans setting, %d\n", TransA); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      if( Diag == CblasUnit ) DI='U';
      else if ( Diag == CblasNonUnit ) DI='N';
      else 
      {
	error("cblas_dtrmm: Illegal Diag setting, %d\n", Diag);
	/* cblas_xerbla(5, "cblas_dtrmm","Illegal Diag setting, %d\n", Diag); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      #ifdef F77_CHAR
         F77_UL = C2F_CHAR(&UL);
         F77_TA = C2F_CHAR(&TA);
         F77_SD = C2F_CHAR(&SD);
         F77_DI = C2F_CHAR(&DI);
      #endif
	 F77_CALL(dtrmm)(F77_SD, F77_UL, F77_TA, F77_DI, &F77_N, &F77_M, &alpha, A, &F77_lda, B, &F77_ldb);
   } 
   else {
     error("cblas_dtrmm: Illegal Order setting, %d\n", Order);
     /* cblas_xerbla(1, "cblas_dtrmm", "Illegal Order setting, %d\n", Order); */
   }
   CBLAS_CallFromC = 0;
   RowMajorStrg = 0;
   return;
}
