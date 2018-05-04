/*  
    Jasjeet S. Sekhon <sekhon@berkeley.edu>
    HTTP://sekhon.berkeley.edu/
    UC Berkeley
    
    2013/10/28
    Under the GNU Public License Version 3

    A *lot* of work and trail-and-error has gone into these functions
    to ensure that they are reliable and fast when run either serially
    or in parallel.  Parallel execution is especially tricky because
    an algorithm which may be fast in serial mode can cause odd
    bottlenecks when run in parallel (such as a cache-bottleneck when
    executing SSE3 instructions via BLAS).  Also, the loops and other
    structures in these functions have been written so that g++ does a
    good job of optimizing them.  Indeed, for these functions icc is
    no faster than g++.

    Details of individuals function are provided right before they are
    defined, but a general description of them is offered here.

    'EstFuncC': This function is directly called from R (by the
    est.func() function in Match().  And it estimates the causal
    effect and properly weights the observations by the number of ties
    (and the weights on input).

    There are *four* different function to perform the actual
    matching.  There are four versions because of speed
    considerations.  That is, they make use of special cases (such as
    the absence of observational specific weights), to make speed
    gains.  And for clarity of code it was determined to make them
    four separate functions instead of adding a lot of if-then
    statements in one big function.  These functions are long, but it
    was generally found to be faster to do it this way.  Odd things
    occur with the optimizations algorithms in various compilers when
    cross-function optimizations are requested---some optimizations
    are performed and some are not.  It appears to significantly help
    gcc to not breakup the functions (at least not to break them up
    the way I was doing do).

    The core matching functions (note that the slow R equivalent is
    RmatchLoop() which is located in the Matching.R file.

    'FasterMatchC': for GenMatch(), when there are *no* observation
    specific weights, we are matching with replacement, keeping ties
    and are using no caliper, no exact matching and no restriction
    matrix.  Also note that GenMatch() assumes that the Weight.matrix
    is a diagonal matrix.  Note that this is the most commonly used
    matching function for GenMatch, and the fastest.  This function
    cannot be used with Match(), because it does not reorder the
    indexes at the end so to be usable to actually estimate causal
    effects.  This makes the function faster than its equivalent which
    is called from Match(), 'MatchLoopCfast'.

    'FastMatchC': for GenMatch(), when there are observation specific
    weights, but we are matching with replacement, keeping ties and
    are using no caliper, no exact matching and no restriction matrix.
    Also note that GenMatch() assumes that the Weight.matrix is a
    diagonal matrix.  This function cannot be used with Match(),
    because it does not reorder the indexes at the end so to be usable
    to actually estimate causal effects.  This makes the function
    faster than its equivalent which is called from Match(),
    'MatchLoopC'.

    'MatchLoopC': This is the most featured matching function.  It is
    called by Match() when we have observational specific weights, and
    by GenMatch() when there are observational specific weights and
    one or more of no replacement, no ties, exact matching, caliper
    matching or use of the restriction matrix.

    'MatchLoopCfast': This function is called by Match() when there
    are no observation specific weights, and by GenMatch() when there
    are no observational specific weights but one or more of no
    replacement, no ties, exact matching, caliper matching or use of
    the restriction matrix.
*/

/* Note: We are using the include cblas header which will direct to the
   central R BLAS */

#include "scythematrix.h"

using namespace SCYTHE;
using namespace std;

/* #include <stdio.h> */
#include <stdlib.h>
#include <cmath>
#include <vector>
/* Now included from  "scythematrix.h"
#include <R.h>
*/
#include <Rdefines.h>

#include "matching.h"

#ifdef __NBLAS__
/* #if defined(__darwin__) || defined(__APPLE__)
#include <vecLib/cblas.h> 
#else */
#define INTERNAL_CBLAS 
#include <R_ext/Applic.h>
/* #endif  */
#endif

/* 
   __GenMatchBLAS_ needlessly calculates Distance for observations who were assigned to the same treatment
#ifdef __NBLAS__
#if !defined(__darwin__) & !defined(__APPLE__)
#define __GenMatchBLAS__
#endif
#endif
*/

extern "C"
{
#ifdef INTERNAL_CBLAS
#include "cblas.h"

  void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
		   const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
		   const int K, const double alpha, const double  *A,
		   const int lda, const double  *B, const int ldb,
		   const double beta, double  *C, const int ldc);

  void cblas_dscal( const int N, const double alpha, double *X, 
		    const int incX);

  double cblas_dasum( const int N, const double *X, const int incX);
#endif

/*---------------------------------------------------------------------------
   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :   use the median() macro defined below to get the median. 

                Reference:

                http://www.eso.org/~ndevilla/median/

                  Author: Wirth, Niklaus 
                   Title: Algorithms + data structures = programs 
               Publisher: Englewood Cliffs: Prentice-Hall, 1976 
    Physical description: 366 p. 
                  Series: Prentice-Hall Series in Automatic Computation 

 ---------------------------------------------------------------------------*/


  double kth_smallest(double *a, long n, long k)
  {
    long i,j,l,m ;
    double x, tmp;
    
    l=0 ; m=n-1 ;
    while (l<m) {
      x=a[k] ;
      i=l ;
      j=m ;
      do {
	while (a[i]<x) i++ ;
	while (x<a[j]) j-- ;
	if (i<=j) {
	  tmp=a[i];
	  a[i]=a[j];
	  a[j]=tmp;
	      i++ ; j-- ;
	}
      } while (i<=j) ;
      if (j<k) l=i ;
      if (k<i) m=j ;
    }
    return a[k] ;
  } // end of kth_smallest


/*---------------------------------------------------------------------------
   Function : EstFuncC

   IN       : I_N: the number of observations

              I_All: what we are estimating---i.e.,the 'estmand' argument to Match()

              I_length: the number of matched pairs which were found by Match()

              I_Y: the outcome---i.e,. the 'Y' argument to Match()

              I_Tr: the treatment indicator--i.e., the 'Tr' argument to Match()

              I_weight: the observation specific weights---i.e, the 'weight' argument 
                        to Match()

              I_indx: the return object from one of the core Match()
                      functions such as 'RmatchLoop', 'FasterMatchC', 'FastMatchC', 
                      'MatchLoopC', and 'MatchLoopCfast'.


   OUT      : A matrix of Yobs by 3 or 4 (if bias correction is being done).  
              The first column is the causal effect ('YCAUS').
	       If BiasAdjustment==TRUE, the second column will be the
	       BiasAdjustment adjusted yby the weights.  And the next
	       two columns are the weights ('Kcount' and 'KKcount').

   JOB      : To estimate the causal effect and to properly weight the
              observations by the number of ties (and the weights on input).  
	       This function is directly called from R (by the est.func() function 
	       in Match().
 ---------------------------------------------------------------------------*/
  
  SEXP EstFuncC (SEXP I_N, SEXP I_All, SEXP I_length, SEXP I_Y, SEXP I_Tr, SEXP I_weight, SEXP I_indx)
	      
  {
    SEXP ret;
    
    long N, All, i, l, j, k, SumFoo, length;
    double SumFooWeight, SumIndx3;
    long ic=3;

    N = asInteger(I_N);
    All = asInteger(I_All);
    length = asInteger(I_length);

    // In vars
    Matrix Tr = Matrix(N, 1);
    Matrix weight = Matrix(N, 1);
    Matrix indx = Matrix(length, ic);
    Matrix Y = Matrix(N, 1);
    
    // Extra vars
    Matrix Kcount = Matrix(N, 1);
    Matrix KKcount = Matrix(N, 1);
    Matrix YCAUS = Matrix(N, 1);
    Matrix foo = Matrix(N,1);
    Matrix WFW = Matrix(N,1);
    Matrix WWFW = Matrix(N, 1);

    k=0;
    //rows and colums are fliped!! j,i != i,j
    for(j=0;j<ic; j++)
      {
	for(i=0;i<length; i++)
	  {
	    indx[M(i,j,ic)] = REAL(I_indx)[k];
	    
	    // count from offest 0 instead of R's offset 1
	    if(j<2)
	      {
		indx[M(i,j,ic)] = indx[M(i,j,ic)] - 1;	
	      }
	    k++;
	  }
      }
    
    for(i=0;i<N; i++)
      {
	Y[i] = REAL(I_Y)[i];
	Tr[i] = REAL(I_Tr)[i];
	weight[i] = REAL(I_weight)[i];
      }
    
    
    //Master Loop
    for (i=0; i<N; i++)
      {
	if ( ((int) Tr[i]==1 & All!=1) | (All==1))
	  {
	    
	    SumFooWeight = 0;
	    SumFoo = 0;
	    SumIndx3 = 0;
	    for (l=0; l<length; l++)
	      {
		if((int) indx[M(l,0,ic)]==i)
		  {
		    SumFoo++;
		    
		    k = (int) indx[M(l,1,ic)];
		    
		    foo[k] = 1.0;
		    
		    SumIndx3 = SumIndx3 + indx[M(l,2,ic)];
		    SumFooWeight = SumFooWeight + weight[k];
		  } // end of if
	      }//end l for
	  } //Tr[i]
	
	if (SumFoo > 0)
	  {
	    for (l=0; l<length; l++)
	      {
		
		if((int) indx[M(l,0,ic)]==i)
		  {
         	    k = (int) indx[M(l,1,ic)];
		    
		    WFW[k]  = weight[k]*foo[k]/SumFooWeight;
		    WWFW[k] = weight[k]*WFW[k]/SumFooWeight;
		    
		    Kcount[k] = Kcount[k] + weight[i] * WFW[k];
		    
		    KKcount[k] = KKcount[k] + weight[i]*WWFW[k];
		    
		    //if ( ((int) Tr[i]==1 & All!=1) | (All==1))
		      YCAUS[i] = Y[k]*indx[M(l,2,ic)]/SumIndx3 + YCAUS[i];
		  } // endof if
	      } // l for loop
		
		if ((int) Tr[i]==1)
		  {
		    YCAUS[i] = Y[i] - YCAUS[i];
		  } 
		else 
		  {
		    YCAUS[i] = YCAUS[i] - Y[i];
		  }
	      } //if SumFoo
      }// end of master N loop

	    
    PROTECT(ret=allocMatrix(REALSXP, N, 3));
    
    // stack up cbind(YCAUS, Kcount, KKcount)
    k = 0;
    for( i = 0; i < N; i++ )
      {
	REAL(ret)[k] = YCAUS[i];
	k++;
      }
    for( i = 0; i < N; i++ )
      {
	REAL(ret)[k] = Kcount[i];
	k++;
      }
    for( i = 0; i < N; i++ )
      {
	REAL(ret)[k] = KKcount[i];
	k++;
      }
    UNPROTECT(1);
    return(ret); 
  } //end of EstFuncC


/*---------------------------------------------------------------------------
   Function : FasterMatchC

   In       : I_N: the number of observations

              I_xvars: the number of columns in the 'X' matrix---i.e.,
                       the number of variables to match on.

              I_All: what we are estimating---i.e.,the 'estmand' argument to GenMatch()

              I_M: 1-to-M matching---i.e., this is the 'M' option to GenMatch()

	       I_cdd: Distance tolerance

              I_ww: Weight.matrix, GenMatch assumes that it is diagonal

              I_Tr: the treatment indicator--i.e., the 'Tr' argument to GenMatch()

              I_X: The matrix containing the variable we are going to
              match on, the 'X' option to GenMatch().


   OUT: A matrix of number of matches rows and 3 columns. The first
   column contains the treatment observations being matched, the
   second the control observations and the third the weights (adjusted
   for ties).

   JOB: For GenMatch(), when there are *no* observation specific
   weights, we are matching with replacement, keeping ties and are
   using no caliper, no exact matching and no restriction matrix.
   Also note that GenMatch() assumes that the Weight.matrix is a
   diagonal matrix.  Note that this is the most commonly used matching
   function for GenMatch, and the fastest.  This function cannot be
   used with Match(), because it does not reorder the indexes at the
   end so to be usable to actually estimate causal effects.  This
   makes the function faster than its equivalent which is called from
   Match(), 'MatchLoopCfast'.

 ---------------------------------------------------------------------------*/

  SEXP FasterMatchC(SEXP I_N, SEXP I_xvars, SEXP I_All, SEXP I_M, SEXP I_cdd,
		    SEXP I_ww, SEXP I_Tr, SEXP I_X)
  {
    SEXP ret;

    long N, xvars, All, M; 
    double cdd, Distmax, Wi_ratio, dfoo;
    
    long i, j, k;
    
    N = asInteger(I_N);
    xvars = asInteger(I_xvars);
    All = asInteger(I_All);
    M = asInteger(I_M);
    cdd = asReal(I_cdd);
    
    Matrix ww = Matrix(xvars, xvars);
    Matrix Tr = Matrix(N, 1);
    Matrix X = Matrix(N, xvars);
    
    k=0;
    //rows and colums are fliped!! j,i != i,j
    for(j=0;j<xvars; j++)
    {
      for(i=0;i<xvars; i++)
      {
        //ww(i,j) = REAL(I_ww)[k];
        ww[M(i,j,xvars)] = REAL(I_ww)[k];
        k++;
      }
    }
    
    for(i=0;i<N; i++)
    {
      Tr[i] = REAL(I_Tr)[i];
    }
    
    //rows and colums are fliped!! j,i != i,j
    k=0;
    for(j=0;j<xvars; j++)
    {
      for(i=0;i<N; i++)
      {
        //X(i,j) = REAL(I_X)[k];
        X[M(i,j,xvars)] = REAL(I_X)[k];
        k++;
      }
    }
    
    Matrix INN = Matrix::seqa(1, 1, N);
    // TT is just equal to INN

    Matrix ZX(N, xvars), Dist(N, 1);

    /* set misc */
    int TREATi = 0;
    
    Matrix DistPot;
    Matrix ACTMAT, POTMAT;

    // Matrices to avoid rbind
    int NM = N*M*100;
    Matrix I  = Matrix(NM,1);
    Matrix IM = Matrix(NM,1);
    Matrix W  = Matrix(NM,1);
    
    //These are larger than needed; it is just easier this way
    int *order_DistPot = (int *) malloc(N*sizeof(int));  
    //double *S = (double *) malloc(N*sizeof(double));  
    
    int MatchCount=0;
    int MCindx=0;
    int overFirstNM=1;
    for(i=0; i < N; i++)
    {
      // treatment indicator for observation to be matched        
      TREATi = (int) Tr[i];
      
      // proceed with all observations if All==1
      // but only with treated observations if All=0        
      if ( (TREATi==1 & All!=1) | (All==1) )
      {

#ifdef __GenMatchBLAS__
        // this loop is equivalent to the matrix X less the matrix A, which is
        // the product of N rows by 1 column of 1.0 and row R of X.
	// I.E.,         xx = X(i,_); DX = (X - (index_onesN * xx)) 
        double *dest = ZX.data;
        double *src  = X.data;
        double *row  = X.data + (i * xvars);
        for (int jj = 0; jj < N; ++jj) {
          for (int kk = 0; kk < xvars; ++kk, ++dest, ++src) {
            *dest = *src - row[kk];
          }
        }

	//JSS, dgemm is slower and I assume that dtrmm and dsymm are also
	//http://docs.sun.com/source/819-3691/dscal.html
	for (int kk=0; kk < xvars; kk++)
	  {
	    cblas_dscal(N, ww.data[M(kk,kk,xvars)], ZX.data+kk, xvars);
	  }

	ZX.multi_scalar(ZX);

	//http://docs.sun.com/source/819-3691/dasum.html 
	for (int jj=0; jj < N; jj++)
	  {
	    Dist.data[jj] = cblas_dasum(xvars, ZX.data+M(jj, 0, xvars), 1);
	  }

#else
	// Don't calculate Distance for observations with of the same treatment assignment
	// Original No BLAS version from Matching 197, this version from 4.4-52
        for (int jj = 0; jj < N; jj++) {
	  if( abs(TREATi-Tr[jj]) > TOL )
	    {
	      Dist.data[jj] = 0.0;
	      for (int kk=0; kk < xvars; kk++)
		{
		  ZX.data[M(jj,kk,xvars)] = X.data[M(jj,kk,xvars)] - X.data[M(i,kk,xvars)];
		  dfoo  = ZX.data[M(jj,kk,xvars)] * ww.data[M(kk,kk,xvars)];
		  Dist.data[jj] += dfoo*dfoo;
		}
	    } // if TR
	} //end of jj loop
#endif /* end __GenMatchBLAS__ */
        
        // Dist distance to observation to be matched
        // is N by 1 vector	    
        
        // set of potential matches (all observations with other treatment)
        // JSS, note:logical vector
        POTMAT = EqualityTestScalar(Tr, 1-TREATi);
        
	// X's for potential matches
	DistPot = selif(Dist, POTMAT);
	//DistPot_size is a constant!! Fix it
	long DistPot_size = size(DistPot);
	
	//rsort_with_index(DistPot.data, order_DistPot, DistPot_size);
	//R_rsort(DistPot.data, DistPot_size);
	//rPsort(DistPot.data, DistPot_size, M);
	
	Distmax = kth_smallest(DistPot.data, DistPot_size, (M-1));        
        
        // selection of actual matches 
        // logical index
        ACTMAT = LessEqualTestScalar(Dist,  (Distmax+cdd));
        ACTMAT = VectorAnd(POTMAT, ACTMAT);
        
        // distance to actual matches.  This is NEVER used.
        // ACTDIST = selif(Dist, ACTMAT);

	// counts how many times each observation is matched.
	double Mi = sum(ACTMAT);
	Wi_ratio = 1/Mi;

	// collect results
	MatchCount = MatchCount + (int) Mi;
	
	if(MatchCount > NM)
	  {
	    
	    NM = NM+N*M*100;
	    
	    if(overFirstNM > 0)
	      {
		Rprintf("Increasing memory because of ties: allocating a matrix of size 3 times %d doubles.\n", NM);
		Rprintf("I would be faster with the ties=FALSE option.\n");
		warning("Increasing memory because of ties.  I would be faster with the ties=FALSE option.");
	      }
	    
	    int OldMatchCount = MatchCount - (int) Mi;
	    
	    Matrix tI_tmp  = Matrix(NM,1);
	    Matrix tIM_tmp = Matrix(NM,1);
	    Matrix tW_tmp  = Matrix(NM,1);

	    memcpy(tI_tmp.data, I.data, OldMatchCount*sizeof(double));
	    memcpy(tIM_tmp.data, IM.data, OldMatchCount*sizeof(double));
	    memcpy(tW_tmp.data, W.data, OldMatchCount*sizeof(double));
	    /*
	    cblas_dcopy(OldMatchCount, I.data, 0, tI_tmp.data, 0);
	    cblas_dcopy(OldMatchCount, IM.data, 0, tIM_tmp.data, 0);
	    cblas_dcopy(OldMatchCount, W.data, 0, tW_tmp.data, 0);
	    */

	    I = tI_tmp;
	    IM = tIM_tmp;
	    W = tW_tmp;
	  }
	
	//foo1 = ones(ACTMATsum, 1)*(i+1);
	//memcpy(tI.data+MCindx, foo1.data, foo1.size*sizeof(double));
	//memcpy(tW.data+MCindx, Wi.data, Wi.size*sizeof(double));
	
	for (j=0; j < (int) Mi; j++)
	  {
	    I.data[MCindx+j] = i+1;
	    W.data[MCindx+j] = Wi_ratio;
	  }
	
	k=0;
	for (j=0; j<N; j++)
	  {
	    if(ACTMAT.data[j] > (1-TOL))
	      {
		IM.data[MCindx+k] = j+1;
		k++;
	      }
	  }
	MCindx = MCindx+k;
      } // end of (TREATi==1 & All!=1) | (All==1) )
    } //END OF i MASTER LOOP!

    // subset matrices to get back to the right dims
    long orig_rowsize = I.rowsize;
    if (MatchCount > 0)
      {
	if (orig_rowsize != MatchCount)
	  {
	    I.rowsize = MatchCount;
	    IM.rowsize = MatchCount;
	    W.rowsize = MatchCount;
	  }
      }
    else
      {
	I=Matrix(1, 1);
	IM=Matrix(1, 1);
	W=Matrix(1, 1);
      }

    /*ATT is okay already */
    /* ATE*/
    if(All==1)
    {
      long tl  = MatchCount;
      Matrix I2  = Matrix::zeros(tl, 1);
      Matrix IM2 = Matrix::zeros(tl, 1);
      Matrix trt = Matrix::zeros(tl, 1);
      
      for(i=0; i<tl; i++)
      {
        k =(int) I[i] -1 ;
        trt[i] = Tr[k];
      }
      
      for (i=0; i<tl; i++)
      {
        if (trt[i]==1)
	  {
	    I2[i] = I[i];
	    IM2[i] = IM[i];
	  } 
        else
	  {
	    I2[i] = IM[i];
	    IM2[i] = I[i];		
	  }
      }
      
      I = I2;
      IM = IM2;
    } 
    else if(All==2)     /* ATC */
      {
	Matrix Itmp = I;
	Matrix IMtmp = IM;
	
	I = IMtmp;
	IM = Itmp;
      }
    
    /* Free Memory */
    //free(S);
    free(order_DistPot);

    //Do we need to readjust to free memory?.  Don't think so.  Check for memory leaks.
    /*
    if (orig_rowsize != MatchCount)
      {
	I.rowsize = orig_rowsize*1000;
	IM.rowsize = orig_rowsize*1000;
	W.rowsize = orig_rowsize*1000;
      }
    */
    
    PROTECT(ret=allocMatrix(REALSXP, MatchCount, 3));
    /* Loop through the data and display the same in matrix format */
    k = 0;
    for( j = 0; j < MatchCount; j++, k++)
      {
        REAL(ret)[k] = I[j];
      }
    for( j = 0; j < MatchCount; j++, k++)
      {
        REAL(ret)[k] = IM[j];
      }
    for( j = 0; j < MatchCount; j++, k++)
      {
        REAL(ret)[k] = W[j];
      }
    UNPROTECT(1);

    return(ret);    
  } //end of FasterMatchC


/*---------------------------------------------------------------------------
   Function : FastMatchC

   In       : I_N: the number of observations

              I_xvars: the number of columns in the 'X' matrix---i.e.,
                       the number of variables to match on.

              I_All: what we are estimating---i.e.,the 'estmand' argument to GenMatch()

              I_M: 1-to-M matching---i.e., this is the 'M' option to GenMatch()

	       I_cdd: Distance tolerance

              I_ww: Weight.matrix, GenMatch assumes that it is diagonal

              I_Tr: the treatment indicator--i.e., the 'Tr' argument to GenMatch()

              I_X: The matrix containing the variable we are going to
              match on, the 'X' option to GenMatch().

              I_weight: the observation specific weights---i.e, the 'weight' argument 
                        to GenMatch()


   OUT: A matrix of number of matches rows and 3 columns. The first
   column contains the treatment observations being matched, the
   second the control observations and the third the weights (adjusted
   for ties).

   JOB: For GenMatch(), when there are observation specific weights,
   but we are matching with replacement, keeping ties and are using no
   caliper, no exact matching and no restriction matrix.  Also note
   that GenMatch() assumes that the Weight.matrix is a diagonal
   matrix.  This function cannot be used with Match(), because it does
   not reorder the indexes at the end so to be usable to actually
   estimate causal effects.  This makes the function faster than its
   equivalent which is called from Match(), 'MatchLoopC'.

 ---------------------------------------------------------------------------*/


  SEXP FastMatchC(SEXP I_N, SEXP I_xvars, SEXP I_All, SEXP I_M, SEXP I_cdd,
                  SEXP I_ww, SEXP I_Tr, SEXP I_X, SEXP I_weight)
  {
    SEXP ret;
    
    long N, xvars, All, M; 
    double cdd, dfoo;
    
    long i, j, k;
    
    N = asInteger(I_N);
    xvars = asInteger(I_xvars);
    All = asInteger(I_All);
    M = asInteger(I_M);
    cdd = asReal(I_cdd);
    
    Matrix ww = Matrix(xvars, xvars);

    Matrix Tr = Matrix(N, 1);
    Matrix X = Matrix(N, xvars);
    Matrix weight = Matrix(N, 1);
    
    k=0;
    //rows and colums are fliped!! j,i != i,j
    for(j=0;j<xvars; j++)
    {
      for(i=0;i<xvars; i++)
      {
        //ww(i,j) = REAL(I_ww)[k];
        ww[M(i,j,xvars)] = REAL(I_ww)[k];
        k++;
      }
    }
    
    for(i=0;i<N; i++)
    {
      Tr[i] = REAL(I_Tr)[i];
    }
    
    //rows and colums are fliped!! j,i != i,j
    k=0;
    for(j=0;j<xvars; j++)
    {
      for(i=0;i<N; i++)
      {
        //X(i,j) = REAL(I_X)[k];
        X[M(i,j,xvars)] = REAL(I_X)[k];
        k++;
      }
    }
    
    for(i=0;i<N; i++)
    {
      weight[i] = REAL(I_weight)[i];
    }
    
    Matrix INN = Matrix::seqa(1, 1, N);
    // TT is just equal to INN

    Matrix ZX(N, xvars), Dist(N, 1);

    Matrix foo1;
    Matrix foo2;

    /* set misc */
    int TREATi = 0, ACTMATsum = 0;
    
    Matrix DistPot, weightPot, tt,
      weightPot_sort, weightPot_sum, Wi, I, IM, W;
    Matrix ACTMAT, POTMAT;
    
    //These are larger than needed; it is just easier this way
    int *order_DistPot = (int *) malloc(N*sizeof(int));  
    double *S = (double *) malloc(N*sizeof(double));  
    
    int first=1;
    for(i=0; i < N; i++)
    {
      // treatment indicator for observation to be matched        
      TREATi = (int) Tr[i];
      
      // proceed with all observations if All==1
      // but only with treated observations if All=0        
      if ( (TREATi==1 & All!=1) | (All==1) )
      {
#ifdef __GenMatchBLAS__
        // this loop is equivalent to the matrix X less the matrix A, which is
        // the product of N rows by 1 column of 1.0 and row R of X.
	// I.E.,         xx = X(i,_); DX = (X - (index_onesN * xx)) 
        double *dest = ZX.data;
        double *src  = X.data;
        double *row  = X.data + (i * xvars);
        for (int jj = 0; jj < N; ++jj) {
          for (int kk = 0; kk < xvars; ++kk, ++dest, ++src) {
            *dest = *src - row[kk];
          }
        }

	//JSS, dgemm is slower and I assume that dtrmm and dsymm are also
	//http://docs.sun.com/source/819-3691/dscal.html
	for (int kk=0; kk < xvars; kk++)
	  {
	    cblas_dscal(N, ww.data[M(kk,kk,xvars)], ZX.data+kk, xvars);
	  }

	ZX.multi_scalar(ZX);

	//http://docs.sun.com/source/819-3691/dasum.html 
	for (int jj=0; jj < N; jj++)
	  {
	    Dist.data[jj] = cblas_dasum(xvars, ZX.data+M(jj, 0, xvars), 1);
	  }
#else
	// Don't calculate Distance for observations with of the same treatment assignment
	// Original No BLAS version from Matching 197, this version from 4.4-52
        for (int jj = 0; jj < N; jj++) {
	  if( abs(TREATi-Tr[jj]) > TOL )
	    {
	      Dist.data[jj] = 0.0;
	      for (int kk=0; kk < xvars; kk++)
		{
		  ZX.data[M(jj,kk,xvars)] = X.data[M(jj,kk,xvars)] - X.data[M(i,kk,xvars)];
		  dfoo  = ZX.data[M(jj,kk,xvars)] * ww.data[M(kk,kk,xvars)];
		  Dist.data[jj] += dfoo*dfoo;
		}
	    } // if TR
	} //end of jj loop
#endif /* end __GenMatchBLAS__ */
        
        // Dist distance to observation to be matched
        // is N by 1 vector	    
        
        // set of potential matches (all observations with other treatment)
        // JSS, note:logical vector
        POTMAT = EqualityTestScalar(Tr, 1-TREATi);
        
        // X's for potential matches
        DistPot = selif(Dist, POTMAT);
        weightPot = selif(weight, POTMAT);
        
        long weightPot_size = size(weightPot);
        
        for(j=0; j< weightPot_size; j++)
	  {
	    // assume that weightPot_size = size(DistPot)
	    order_DistPot[j] = j;
	    S[j] = (double) DistPot[j];
	  }
        
        rsort_with_index (S, order_DistPot, weightPot_size);
        
        weightPot_sort = Matrix(weightPot_size, 1);
        for(j=0; j < weightPot_size; j++)
	  {
	    weightPot_sort[j] = weightPot[order_DistPot[j]];
	  }
        weightPot_sum = cumsum(weightPot_sort);
        
        tt = Matrix::seqa(1, 1, rows(weightPot_sum));
        
        foo1 = GreaterEqualTestScalar(weightPot_sum, M);
        foo2 = selif(tt, foo1);
        
        long MMM = (long) min(foo2) - 1;
        
        // distance at last match
        double Distmax = S[MMM];
        
        // selection of actual matches 
        // logical index
        ACTMAT = LessEqualTestScalar(Dist,  (Distmax+cdd));
        ACTMAT = VectorAnd(POTMAT, ACTMAT);
        
        // distance to actual matches.  This is NEVER used.
        // ACTDIST = selif(Dist, ACTMAT);
        
        // counts how many times each observation is matched.
        double Mi = sum(multi_scalar(weight, ACTMAT));
        
	//previously in a __NBLAS__ wrapper, but it should always be used
        foo1 = weight;
        foo1.multi_scalar(weight);
        foo1.multi_scalar(ACTMAT);
	
        Wi = selif(weight, ACTMAT);
        Wi = weight[i]*Wi/Mi;
        
        ACTMATsum = (int) sumc(ACTMAT)[0];
        
        // collect results
        if (first==1)
	  {
	    I = Matrix::ones(ACTMATsum, 1)*(i+1);
	    IM = selif(INN, ACTMAT);
	    W = Wi;
	    first = 0;
	  }// end of first==1 
        else 
	  {
	    I = rbind(I, Matrix::ones(ACTMATsum, 1)*(i+1));
	    IM = rbind(IM, selif(INN, ACTMAT));
	    W = rbind(W, Wi);
	  } // end of i else
        
      } // end of (TREATi==1 & All!=1) | (All==1) )
    } //END OF i MASTER LOOP!
    
    /*ATT is okay already */
    /* ATE*/
    if(All==1)
    {
      long tl  = rows(I);
      Matrix I2  = Matrix::zeros(tl, 1);
      Matrix IM2 = Matrix::zeros(tl, 1);
      Matrix trt = Matrix::zeros(tl, 1);
      
      for(i=0; i<tl; i++)
      {
        k =(int) I[i] -1 ;
        trt[i] = Tr[k];
      }
      
      for (i=0; i<tl; i++)
      {
        if (trt[i]==1)
	      {
          I2[i] = I[i];
          IM2[i] = IM[i];
	      } 
        else
	      {
          I2[i] = IM[i];
          IM2[i] = I[i];		
	      }
      }
      
      I = I2;
      IM = IM2;
    } 
    else if(All==2)     /* ATC */
    {
      Matrix Itmp = I;
      Matrix IMtmp = IM;
      
      I = IMtmp;
      IM = Itmp;
    }
    
    /* Free Memory */
    free(S);
    free(order_DistPot);

    i = I.rowsize;
    PROTECT(ret=allocMatrix(REALSXP, i, 3));
    /* Loop through the data and display the same in matrix format */
    k = 0;
    for( j = 0; j < i; j++, k++)
      {
        REAL(ret)[k] = I[j];
      }
    for( j = 0; j < i; j++, k++)
      {
        REAL(ret)[k] = IM[j];
      }
    for( j = 0; j < i; j++, k++)
      {
        REAL(ret)[k] = W[j];
      }
    UNPROTECT(1);

    return(ret);    
  } //end of FastMatchC


/*---------------------------------------------------------------------------
   Function : MatchLoopC

   In       : I_N: the number of observations

              I_xvars: the number of columns in the 'X' matrix---i.e.,
                       the number of variables to match on.

              I_All: what we are estimating---i.e.,the 'estmand' argument to GenMatch()

              I_M: 1-to-M matching---i.e., this is the 'M' option to GenMatch()

	       I_cdd: Distance tolerance

	       I_caliper: indicator for if we are using a caliper

	       I_replace: indicator for if we doing matching with replacement?

	       I_ties: indicator for if are we ignoring ties?

              I_ww: Weight.matrix, GenMatch assumes that it is diagonal

              I_Tr: the treatment indicator--i.e., the 'Tr' argument to GenMatch()

              I_X: The matrix containing the variable we are going to
              match on, the 'X' option to GenMatch().

              I_weight: the observation specific weights---i.e, the 'weight' argument 
                        to GenMatch()

	       I_CaliperVec: A vector which contains the caliper which we must fit in
            
              I_Xorig: The X matrix unscaled (we need this to make caliper comparisons)

              I_restrict_trigger: an indicator for if we are using the restriction matrix

              I_restrict_nrow: the number of restriction we have

              I_restrict: the actual restriction matrix

              I_DaigWeightMatrixFlag: is our I_ww matrix diagonal?


   OUT: A matrix of number of matches rows and 5 columns. The first
   column contains the treatment observations being matched, the
   second the control observations and the third the weights (adjusted
   for ties).

   JOB: This is the most featured matching function.  It is called by
   Match() when we have observational specific weights, and by
   GenMatch() when there are observational specific weights and one or
   more of no replacement, no ties, exact matching, caliper matching
   or use of the restriction matrix.

 ---------------------------------------------------------------------------*/
  
  SEXP MatchLoopC(SEXP I_N, SEXP I_xvars, SEXP I_All, SEXP I_M, SEXP I_cdd,
                  SEXP I_caliper, SEXP I_replace, SEXP I_ties,
                  SEXP I_ww, SEXP I_Tr, SEXP I_X, SEXP I_weight,
                  SEXP I_CaliperVec, SEXP I_Xorig,
                  SEXP I_restrict_trigger, SEXP I_restrict_nrow, SEXP I_restrict,
		  SEXP I_DiagWeightMatrixFlag)
  {
    SEXP ret;
    
    long N, xvars, All, M, caliper, replace, ties, restrict_trigger, restrict_nrow, DiagWeightMatrixFlag,
      sum_caliper_drops=0, replace_count=0;
    double cdd, diff, dfoo;
    
    long i, j, k, r, c;
    
    N = asInteger(I_N);
    xvars = asInteger(I_xvars);
    All = asInteger(I_All);
    M = asInteger(I_M);
    cdd = asReal(I_cdd);
    caliper = (long) asReal(I_caliper);
    replace = asInteger(I_replace);
    ties = asInteger(I_ties);
    restrict_nrow = asInteger(I_restrict_nrow);
    restrict_trigger = asInteger(I_restrict_trigger);
    DiagWeightMatrixFlag = asInteger(I_DiagWeightMatrixFlag);
    
    Matrix ww = Matrix(xvars, xvars);
    Matrix Tr = Matrix(N, 1);
    Matrix X = Matrix(N, xvars);
    Matrix weight = Matrix(N, 1);
    
    k=0;
    //rows and colums are fliped!! j,i != i,j
    for(j=0;j<xvars; j++)
      {
	for(i=0;i<xvars; i++)
	  {
	    //ww(i,j) = REAL(I_ww)[k];
	    ww[M(i,j,xvars)] = REAL(I_ww)[k];
	    k++;
	  }
      }
    
    for(i=0;i<N; i++)
      {
	Tr[i] = REAL(I_Tr)[i];
      }
    
    //rows and colums are fliped!! j,i != i,j
    k=0;
    for(j=0;j<xvars; j++)
      {
	for(i=0;i<N; i++)
	  {
	    //X(i,j) = REAL(I_X)[k];
	    X[M(i,j,xvars)] = REAL(I_X)[k];
	    k++;
	  }
      }
    
    for(i=0;i<N; i++)
      {
	weight[i] = REAL(I_weight)[i];
      }
    
    Matrix IMi;
    Matrix Xorig;
    Matrix CaliperVec;
    if(caliper==1)
    {
      Xorig = Matrix(N, xvars);
      CaliperVec = Matrix(xvars, 1);
      
      for (i=0; i<xvars; i++)
      {
        CaliperVec[i] = REAL(I_CaliperVec)[i];
      }
      
      //rows and colums are fliped!! j,i != i,j
      k=0;
      for(j=0;j<xvars; j++)
      {
        for(i=0;i<N; i++)
	      {
          //X(i,j) = REAL(I_X)[k];
          Xorig[M(i,j,xvars)] = REAL(I_Xorig)[k];
          k++;
	      }
      }	
    } // end of caliper==1
    
    Matrix restrict; 
    if (restrict_trigger==1)
    {
      restrict = Matrix(restrict_nrow, 3);
      
      k=0;
      for(j=0;j<3; j++)
      {
        for(i=0;i<restrict_nrow; i++)
	      {
		restrict[M(i,j,3)] = REAL(I_restrict)[k];
		k++;
	      }
      }	
    } /* if (restrict_trigger==1) */

    // required for replace=0, to keep track of obs already used
    long *ReplaceVector;
    if(replace==0)
      {
	ReplaceVector = (long *) malloc(N*sizeof(long));
      }

    //Start random number generator if we are breaking ties
    if (ties==0)
      GetRNGstate();

    Matrix INN = Matrix::seqa(1, 1, N);
    // TT is just equal to INN

#ifdef __NBLAS__
    Matrix DX, ZX(N, xvars), Dist(N, 1);
    if (DiagWeightMatrixFlag!=1) {
      /* ctor zeros data so this does not add any computational time and maintains the dims outside 
	 the scope of this if statement */
      DX = Matrix::zeros(N, xvars);
    }
#else
    Matrix index_onesN = Matrix::ones(N, 1);
    Matrix xx;
    Matrix DX;
#endif

    Matrix Kcount  = Matrix::zeros(N,1);
    Matrix KKcount = Matrix::zeros(N,1);
    Matrix foo1;
    Matrix foo2;

    /* set misc */
    int TREATi = 0, ACTMATsum = 0;

    Matrix DistPot, weightPot, tt, 
      weightPot_sort, weightPot_sum, Wi, I, IM, IMt, W;
    Matrix ACTMAT, POTMAT(N, 1);

    //These are larger than needed; it is just easier this way
    int *order_DistPot = (int *) malloc(N*sizeof(int));  
    double *S = (double *) malloc(N*sizeof(double));  

    int first=1;
    for(i=0; i < N; i++)
      {
	// treatment indicator for observation to be matched        
	TREATi = (int) Tr[i];
	
	// proceed with all observations if All==1
	// but only with treated observations if All=0        
	if ( (TREATi==1 & All!=1) | (All==1) )
	  {
#ifdef __NBLAS__
	    if (DiagWeightMatrixFlag!=1)
	      {
		// this loop is equivalent to the matrix X less the matrix A, which is
		// the product of N rows by 1 column of 1.0 and row R of X.
		double *dest = ZX.data;
		double *src  = X.data;
		double *row  = X.data + (i * xvars);
		for (int jj = 0; jj < N; ++jj) {
		  for (int kk = 0; kk < xvars; ++kk, ++dest, ++src) {
		    *dest = *src - row[kk];
		  }
		}

		/* note that xvars must be > 1 */
		//JSS
		// do the second multiplication with dgemm, multiplying the matrix
		// above, D, by the transpose of matrix W.
		cblas_dgemm(CblasRowMajor,// column major
			    CblasNoTrans, // A not transposed
			    CblasTrans,   // B transposed
			    xvars,        // M
			    N,            // N
			    xvars,        // K
			    1.0,          // alpha, (alpha * A * B)
			    ww.data,      // A
			    xvars,        // leading dimension for A
			    ZX.data,      // B
			    xvars,        // leading dimension for B
			    0.0,          // beta, (beta * C)
			    DX.data,      // C
			    N);           // leading dimension for C
		
		DX.multi_scalar(DX);
		
		std::swap(DX.colsize, DX.rowsize);
		    
		Dist = sumc(DX);
		
		std::swap(Dist.colsize, Dist.rowsize); // transpose 1 x N -> N x 1
		std::swap(DX.colsize, DX.rowsize);
	      } else 
	      {
#ifdef __GenMatchBLAS__
		// this loop is equivalent to the matrix X less the matrix A, which is
		// the product of N rows by 1 column of 1.0 and row R of X.
		double *dest = ZX.data;
		double *src  = X.data;
		double *row  = X.data + (i * xvars);
		for (int jj = 0; jj < N; ++jj) {
		  for (int kk = 0; kk < xvars; ++kk, ++dest, ++src) {
		    *dest = *src - row[kk];
		  }
		}

		//http://docs.sun.com/source/819-3691/dscal.html
		for (int kk=0; kk < xvars; kk++)
		  {
		    cblas_dscal(N, ww.data[M(kk,kk,xvars)], ZX.data+kk, xvars);
		  }
		
		ZX.multi_scalar(ZX);
		
		//http://docs.sun.com/source/819-3691/dasum.html 
		for (int jj=0; jj < N; jj++)
		  {
		    Dist.data[jj] = cblas_dasum(xvars, ZX.data+M(jj, 0, xvars), 1);
		  }
#else
		// Don't calculate Distance for observations with of the same treatment assignment
		// Original No BLAS version from Matching 197, this version from 4.4-52
		for (int jj = 0; jj < N; jj++) {
		  if( abs(TREATi-Tr[jj]) > TOL )
		    {
		      Dist.data[jj] = 0.0;
		      for (int kk=0; kk < xvars; kk++)
			{
			  ZX.data[M(jj,kk,xvars)] = X.data[M(jj,kk,xvars)] - X.data[M(i,kk,xvars)];
			  dfoo  = ZX.data[M(jj,kk,xvars)] * ww.data[M(kk,kk,xvars)];
			  Dist.data[jj] += dfoo*dfoo;
			}
		    } // if TR
		} //end of jj loop
#endif /* end of __GenMatchBLAS__ ifdef */
	      } /* end of if (DiagWeightMatrixFlag!=1) */
#else
	    // covariate value for observation to be matched                        
	    xx = X(i,_);
	    
	    
	    DX = (X - (index_onesN * xx)) * t(ww);
	    
	    if (xvars>1)
	      {
		//JSS
		foo1 = t(multi_scalar(DX, DX));
		Dist = t(sumc(foo1));
		
	      } 
	    else 
	      {
		Dist = multi_scalar(DX, DX);
	      } // end of xvars
#endif /* end __NBLAS__ */
	    
	    // Dist distance to observation to be matched
	    // is N by 1 vector	    
	    
	    if (restrict_trigger==1)
	      {
		for(j=0; j<restrict_nrow; j++)
		  {
		    if ( ((long) restrict[M(j,0,3)])-1 ==i )
		      {
			
			if (restrict[M(j,2,3)] < 0) {
			  Dist[ ((long) restrict[M(j,1,3)])-1 ] = DOUBLE_XMAX;
			}
			else {
			  Dist[ ((long) restrict[M(j,1,3)])-1 ] = restrict[M(j,2,3)];
			}
		      }
		    else if ( ((long) restrict[M(j,1,3)])-1 ==i ) 
		      {
			
			if (restrict[M(j,2,3)] < 0) {
			  Dist[ ((long) restrict[M(j,0,3)])-1 ] = DOUBLE_XMAX;
			}
			else {
			  Dist[ ((long) restrict[M(j,0,3)])-1 ] = restrict[M(j,2,3)];
			}
		      }
		  }
	      } /* if (restrict_trigger==1) */

	    // Don't match observations we have already matched
	    if(replace_count > 0)
	      {
		for(j=0; j<replace_count; j++)
		  {
		    Dist[ ( ReplaceVector[j] - 1) ] = DOUBLE_XMAX;
		  }
	      } // end of replace_count

            // set of potential matches (all observations with other treatment)
            // JSS, note:logical vector
	    POTMAT = EqualityTestScalar(Tr, 1-TREATi);
	    
	    if (caliper==1)
	      {
		for (j=0; j<N; j++)
		  {
		    if((int) POTMAT[j]==1)
		      {
			for (k=0; k<xvars; k++)
			  {
			    diff = abs(Xorig[M(i, k, xvars)] - Xorig[M(j,k,xvars)]); 
			    if (diff > CaliperVec[k])
			      {
				Dist[j] = DOUBLE_XMAX;
				break;
			      }
			  }
		      }
		  } 
	      }//end of if caliper
	    
            // X's for potential matches
            DistPot = selif(Dist, POTMAT);
            weightPot = selif(weight, POTMAT);

	    long weightPot_size = size(weightPot);

	    for(j=0; j< weightPot_size; j++)
	      {
		// assume that weightPot_size = size(DistPot)
		order_DistPot[j] = j;
		S[j] = (double) DistPot[j];
	      }
	    
	    rsort_with_index (S, order_DistPot, weightPot_size);
	    
	    weightPot_sort = Matrix(weightPot_size, 1);
	    for(j=0; j < weightPot_size; j++)
	      {
		weightPot_sort[j] = weightPot[order_DistPot[j]];
	      }
            weightPot_sum = cumsum(weightPot_sort);
	    
	    tt = Matrix::seqa(1, 1, rows(weightPot_sum));
	    
	    foo1 = GreaterEqualTestScalar(weightPot_sum, M);
	    foo2 = selif(tt, foo1);
	    
	    long MMM = (long) min(foo2) - 1;

	    // distance at last match
            double Distmax = S[MMM];

	    if (restrict_trigger==1 | caliper==1)
	      {
		if ( (Distmax+cdd) > DOUBLE_XMAX_CHECK)
		  {
		    sum_caliper_drops++;
		    continue;
		  }
	      } 

            // selection of actual matches 
            // logical index
	    ACTMAT = LessEqualTestScalar(Dist,  (Distmax+cdd));
	    ACTMAT = VectorAnd(POTMAT, ACTMAT);

	    if (ties==0)
	      {
		int Mii = (int) sum(ACTMAT);
		//Do we have ties?
		if (Mii > M)
		  {
		    IMt = selif(INN, ACTMAT);
		    int nties_broken = 0;
		    int ntiesToBreak = Mii - M;
		    while (nties_broken < ntiesToBreak) {
		      int idrop = (int) ( unif_rand()*(double) Mii);
		      k = (int) IMt[idrop];
		      if (k > (-1+TOL) )
			{
			  ACTMAT[k - 1] = 0;
			  IMt[idrop] = -1;
			  nties_broken++;
			} // end if
		    }// end of while loop
		  }
	      }// end of ties loop	    
	    
	    // distance to actual matches.  This is NEVER used.
	    // ACTDIST = selif(Dist, ACTMAT);
	    
            // counts how many times each observation is matched.
	    double Mi = sum(multi_scalar(weight, ACTMAT));
	    
            Kcount = Kcount + 
	      weight[i] * multi_scalar(weight, ACTMAT)/Mi;
	    
	    //previously in a __NBLAS__ wrapper, but it should always be used
	    foo1 = weight;
	    foo1.multi_scalar(weight);
	    foo1.multi_scalar(ACTMAT);
	    
	    KKcount = KKcount + (weight[i]* foo1)/(Mi*Mi);
	    
            Wi = selif(weight, ACTMAT);
	    Wi = weight[i]*Wi/Mi;
	    
	    ACTMATsum = (int) sumc(ACTMAT)[0];

	    //if no replacement
	    if(replace==0)
	      {
		IMt = selif(INN, ACTMAT);
		for (j=0; j<IMt.rowsize; j++)
		  {
		    ReplaceVector[replace_count] = (long) IMt.data[j];
		    replace_count++;
		  }
	      }//end of replace==0

	    // collect results
	    if (first==1)
	      {
		I = Matrix::ones(ACTMATsum, 1)*(i+1);
		if(replace==0)
		  {
		    IM = IMt;
		  } 
		else 
		  {
		    IM = selif(INN, ACTMAT);
		  }
		W = Wi;
		first = 0;
	      }// end of first==1 
	    else 
	      {
		I = rbind(I, Matrix::ones(ACTMATsum, 1)*(i+1));
		if(replace==0)
		  {
		    IM = rbind(IM, IMt);
		  } 
		else
		  {
		    IM = rbind(IM, selif(INN, ACTMAT));
		  }
		W = rbind(W, Wi);
	      } // end of i else
	    
	  } // end of (TREATi==1 & All!=1) | (All==1) )
      } //END OF i MASTER LOOP!

    //Stop random number generator if we are breaking ties
    if (ties==0)
      PutRNGstate();
    
    Matrix rr = cbind(I, IM);
    rr = cbind(rr, W);
    
    long tl  = rows(I);
    /*ATT is okay already */
    /* ATE*/
    if((All==1) & (I[0]!=0))
      {
	Matrix I2  = Matrix::zeros(tl, 1);
	Matrix IM2 = Matrix::zeros(tl, 1);
	Matrix trt = Matrix::zeros(tl, 1);
	
	for(i=0; i<tl; i++)
	  {
	    k =(int) I[i] -1 ;
	    trt[i] = Tr[k];
	  }
	
	for (i=0; i<tl; i++)
	  {
	    if (trt[i]==1)
	      {
		I2[i] = I[i];
		IM2[i] = IM[i];
	      } 
	    else
	      {
		I2[i] = IM[i];
		IM2[i] = I[i];		
	      }
	  }
	
	I = I2;
	IM = IM2;
      } 
    else if(All==2)     /* ATC */
      {
	Matrix Itmp = I;
	Matrix IMtmp = IM;
	
	I = IMtmp;
	IM = Itmp;
      }
    
    rr = cbind(rr, I);
    rr = cbind(rr, IM);
    
    if (caliper==1)
      {
	Matrix scalar_returns = Matrix::zeros(tl, 1);
	scalar_returns[0] = sum_caliper_drops;
	rr = cbind(rr, scalar_returns);
      }
    
    // rr key
    // 1] I (unadjusted); 2] IM (unadjusted); 3] weight; 4] I (adjusted); 5] IM (adjusted);
    // scalar returns [0]: caliper drops
    
    /* Free Memory */
    free(S);
    free(order_DistPot);
    if(replace==0)
      {
	free(ReplaceVector);
      }
    
    r = rows(rr);
    c = cols(rr);
    
    PROTECT(ret=allocMatrix(REALSXP, r, c));
    /* Loop through the data and display the same in matrix format */
    k = 0;
    for( i = 0; i < c; i++ )
      {
	for( j = 0; j < r; j++ )
	  {
	    // REAL(ret)[k] = rr(j,i);
	    // REAL(ret)[k] = rr[j*c+i];
	    /* Use Macro to Index */
	    REAL(ret)[k] = rr[M(j, i, c)];
	    k++;
	  }
      }
    UNPROTECT(1);
    return(ret);    
  } //end of MatchLoopC


/*---------------------------------------------------------------------------
   Function : MatchLoopCfast

   In       : I_N: the number of observations

              I_xvars: the number of columns in the 'X' matrix---i.e.,
                       the number of variables to match on.

              I_All: what we are estimating---i.e.,the 'estmand' argument to GenMatch()

              I_M: 1-to-M matching---i.e., this is the 'M' option to GenMatch()

	       I_cdd: Distance tolerance

	       I_caliper: indicator for if we are using a caliper

	       I_replace: indicator for if we doing matching with replacement?

	       I_ties: indicator for if are we ignoring ties?

              I_ww: Weight.matrix, GenMatch assumes that it is diagonal

              I_Tr: the treatment indicator--i.e., the 'Tr' argument to GenMatch()

              I_X: The matrix containing the variable we are going to
              match on, the 'X' option to GenMatch().

              I_weight: the observation specific weights---i.e, the 'weight' argument 
                        to GenMatch()

	       I_CaliperVec: A vector which contains the caliper which we must fit in
            
              I_Xorig: The X matrix unscaled (we need this to make caliper comparisons)

              I_restrict_trigger: an indicator for if we are using the restriction matrix

              I_restrict_nrow: the number of restriction we have

              I_restrict: the actual restriction matrix

              I_DaigWeightMatrixFlag: is our I_ww matrix diagonal?


   OUT: A matrix of number of matches rows and 5 columns. The first
   column contains the treatment observations being matched, the
   second the control observations and the third the weights (adjusted
   for ties).

   JOB: This function is called by Match() when there are no
   observation specific weights, and by GenMatch() when there are no
   observational specific weights but one or more of no replacement,
   no ties, exact matching, caliper matching or use of the restriction
   matrix.

 ---------------------------------------------------------------------------*/

  SEXP MatchLoopCfast(SEXP I_N, SEXP I_xvars, SEXP I_All, SEXP I_M, SEXP I_cdd,
		      SEXP I_caliper, SEXP I_replace, SEXP I_ties,
		      SEXP I_ww, SEXP I_Tr, SEXP I_X, 
		      SEXP I_CaliperVec, SEXP I_Xorig,
		      SEXP I_restrict_trigger, SEXP I_restrict_nrow, SEXP I_restrict,
		      SEXP I_DiagWeightMatrixFlag)
  {
    SEXP ret;
    
    long N, xvars, All, M, caliper, replace, ties, restrict_trigger, restrict_nrow, DiagWeightMatrixFlag,
      sum_caliper_drops=0, replace_count=0;
    double cdd, diff, Distmax, Wi_ratio, dfoo;
    
    long i, j, k, r, c;
    
    N = asInteger(I_N);
    xvars = asInteger(I_xvars);
    All = asInteger(I_All);
    M = asInteger(I_M);
    cdd = asReal(I_cdd);
    caliper = (long) asReal(I_caliper);
    replace = asInteger(I_replace);
    ties = asInteger(I_ties);
    restrict_nrow = asInteger(I_restrict_nrow);
    restrict_trigger = asInteger(I_restrict_trigger);
    DiagWeightMatrixFlag = asInteger(I_DiagWeightMatrixFlag);
    
    Matrix ww = Matrix(xvars, xvars);
    Matrix Tr = Matrix(N, 1);
    Matrix X = Matrix(N, xvars);
    
    k=0;
    //rows and colums are fliped!! j,i != i,j
    for(j=0;j<xvars; j++)
      {
	for(i=0;i<xvars; i++)
	  {
	    //ww(i,j) = REAL(I_ww)[k];
	    ww[M(i,j,xvars)] = REAL(I_ww)[k];
	    k++;
	  }
      }
    
    for(i=0;i<N; i++)
      {
	Tr[i] = REAL(I_Tr)[i];
      }
    
    //rows and colums are fliped!! j,i != i,j
    k=0;
    for(j=0;j<xvars; j++)
      {
	for(i=0;i<N; i++)
	  {
	    //X(i,j) = REAL(I_X)[k];
	    X[M(i,j,xvars)] = REAL(I_X)[k];
	    k++;
	  }
      }
    
    Matrix IMi;
    Matrix Xorig;
    Matrix CaliperVec;
    if(caliper==1)
    {
      Xorig = Matrix(N, xvars);
      CaliperVec = Matrix(xvars, 1);
      
      for (i=0; i<xvars; i++)
      {
        CaliperVec[i] = REAL(I_CaliperVec)[i];
      }
      
      //rows and colums are fliped!! j,i != i,j
      k=0;
      for(j=0;j<xvars; j++)
      {
        for(i=0;i<N; i++)
	      {
          //X(i,j) = REAL(I_X)[k];
          Xorig[M(i,j,xvars)] = REAL(I_Xorig)[k];
          k++;
	      }
      }	
    } // end of caliper==1
    
    Matrix restrict; 
    if (restrict_trigger==1)
    {
      restrict = Matrix(restrict_nrow, 3);
      
      k=0;
      for(j=0;j<3; j++)
      {
        for(i=0;i<restrict_nrow; i++)
	      {
		restrict[M(i,j,3)] = REAL(I_restrict)[k];
		k++;
	      }
      }	
    } /* if (restrict_trigger==1) */

    // required for replace=0, to keep track of obs already used
    long *ReplaceVector;
    if(replace==0)
      {
	ReplaceVector = (long *) malloc(N*sizeof(long));
      }

    //Start random number generator if we are breaking ties
    if (ties==0)
      GetRNGstate();

    Matrix INN = Matrix::seqa(1, 1, N);
    // TT is just equal to INN

#ifdef __NBLAS__
    Matrix DX, ZX(N, xvars), Dist(N, 1);
    if (DiagWeightMatrixFlag!=1) {
      /* ctor zeros data so this does not add any computational time and maintains the dims outside 
	 the scope of this if statement */
      DX = Matrix::zeros(N, xvars);
    }
#else
    Matrix index_onesN = Matrix::ones(N, 1);
    Matrix xx;
    Matrix DX, Dist;

    Matrix foo1;
#endif

    /* set misc */
    int TREATi = 0;

    Matrix IMt;
    Matrix ACTMAT, POTMAT(N, 1);

    // Temporary size for these matrices to avoid rbind
    int NM = N*M*100;
    if (ties==0)
      NM = N*M;

    Matrix I  = Matrix(NM,1);
    Matrix IM = Matrix(NM,1);
    Matrix W  = Matrix(NM,1);

    //These are larger than needed; it is just easier this way
    int *order_DistPot = (int *) malloc(N*sizeof(int));  
    //double *S = (double *) malloc(N*sizeof(double));  
    Matrix DistPot=Matrix(N, 1);

    int MatchCount=0;
    int MCindx=0;
    int overFirstNM=1;

    for(i=0; i < N; i++)
      {
	// treatment indicator for observation to be matched        
	TREATi = (int) Tr[i];
	
	// proceed with all observations if All==1
	// but only with treated observations if All=0        
	if ( (TREATi==1 & All!=1) | (All==1) )
	  {
#ifdef __NBLAS__
	    if (DiagWeightMatrixFlag!=1)
	      {
		// this loop is equivalent to the matrix X less the matrix A, which is
		// the product of N rows by 1 column of 1.0 and row R of X.
		double *dest = ZX.data;
		double *src  = X.data;
		double *row  = X.data + (i * xvars);
		for (int jj = 0; jj < N; ++jj) {
		  for (int kk = 0; kk < xvars; ++kk, ++dest, ++src) {
		    *dest = *src - row[kk];
		  }
		}

		/* note that xvars must be > 1 */
		//JSS
		// do the second multiplication with dgemm, multiplying the matrix
		// above, D, by the transpose of matrix W.
		cblas_dgemm(CblasRowMajor,// column major
			    CblasNoTrans, // A not transposed
			    CblasTrans,   // B transposed
			    xvars,        // M
			    N,            // N
			    xvars,        // K
			    1.0,          // alpha, (alpha * A * B)
			    ww.data,      // A
			    xvars,        // leading dimension for A
			    ZX.data,      // B
			    xvars,        // leading dimension for B
			    0.0,          // beta, (beta * C)
			    DX.data,      // C
			    N);           // leading dimension for C
		
		DX.multi_scalar(DX);
		
		std::swap(DX.colsize, DX.rowsize);
		    
		Dist = sumc(DX);
		
		std::swap(Dist.colsize, Dist.rowsize); // transpose 1 x N -> N x 1
		std::swap(DX.colsize, DX.rowsize);
	      } else 
	      {
#ifdef __GenMatchBLAS__
		// this loop is equivalent to the matrix X less the matrix A, which is
		// the product of N rows by 1 column of 1.0 and row R of X.
		double *dest = ZX.data;
		double *src  = X.data;
		double *row  = X.data + (i * xvars);
		for (int jj = 0; jj < N; ++jj) {
		  for (int kk = 0; kk < xvars; ++kk, ++dest, ++src) {
		    *dest = *src - row[kk];
		  }
		}

		//http://docs.sun.com/source/819-3691/dscal.html
		for (int kk=0; kk < xvars; kk++)
		  {
		    cblas_dscal(N, ww.data[M(kk,kk,xvars)], ZX.data+kk, xvars);
		  }
		
		ZX.multi_scalar(ZX);
		
		//http://docs.sun.com/source/819-3691/dasum.html 
		for (int jj=0; jj < N; jj++)
		  {
		    Dist.data[jj] = cblas_dasum(xvars, ZX.data+M(jj, 0, xvars), 1);
		  }
#else
		// Don't calculate Distance for observations with of the same treatment assignment
		// Original No BLAS version from Matching 197, this version from 4.4-52
		for (int jj = 0; jj < N; jj++) {
		  if( abs(TREATi-Tr[jj]) > TOL )
		    {
		      Dist.data[jj] = 0.0;
		      for (int kk=0; kk < xvars; kk++)
			{
			  ZX.data[M(jj,kk,xvars)] = X.data[M(jj,kk,xvars)] - X.data[M(i,kk,xvars)];
			  dfoo  = ZX.data[M(jj,kk,xvars)] * ww.data[M(kk,kk,xvars)];
			  Dist.data[jj] += dfoo*dfoo;
			}
		    } // if TR
		} //end of jj loop
#endif /* end of __GenMatchBLAS__ ifdef */
	      } /* end of if (DiagWeightMatrixFlag!=1) */
#else
	    // covariate value for observation to be matched                        
	    xx = X(i,_);
	    
	    
	    DX = (X - (index_onesN * xx)) * t(ww);
	    
	    if (xvars>1)
	      {
		//JSS
		foo1 = t(multi_scalar(DX, DX));
		Dist = t(sumc(foo1));
		
	      } 
	    else 
	      {
		Dist = multi_scalar(DX, DX);
	      } // end of xvars
#endif /* end __NBLAS__ */
	    
	    // Dist distance to observation to be matched
	    // is N by 1 vector	    

	    if (restrict_trigger==1)
	      {
		for(j=0; j<restrict_nrow; j++)
		  {
		    if ( ((long) restrict[M(j,0,3)])-1 ==i )
		      {
			
			if (restrict[M(j,2,3)] < 0) {
			  Dist[ ((long) restrict[M(j,1,3)])-1 ] = DOUBLE_XMAX;
			}
			else {
			  Dist[ ((long) restrict[M(j,1,3)])-1 ] = restrict[M(j,2,3)];
			}
		      }
		    else if ( ((long) restrict[M(j,1,3)])-1 ==i ) 
		      {
			
			if (restrict[M(j,2,3)] < 0) {
			  Dist[ ((long) restrict[M(j,0,3)])-1 ] = DOUBLE_XMAX;
			}
			else {
			  Dist[ ((long) restrict[M(j,0,3)])-1 ] = restrict[M(j,2,3)];
			}
		      }
		  }
	      } /* if (restrict_trigger==1) */

	    // Don't match observations we have already matched
	    if(replace_count > 0)
	      {
		for(j=0; j<replace_count; j++)
		  {
		    Dist[ ( ReplaceVector[j] - 1) ] = DOUBLE_XMAX;
		  }
	      } // end of replace_count

            // set of potential matches (all observations with other treatment)
            // JSS, note:logical vector
	    POTMAT = EqualityTestScalar(Tr, 1-TREATi);

	    if (caliper==1)
	      {
		for (j=0; j<N; j++)
		  {
		    if((int) POTMAT[j]==1)
		      {
			for (k=0; k<xvars; k++)
			  {
			    diff = abs(Xorig[M(i, k, xvars)] - Xorig[M(j,k,xvars)]); 
			    if (diff > CaliperVec[k])
			      {
				Dist[j] = DOUBLE_XMAX;
				break;
			      }
			  }
		      }
		  } 
	      }//end of if caliper

            // X's for potential matches
            DistPot = selif(Dist, POTMAT);
	    //DistPot_size is a constant!! Fix it
	    long DistPot_size = size(DistPot);

	    //rsort_with_index(DistPot.data, order_DistPot, DistPot_size);
	    //R_rsort(DistPot.data, DistPot_size);
	    //rPsort(DistPot.data, DistPot_size, M);

	    Distmax = kth_smallest(DistPot.data, DistPot_size, (M-1));

	    if (restrict_trigger==1 | caliper==1)
	      {
		if ( (Distmax+cdd) > DOUBLE_XMAX_CHECK)
		  {
		    sum_caliper_drops++;
		    continue;
		  }
	      } 

            // selection of actual matches 
            // logical index
	    ACTMAT = LessEqualTestScalar(Dist,  (Distmax+cdd));
	    ACTMAT = VectorAnd(POTMAT, ACTMAT);

	    if (ties==0)
	      {
		int Mii = (int) sum(ACTMAT);
		//Do we have ties?
		if (Mii > M)
		  {
		    IMt = selif(INN, ACTMAT);
		    int nties_broken = 0;
		    int ntiesToBreak = Mii - M;
		    while (nties_broken < ntiesToBreak) {
		      int idrop = (int) ( unif_rand()*(double) Mii);
		      k = (int) IMt[idrop];
		      if (k > (-1+TOL) )
			{
			  ACTMAT[k - 1] = 0;
			  IMt[idrop] = -1;
			  nties_broken++;
			} // end if
		    }// end of while loop
		  }
	      }// end of ties loop	    
	    
	    // distance to actual matches.  This is NEVER used.
	    // ACTDIST = selif(Dist, ACTMAT);
	    
            // counts how many times each observation is matched.
	    double Mi = sum(ACTMAT);
	    //ACTMATsum = (int) sumc(ACTMAT)[0];

	    Wi_ratio = 1/Mi;
	    //Wi = Matrix::ones(ACTMATsum, 1)*1/Mi;

	    //if no replacement
	    if(replace==0)
	      {
		IMt = selif(INN, ACTMAT);
		for (j=0; j<IMt.rowsize; j++)
		  {
		    ReplaceVector[replace_count] = (long) IMt.data[j];
		    replace_count++;
		  }
	      }//end of replace==0

	    // collect results
	    MatchCount = MatchCount + (int) Mi;

	    if(MatchCount > NM)
	      {
		
		NM = NM+N*M*100;

		if(overFirstNM > 0)
		  {
		    Rprintf("Increasing memory because of ties: allocating a matrix of size 3 times %d doubles.\n", NM);
		    Rprintf("I would be faster with the ties=FALSE option.\n");
		    warning("Increasing memory because of ties.  I would be faster with the ties=FALSE option.");
		  }
		else 
		  {
		    Rprintf("Increasing memory because of ties: allocating a matrix of size 3 times %d doubles.", NM);
		  }
		
		int OldMatchCount = MatchCount - (int) Mi;
		
		Matrix tI_tmp  = Matrix(NM,1);
		Matrix tIM_tmp = Matrix(NM,1);
		Matrix tW_tmp  = Matrix(NM,1);
		
		memcpy(tI_tmp.data, I.data, OldMatchCount*sizeof(double));
		memcpy(tIM_tmp.data, IM.data, OldMatchCount*sizeof(double));
		memcpy(tW_tmp.data, W.data, OldMatchCount*sizeof(double));
		
		I = tI_tmp;
		IM = tIM_tmp;
		W = tW_tmp;
	      }

	    for (j=0; j < (int) Mi; j++)
	      {
		I.data[MCindx+j] = i+1;
		W.data[MCindx+j] = Wi_ratio;
	      }
	    
	    if(replace==0)
	      {
		memcpy(IM.data+MCindx, IMt.data, IMt.size*sizeof(double));
		MCindx = MCindx+IMt.size;
	      }
	    else
	      {
		k=0;
		for (j=0; j<N; j++)
		  {
		    if(ACTMAT.data[j] > (1-TOL))
		      {
			IM.data[MCindx+k] = j+1;
			k++;
		      }
		  }
		MCindx = MCindx+k;
	      }
	  } // end of (TREATi==1 & All!=1) | (All==1) )
      } //END OF i MASTER LOOP!

    //Stop random number generator if we are breaking ties
    if (ties==0)
      PutRNGstate();

    // subset matrices to get back to the right dims
    long orig_rowsize = I.rowsize;
    if (MatchCount > 0)
      {
	if (orig_rowsize != MatchCount)
	  {
	    I.rowsize = MatchCount;
	    IM.rowsize = MatchCount;
	    W.rowsize = MatchCount;
	  }
      }
    else
      {
	I=Matrix(1, 1);
	IM=Matrix(1, 1);
	W=Matrix(1, 1);
      }

    Matrix rr = cbind(I, IM);
    rr = cbind(rr, W);

    long tl  = rows(I);

    /*ATT is okay already */
    /* ATE*/
    if((All==1) & (I[0]!=0))
      {
	Matrix I2  = Matrix::zeros(tl, 1);
	Matrix IM2 = Matrix::zeros(tl, 1);
	Matrix trt = Matrix::zeros(tl, 1);

	for(i=0; i<tl; i++)
	  {
	    k =(int) I[i] -1 ;
	    trt[i] = Tr[k];
	  }
	for (i=0; i<tl; i++)
	  {
	    if (trt[i]==1)
	      {
		I2[i] = I[i];
		IM2[i] = IM[i];
	      } 
	    else
	      {
		I2[i] = IM[i];
		IM2[i] = I[i];		
	      }
	  }
	I = I2;
	IM = IM2;
      } 
    else if(All==2)     /* ATC */
      {
	Matrix Itmp = I;
	Matrix IMtmp = IM;
	
	I = IMtmp;
	IM = Itmp;
      }

    rr = cbind(rr, I);
    rr = cbind(rr, IM);

    if (caliper==1)
      {
	Matrix scalar_returns = Matrix::zeros(tl, 1);
	scalar_returns[0] = sum_caliper_drops;
	rr = cbind(rr, scalar_returns);
      }
    
    // rr key
    // 1] I (unadjusted); 2] IM (unadjusted); 3] weight; 4] I (adjusted); 5] IM (adjusted);
    // scalar returns [0]: caliper drops

    /* Free Memory */
    //free(S);
    free(order_DistPot);
    if(replace==0)
      {
	free(ReplaceVector);
      }
    
    r = rows(rr);
    c = cols(rr);

    PROTECT(ret=allocMatrix(REALSXP, r, c));
    /* Loop through the data and display the same in matrix format */
    k = 0;
    for( i = 0; i < c; i++ )
      {
	for( j = 0; j < r; j++ )
	  {
	    // REAL(ret)[k] = rr(j,i);
	    // REAL(ret)[k] = rr[j*c+i];
	    /* Use Macro to Index */
	    REAL(ret)[k] = rr[M(j, i, c)];
	    k++;
	  }
      }
    UNPROTECT(1);
    return(ret);    
  } //end of MatchLoopCfast


/*---------------------------------------------------------------------------
   Function : VarCalcMatchC

   In       : I_N: the number of observations

              I_xvars: the number of columns in the 'X' matrix---i.e.,
                       the number of variables to match on.

              I_M: 1-to-Var.calc this is the 'Var.calc' option to Match()

	       I_cdd: Distance tolerance

	       I_caliper: indicator for if we are using a caliper

              I_ww: Weight.matrix, GenMatch assumes that it is diagonal

              I_Tr: the treatment indicator--i.e., the 'Tr' argument to GenMatch()

              I_X: The matrix containing the variable we are going to
              match on, the 'X' option to GenMatch().

              I_weight: the observation specific weights---i.e, the 'weight' argument 
                        to GenMatch()

	       I_CaliperVec: A vector which contains the caliper which we must fit in
            
              I_Xorig: The X matrix unscaled (we need this to make caliper comparisons)

              I_restrict_trigger: an indicator for if we are using the restriction matrix

              I_restrict_nrow: the number of restriction we have

              I_restrict: the actual restriction matrix

              I_DaigWeightMatrixFlag: is our I_ww matrix diagonal?


   OUT: Sigs, the key component of the variance estimate for AI SEs when Var.calc > 0

   JOB: This function is modeled after MatchLoopCfast, but it is used
   to calculated AI SEs when Var.calc > 0. This requires matching treated
   to treated, and controls to controls.

 ---------------------------------------------------------------------------*/

  SEXP VarCalcMatchC(SEXP I_N, SEXP I_xvars, SEXP I_M, SEXP I_cdd,
		     SEXP I_caliper, 
		     SEXP I_ww, SEXP I_Tr, SEXP I_X, 
		     SEXP I_CaliperVec, SEXP I_Xorig,
		     SEXP I_restrict_trigger, SEXP I_restrict_nrow, SEXP I_restrict,
		     SEXP I_DiagWeightMatrixFlag, 
		     SEXP I_Y, 
		     SEXP I_weightFlag, SEXP I_weight)
  {
    SEXP ret;
    
    long N, xvars, M, caliper, restrict_trigger, restrict_nrow, DiagWeightMatrixFlag,
      sum_caliper_drops=0;
    double cdd, diff, Distmax, dfoo;
    
    long i, j, k;
    
    N = asInteger(I_N);
    xvars = asInteger(I_xvars);
    M = asInteger(I_M);
    cdd = asReal(I_cdd);
    caliper = (long) asReal(I_caliper);
    restrict_nrow = asInteger(I_restrict_nrow);
    restrict_trigger = asInteger(I_restrict_trigger);
    DiagWeightMatrixFlag = asInteger(I_DiagWeightMatrixFlag);
    
    Matrix ww = Matrix(xvars, xvars);
    Matrix Tr = Matrix(N, 1);
    Matrix X = Matrix(N, xvars);
    Matrix Y = Matrix(N, 1);
    
    k=0;
    //rows and colums are fliped!! j,i != i,j
    for(j=0;j<xvars; j++)
      {
	for(i=0;i<xvars; i++)
	  {
	    //ww(i,j) = REAL(I_ww)[k];
	    ww[M(i,j,xvars)] = REAL(I_ww)[k];
	    k++;
	  }
      }
    
    for(i=0;i<N; i++)
      {
	Tr[i] = REAL(I_Tr)[i];
      }
    
    //rows and colums are fliped!! j,i != i,j
    k=0;
    for(j=0;j<xvars; j++)
      {
	for(i=0;i<N; i++)
	  {
	    //X(i,j) = REAL(I_X)[k];
	    X[M(i,j,xvars)] = REAL(I_X)[k];
	    k++;
	  }
      }
    
    Matrix IMi;
    Matrix Xorig;
    Matrix CaliperVec;
    if(caliper==1)
    {
      Xorig = Matrix(N, xvars);
      CaliperVec = Matrix(xvars, 1);
      
      for (i=0; i<xvars; i++)
      {
        CaliperVec[i] = REAL(I_CaliperVec)[i];
      }
      
      //rows and colums are fliped!! j,i != i,j
      k=0;
      for(j=0;j<xvars; j++)
      {
        for(i=0;i<N; i++)
	      {
          //X(i,j) = REAL(I_X)[k];
          Xorig[M(i,j,xvars)] = REAL(I_Xorig)[k];
          k++;
	      }
      }	
    } // end of caliper==1
    
    Matrix restrict; 
    if (restrict_trigger==1)
    {
      restrict = Matrix(restrict_nrow, 3);
      
      k=0;
      for(j=0;j<3; j++)
      {
        for(i=0;i<restrict_nrow; i++)
	      {
		restrict[M(i,j,3)] = REAL(I_restrict)[k];
		k++;
	      }
      }	
    } /* if (restrict_trigger==1) */

    for(i=0;i<N; i++)
      {
	Y.data[i] = REAL(I_Y)[i];
      }

    Matrix weight, weightPot, tt, weightPot_sort, weightPot_sum, foo1, foo2;
    int weightFlag = asInteger(I_weightFlag);
    if(weightFlag==1)
      {
	weight = Matrix::zeros(N, 1);

	for(i=0;i<N; i++)
	  {
	    weight.data[i] = REAL(I_weight)[i];
	  }	
      }

    Matrix INN = Matrix::seqa(1, 1, N);
    // TT is just equal to INN

#ifdef __NBLAS__
    Matrix DX, ZX(N, xvars), Dist(N, 1);
    if (DiagWeightMatrixFlag!=1) {
      /* ctor zeros data so this does not add any computational time and maintains the dims outside 
	 the scope of this if statement */
      DX = Matrix::zeros(N, xvars);
    }
#else
    Matrix index_onesN = Matrix::ones(N, 1);
    Matrix xx;
    Matrix DX, Dist;
#endif

    /* set misc */
    int TREATi = 0;

    Matrix IMt;
    Matrix ACTMAT, POTMAT(N, 1), Sigs(N,1);

    // Temporary size for these matrices to avoid rbind, in the Var.calc case this IS the max
    int NM = N*M;

    Matrix I  = Matrix(NM,1);
    Matrix IM = Matrix(NM,1);
    Matrix W  = Matrix(NM,1);

    //These are larger than needed; it is just easier this way
    int *order_DistPot = (int *) malloc(N*sizeof(int));  
    double *S = (double *) malloc(N*sizeof(double));  
    Matrix DistPot=Matrix(N, 1);

    for(i=0; i < N; i++)
      {
	// treatment indicator for observation to be matched        
	TREATi = (int) Tr[i];
	
#ifdef __NBLAS__
	if (DiagWeightMatrixFlag!=1)
	  {
	    // this loop is equivalent to the matrix X less the matrix A, which is
	    // the product of N rows by 1 column of 1.0 and row R of X.
	    double *dest = ZX.data;
	    double *src  = X.data;
	    double *row  = X.data + (i * xvars);
	    for (int jj = 0; jj < N; ++jj) {
	      for (int kk = 0; kk < xvars; ++kk, ++dest, ++src) {
		*dest = *src - row[kk];
	      }
	    }
	    
	    /* note that xvars must be > 1 */
	    //JSS
	    // do the second multiplication with dgemm, multiplying the matrix
	    // above, D, by the transpose of matrix W.
	    cblas_dgemm(CblasRowMajor,// column major
			CblasNoTrans, // A not transposed
			CblasTrans,   // B transposed
			xvars,        // M
			N,            // N
			xvars,        // K
			1.0,          // alpha, (alpha * A * B)
			ww.data,      // A
			xvars,        // leading dimension for A
			ZX.data,      // B
			xvars,        // leading dimension for B
			0.0,          // beta, (beta * C)
			DX.data,      // C
			N);           // leading dimension for C
	    
	    DX.multi_scalar(DX);
	    
	    std::swap(DX.colsize, DX.rowsize);
	    
	    Dist = sumc(DX);
	    
	    std::swap(Dist.colsize, Dist.rowsize); // transpose 1 x N -> N x 1
	    std::swap(DX.colsize, DX.rowsize);
	  } else 
	  {
#ifdef __GenMatchBLAS__
	    // this loop is equivalent to the matrix X less the matrix A, which is
	    // the product of N rows by 1 column of 1.0 and row R of X.
	    double *dest = ZX.data;
	    double *src  = X.data;
	    double *row  = X.data + (i * xvars);
	    for (int jj = 0; jj < N; ++jj) {
	      for (int kk = 0; kk < xvars; ++kk, ++dest, ++src) {
		*dest = *src - row[kk];
	      }
	    }
	    
	    //http://docs.sun.com/source/819-3691/dscal.html
	    for (int kk=0; kk < xvars; kk++)
	      {
		cblas_dscal(N, ww.data[M(kk,kk,xvars)], ZX.data+kk, xvars);
	      }
	    
	    ZX.multi_scalar(ZX);
	    
	    //http://docs.sun.com/source/819-3691/dasum.html 
	    for (int jj=0; jj < N; jj++)
	      {
		Dist.data[jj] = cblas_dasum(xvars, ZX.data+M(jj, 0, xvars), 1);
	      }
#else
	    // Don't calculate Distance for observations with of a 
	    // DIFFERENT treatment assignment
	    for (int jj = 0; jj < N; jj++) {
	      if( abs(TREATi-Tr[jj]) < TOL )
		{
		  Dist.data[jj] = 0.0;
		  for (int kk=0; kk < xvars; kk++)
		    {
		      ZX.data[M(jj,kk,xvars)] = X.data[M(jj,kk,xvars)] - X.data[M(i,kk,xvars)];
		      dfoo  = ZX.data[M(jj,kk,xvars)] * ww.data[M(kk,kk,xvars)];
		      Dist.data[jj] += dfoo*dfoo;
		    }
		} // if TR
	    } //end of jj loop
#endif /* end of __GenMatchBLAS__ ifdef */
	  } /* end of if (DiagWeightMatrixFlag!=1) */
#else
	// covariate value for observation to be matched                        
	xx = X(i,_);
	
	
	DX = (X - (index_onesN * xx)) * t(ww);
	
	if (xvars>1)
	  {
	    //JSS
	    foo1 = t(multi_scalar(DX, DX));
	    Dist = t(sumc(foo1));
	    
	  } 
	else 
	  {
	    Dist = multi_scalar(DX, DX);
	  } // end of xvars
#endif /* end __NBLAS__ */

	//Remove self as a potential match
	Dist[i] = DOUBLE_XMAX;
	
	// Dist distance to observation to be matched
	// is N by 1 vector	    
	if (restrict_trigger==1)
	  {
	    for(j=0; j<restrict_nrow; j++)
	      {
		if ( ((long) restrict[M(j,0,3)])-1 ==i )
		  {
		    
		    if (restrict[M(j,2,3)] < 0) {
		      Dist[ ((long) restrict[M(j,1,3)])-1 ] = DOUBLE_XMAX;
		    }
		    else {
		      Dist[ ((long) restrict[M(j,1,3)])-1 ] = restrict[M(j,2,3)];
		    }
		  }
		else if ( ((long) restrict[M(j,1,3)])-1 ==i ) 
		  {
		    
		    if (restrict[M(j,2,3)] < 0) {
		      Dist[ ((long) restrict[M(j,0,3)])-1 ] = DOUBLE_XMAX;
		    }
		    else {
		      Dist[ ((long) restrict[M(j,0,3)])-1 ] = restrict[M(j,2,3)];
		    }
		  }
	      }
	  } /* if (restrict_trigger==1) */

        // set of potential matches (all observations with the *SAME* treatment)
	POTMAT = EqualityTestScalar(Tr, TREATi);
	
	if (caliper==1)
	  {
	    for (j=0; j<N; j++)
	      {
		if((int) POTMAT[j]==1)
		  {
		    for (k=0; k<xvars; k++)
		      {
			diff = abs(Xorig[M(i, k, xvars)] - Xorig[M(j,k,xvars)]); 
			if (diff > CaliperVec[k])
			  {
			    Dist[j] = DOUBLE_XMAX;
			    break;
			  }
		      }
		  }
	      } 
	  }//end of if caliper

	if(weightFlag==0)
	  {
	    // X's for potential matches
	    DistPot = selif(Dist, POTMAT);
	    long DistPot_size = size(DistPot);
	    
	    Distmax = kth_smallest(DistPot.data, DistPot_size, (M-1));
	    
	    if (restrict_trigger==1 | caliper==1)
	      {
		if ( (Distmax+cdd) > DOUBLE_XMAX_CHECK)
		  {
		    sum_caliper_drops++;
		    continue;
		  }
	      } 
	    
	    // selection of actual matches 
	    // logical index
	    ACTMAT = LessEqualTestScalar(Dist,  (Distmax+cdd));
	    ACTMAT = VectorAnd(POTMAT, ACTMAT);
	    
	    // counts how many times each observation is matched.
	    double Mi = sum(ACTMAT);

	    /**********************************************/
	    /* ESTIMATE Sigs */
	    double fm=0, sm=0;
	    double sumweightactmat = Mi + 1.0;
	    
	    //Add self back in
	    ACTMAT.data[i] = 1.0;
	    for(j=0; j<N; j++)
	      {
		if ( abs(ACTMAT.data[j] - 1.0) < TOL)
		  {
		    fm += Y.data[j];
		    sm += (Y.data[j]*Y.data[j]);
		  }
	      }
	    fm = fm/sumweightactmat;
	    sm = sm/sumweightactmat;
	    Sigs[i] = (sm - fm * fm)*sumweightactmat/(sumweightactmat-1.0);
	  } else {
	  // X's for potential matches
	  DistPot = selif(Dist, POTMAT);
	  weightPot = selif(weight, POTMAT);

	  long weightPot_size = size(weightPot);
	  
	  for(j=0; j< weightPot_size; j++)
	    {
	      // assume that weightPot_size = size(DistPot)
	      order_DistPot[j] = j;
	      S[j] = (double) DistPot[j];
	    }

	  rsort_with_index (S, order_DistPot, weightPot_size);
	  
	  weightPot_sort = Matrix(weightPot_size, 1);
	  for(j=0; j < weightPot_size; j++)
	    {
	      weightPot_sort[j] = weightPot[order_DistPot[j]];
	    }
	  weightPot_sum = cumsum(weightPot_sort);

	  tt = Matrix::seqa(1, 1, rows(weightPot_sum));
	  
	  foo1 = GreaterEqualTestScalar(weightPot_sum, M);
	  foo2 = selif(tt, foo1);
	  
	  long MMM = (long) min(foo2) - 1;
	  
	  // distance at last match
	  double Distmax = S[MMM];
	  
	  if (restrict_trigger==1 | caliper==1)
	    {
	      if ( (Distmax+cdd) > DOUBLE_XMAX_CHECK)
		{
		  sum_caliper_drops++;
		  continue;
		}
	    } 
	  
	  // selection of actual matches 
	  // logical index
	  ACTMAT = LessEqualTestScalar(Dist,  (Distmax+cdd));
	  ACTMAT = VectorAnd(POTMAT, ACTMAT);
	  
	  /**********************************************/
	  /* ESTIMATE Sigs */
	  double fm=0, sm=0, ws=0;
	  
	  //Add self back in
	  ACTMAT.data[i] = 1.0;
	  for(j=0; j<N; j++)
	    {
	      if ( abs(ACTMAT.data[j] - 1.0) < TOL)
		{
		  fm += Y.data[j]*weight[j];
		  sm += (Y.data[j]*Y.data[j]*weight[j]);
		  ws += weight[j];
		}
	    }
	  fm = fm/ws;
	  sm = sm/ws;
	  Sigs[i] = (sm - fm * fm)*ws/(ws-1.0);
	} //end of weightFlag
      } //END OF i MASTER LOOP!
    
    /* Free Memory */
    free(S);
    free(order_DistPot);
    
    PROTECT(ret=allocMatrix(REALSXP, N, 1));
    /* Loop through the data and display the same in matrix format */
    for( j = 0; j < N; j++ )
      {
	REAL(ret)[j] = Sigs[j];
      }
    UNPROTECT(1);
    return(ret);    
  } //end of VarCalcMatchC


} //end of extern "C"

double min_scalar (double a, double b)
{
  if (a < b)
    return(a);
  
  return(b);
} // end of min_scalar

double max_scalar (double a, double b)
{
  if (a > b)
    return(a);

  return(b);
} // end of max_scalar

//previously in a __NBLAS__ wrapper, but it should always be used
Matrix EqualityTestScalar(Matrix a, double s)
{
  for (long i = 0; i < a.size; ++i)
    a.data[i] = (a.data[i] < (s+TOL)) && (a.data[i] > (s-TOL)) ? 1 : 0;
  return a;
} //end of EqualityTestScalar


//previously in a __NBLAS__ wrapper, but it should always be used
Matrix GreaterEqualTestScalar(Matrix a, long s)
{
  for (long i = 0; i < a.size; ++i)
    a.data[i] = (a.data[i] >= (s-TOL)) ? 1 : 0;
  return a;
} //end of GreaterEqualTestScalar


//previously in a __NBLAS__ wrapper, but it should always be used
Matrix LessEqualTestScalar(Matrix a, double s)
{
  for (long i = 0; i < a.size; ++i)
    a.data[i] = (a.data[i] <= (s+TOL)) ? 1 : 0;
  return a;
} //end of LessEqualTestScalar


Matrix VectorAnd(Matrix a, Matrix b)
{
  long nrows = rows(a);
  Matrix ret =  Matrix::zeros(nrows, 1);

  for (long i =0; i< nrows; i++)
    {
      if( (a[i] == 1) &&  (b[i]== 1) )
	{
	  ret[i] = 1;
	}
    }
  return(ret);
} //end of VectorAnd

Matrix EqualityTestMatrix(Matrix a, Matrix s)
{
  long nrows = rows(a);
  long ncols = cols(a);
  Matrix ret =  Matrix::zeros(nrows, ncols);

  long scols = cols(s);

  if (scols==1)
    {
      for (long i =0; i< nrows; i++)
	{
	  for (long j =0; j< ncols; j++)
	    {
	      if( (a[M(i, j, ncols)] < (s[i]+TOL)) &  (a[M(i, j, ncols)] > (s[i]-TOL)) )
		{
		  ret[M(i, j, ncols)] = 1;
		}
	    }
	}
    }
  else if (scols==ncols)
    {
      for (long i =0; i< nrows; i++)
	{
	  for (long j =0; j< ncols; j++)
	    {
	      if( (a[M(i, j, ncols)] < (s[M(i, j, ncols)]+TOL)) &  
		  (a[M(i, j, ncols)] > (s[M(i, j, ncols)]-TOL)) )
		{
		  ret[M(i, j, ncols)] = 1;
		}
	    }
	}
    }
  else 
    {
      Rprintf("ASSERTION in EqualityTestMatrix\n");
    }

  return(ret);
} //end of EqualityTestMatrix

/* cumsum */
Matrix cumsum(Matrix a)
{
  long nrows = rows(a);
  Matrix ret = Matrix::zeros(nrows, 1);
  
  ret[0] = a[0];
  for (long i = 1;  i < nrows; i++) 
  {
    ret[i] = ret[i-1] + a[i];
  }
  
  return ret;
} //end of cumsum

//! Calculate the sum of all of the elements in a  Matrix
double sum (const Matrix & A)
{
  double ret=0;
  long ncols = cols(A);
  long i;
  
  Matrix sumvec = sumc(A);
  
  for (i=0; i<ncols; i++)
  {
    ret = ret + sumvec[i];
  }
  
  return ret;
}

/*DISPLAY This function will display the double matrix stored in an mxArray.
* This function assumes that the mxArray passed as input contains double
* array.
*/
void display(Matrix A)
{
  int i=0, j=0; /* loop index variables */
  int r=0, c=0; /* variables to store the row and column length of the matrix */
  int count=0;
  
  /* Get the size of the matrix */
  r = rows(A);
  c = cols(A);
  
  /* Loop through the data and display the same in matrix format */
  for( i = 0; i < r; i++ ){
    for( j = 0; j < c; j++){
      Rprintf("%e\t",A[count]);
      count++;
    }
    Rprintf("\n");
  }
  Rprintf("\n");
}

//previously in a __NBLAS__ wrapper, but it should always be used
Matrix multi_scalar (Matrix a, Matrix b)
{
  a.multi_scalar(b);
  return a;
} 
