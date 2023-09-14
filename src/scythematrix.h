/* Edited by Jasjeet S. Sekhon <jasjeet_sekhon@berkeley.edu> */
/* https://www.jsekhon.com                             */
/*                                                           */
/* April 26, 2013                                            */
/* get rid of friend-injection for ones, zero, seqa          */
/* January 9, 2012                                           */
/* May 8, 2010, updated header files for Solaris             */
/* remove xpnd to work on Solaris march 31, 2009             */

/* June 26, 2008                                             */
/* __NATE__ additions by Nate Begeman (Apple)                */


//
// This library is the class definition of the Matrix class, part of
// the SCYTHE project.
//
// Scythe C++ Library
// Copyright (C) 2000 Kevin M. Quinn and Andrew D. Martin
//
// This code written by:
//
// Kevin Quinn
// Assistant Professor
// Dept. of Political Science and 
// Center for Statistics and the Social Sciences
// Box 354322
// University of Washington
// Seattle, WA  98195-4322
// quinn@stat.washington.edu
//
// Andrew D. Martin
// Assistant Professor
// Dept. of Political Science
// Campus Box 1063
// Washington University
// St. Louis, MO 63130
// admartin@artsci.wustl.edu
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#ifndef SCYTHE_DOUBLE_MATRIX_H
#define SCYTHE_DOUBLE_MATRIX_H

#include <R.h> /* needed for the error() function */
#include <iostream>
#include <iomanip>
#include <fstream>
 /* Not needed and causes problem in pre 3.3 gcc #include <sstream> */
 /* #include <new> we've moved to malloc */
#include <numeric>
#include <string>
#include <climits>
#include <cmath>

 //explicit include needed for gcc4.3 because of header cleanup
#include <cstdlib>
 //needed for Solaris instead of simply #include <cstdlib>
#include <stdlib.h>

 //http://gcc.gnu.org/gcc-4.3/porting_to.html
 //http://www.cyrius.com/journal/2007/05/10#gcc-4.3-include
#include <cstring>

// Avoid NameSpace Pollution
namespace SCYTHE {
  
  struct all_elements{
  }  const _ = {};
  

class Matrix
{
public:
  int rowsize;    // # of rows in Matrix
  int colsize;    // # of columns in Matrix
  int size;    	  // # of element in Matrix
    		  // in row major order

  double *data;   // array holding the elements of the Matrix
  // BASIC FUNCTIONS

/**********************************************************************/
//  CONSTRUCTOR:  Matrix - Creates empty Matrix object
//  Input:  void
//  Return:  empty Matrix object
//  Usage:  Matrix();
//  Errors:  none
//  Dependencies:  none

    inline Matrix (void) {
			rowsize = 1;
			colsize = 1;
			size = 1;
			//data = new double[1];
			data = (double *) malloc(1*sizeof(double));
			data[0] = 0.0;
		}

/**********************************************************************/
//  CONSTRUCTOR:  Matrix  - Creates Matrix obj. with specificed rows
//  and columns
//  Input:  int rows, int columns
//  Return:  empty Matrix object 
//  Usage:  Matrix(int,int);
//  Errors:  0001
//  Dependencies:  none

   Matrix (const int& rows, const int& cols);

/**********************************************************************/
//  CONSTRUCTOR:  Matrix - Creates Matrix obj. filled with data from
//  array.  ( inputarray is in row major order)
//  Input:  double * inputarray, int rows, int columns
//  Return: Matrix object
//  Usage:  Matrix(double *,int,int);
//  Errors:  0002
//  Dependencies:  none

   Matrix (const double *inputarray, const int& rows, const int& cols);

/**********************************************************************/
//  CONSTRUCTOR:   Matrix - Creates Matrix object from old Matrix
//  Input:  Matrix object
//  Return: Matrix object (filled with same data as old_matrix)
//  Usage:  Matrix(Matrix &);
//  Errors:  none
//  Dependencies:  none

   Matrix (const Matrix & old_matrix);

/**********************************************************************/
//  OPERATOR:  Matrix operator '=' - Allows for the copying of a Matrix
//  Input:  Matrix object referemce
//  Return: Matrix object (filled with same data)
//  Usage:  MatrixA = MatrixB
//  Errors:  none
//  Dependencies:  none

    Matrix & operator = (const Matrix & B);

/**********************************************************************/
//  OPERATOR: Matrix operator [] - allows for retrieval of element using
//  traditional mathematical notation (retrieves ith element)
//  Input:  integer i
//  Return: double (value from Matrix)
//  Usage:  Matrix[i]
//  Errors:  0003
//  Dependencies:  none

  inline double &operator[] (const int& i) {
		if (i >= size || i < 0) {
		  error("Index out of range in [] operator");
		}
		return (data[i]);
  }


/**********************************************************************/
//  OPERATOR: Matrix operator () - Retrieves Matrix element using
//    mathematical notation (retrieves (i,j) element).
//    NOTE: Indexing starts at 0
//  Input:  integer i, integer j
//  Return: double (value from Matrix)
//  Usage:  Matrix(integer, integer);
//  Errors:  0004
//  Dependencies:  none

  inline double &operator () (const int& i, const int& j) {
		if (rowsize < i || colsize < j || i < 0 || j < 0) {
		  error("Index out of range in () operator");
		}
		return data[i * colsize + j];
	}




/**********************************************************************/
//  OPERATOR: Matrix operator () - Retrieves Matrix elements using
//    mathematical notation (retrieves (i,_) element).
//    NOTE: Indexing starts at 0
//  Input:  integer i, all_elements
//  Return: Matrix 
//  Usage:  Matrix(integer, _);
//  Errors:  ????
//  Dependencies:  struct all_elements
  
  Matrix operator () (const int& i, const all_elements& a);



/**********************************************************************/
//  OPERATOR: Matrix operator () - Retrieves Matrix elements using
//    mathematical notation (retrieves (_,j) element).
//    NOTE: Indexing starts at 0
//  Input:  all_elements, integer j
//  Return: Matrix 
//  Usage:  Matrix(_,integer);
//  Errors:  ????
//  Dependencies:  struct all_elements
  
  Matrix operator () (const all_elements& a, const int& j);



/**********************************************************************/
//  OPERATOR: Matrix operator () - Retrieves Matrix elements in row i.
//      Indexing starts at 0
//  Input:  integer i, Matrix J
//  Return: Matrix object (row i elements)
//  Usage:  Matrix(integer, Matrix &);
//  Errors:  0005 - 0007
//  Dependencies:  none

   Matrix operator () (const int& i, const Matrix& J);

/**********************************************************************/
//  OPERATOR: Matrix operator () - Retrieves Matrix elements in column j.
//      Indexing starts at 0
//  Input:  Matrix I, integer j
//  Return: Matrix object (column j elements)
//  Usage:  Matrix(Matrix &, integer);
//  Errors:  0008 - 0010
//  Dependencies:  none

   Matrix operator () (const Matrix& I, const int& j);

/**********************************************************************/
//  OPERATOR: Matrix operator () - Extracts submatrix from existing
//      Matrix.  Get elements i,j from a Matrix where i in I and j in J
//  Indexing starts at 0
//  Input:  Matrix I , Matrix J
//  Return: Matrix object ()
//  Usage:  Matrix(Matrix &, Matrix &);
//  Errors:  0011-0017
//  Dependencies:  none

   Matrix operator () (const Matrix& I, const Matrix& J);


/**********************************************************************/
//  DESTRUCTOR: Matrix object destructor
//  Input:  void
//  Return: void
//  Usage:  ~Matrix();
//  Errors:  
//  Dependencies:  none

   inline ~Matrix (void) {
     //delete[] data;
     free(data);
     data = NULL;
   }

/**********************************************************************/
//  FUNCTION: Print Matrix to screen
//  Input:  integer width, integer setprecision
//  Return: void
//  Usage:  print(int width, int precision);
//  Errors:  
//  Dependencies:  none

   void print (const int width = 10, const int prec = 5);


	 inline void getArray(double *out) const {
		 for (int i = 0; i < size; i++) {
			 out[i] = data[i];
		 }
	 }
   
 /**********************************************************************/
 //  FUNCTION: multiply each scalar element of the matrix
 //  Input:  Matrix I
 //  Return: void
 //  Usage:  multi_scalar(Matrix &);
 //  Errors:  
 //  Dependencies:  none
 
   void multi_scalar (Matrix &I) {
     for (int i = 0; i < rowsize * colsize; ++i)
       data[i] *= I.data[i];
   }
   
/*******************   MORE ADVANCED FUNCTIONS  **********************/
/**********************************************************************/

/**********************************************************************/
//  FUNCTION: c -- concatenates a sequence of doubles into a Matrix
//  Input:  a user determined number of doubles (up to 26)
//  Return:  Matrix object (a column vector)
//  Usage:  c(a,b,c,d)
//  Errors:  
//  Dependencies:

   friend Matrix c (const double& a, const double& b);
   friend Matrix c (const double& a, const double& b, const double& c);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d);  
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e);  
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f);  
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g);  
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h);  
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m, const double& n);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m, const double& n, const double& o);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m, const double& n, const double& o, 
		    const double& p);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m, const double& n, const double& o, 
		    const double& p, const double& q);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m, const double& n, const double& o, 
		    const double& p, const double& q, const double& r);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m, const double& n, const double& o, 
		    const double& p, const double& q, const double& r, 
		    const double& s);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m, const double& n, const double& o, 
		    const double& p, const double& q, const double& r, 
		    const double& s, const double& t);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m, const double& n, const double& o, 
		    const double& p, const double& q, const double& r, 
		    const double& s, const double& t, const double& u);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m, const double& n, const double& o, 
		    const double& p, const double& q, const double& r, 
		    const double& s, const double& t, const double& u,
		    const double& v);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m, const double& n, const double& o, 
		    const double& p, const double& q, const double& r, 
		    const double& s, const double& t, const double& u,
		    const double& v, const double& w);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m, const double& n, const double& o, 
		    const double& p, const double& q, const double& r, 
		    const double& s, const double& t, const double& u,
		    const double& v, const double& w, const double& x);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m, const double& n, const double& o, 
		    const double& p, const double& q, const double& r, 
		    const double& s, const double& t, const double& u,
		    const double& v, const double& w, const double& x,
		    const double& y);
   friend Matrix c (const double& a, const double& b, const double& c, 
		    const double& d, const double& e, const double& f,
		    const double& g, const double& h, const double& i,
		    const double& j, const double& k, const double& l,
		    const double& m, const double& n, const double& o, 
		    const double& p, const double& q, const double& r, 
		    const double& s, const double& t, const double& u,
		    const double& v, const double& w, const double& x,
		    const double& y, const double& z);
  


/**********************************************************************/
//  FUNCTION: Transpose  - computes the transpose of a Matrix
//  Input:  Matrix object
//  Return:  Matrix object 
//  Usage:  t(Matrix & old_matrix)
//  Errors:  
//  Dependencies:  none

   friend Matrix t (const Matrix & old_matrix);

/**********************************************************************/
//  FUNCTION: Ones - creates a Matrix of ones
//  Input:  int rows, int cols
//  Return:  Matrix object 
//  Usage:  ones(int rows, int cols)
//  Errors:  0018
//  Dependencies:  none

   static Matrix ones (const int& rows, const int& cols);

/**********************************************************************/
//  FUNCTION: Zeros - creates a Matrix of zeros
//  Input:  int rows, int cols
//  Return:  Matrix object 
//  Usage:  zeros(int rows, int cols)
//  Errors:  0018
//  Dependencies:  none

   static Matrix zeros (const int& rows, const int& cols);

/**********************************************************************/
//  FUNCTION: Eye - creates an Identity Matrix of size k x k
//  Input:  integer k
//  Return:  Matrix object 
//  Usage:  eye(int k)
//  Errors:  
//  Dependencies:  none

   friend Matrix eye (const int& k);

/**********************************************************************/
//  FUNCTION: Seqa - creates a vector additive sequence Matrix (size x 1)
//  Input:  double start, double incr (increment), int size
//  Return:  Matrix object 
//  Usage:  seqa( double, double, int)
//  Errors:  
//  Dependencies:  none

   static Matrix seqa (const double& start, const double& incr,
		       const int& size);

/**********************************************************************/
//  FUNCTION: sort - sorts all elements of a Matrix using shellsort
//                   and places sorted elements in Matrix the same
//                   size as the original. DOES NOT SORT COLUMN BY COLUMN
//  Input:  Matrix A
//  Return:  Matrix object 
//  Usage:  sort(Matrix)
//  Errors:  
//  Dependencies:  none
//  Notes: from Sedgewick, 1992, p. 109 with modifications (0 indexing)

   friend Matrix sort (const Matrix& A);


/**********************************************************************/
//  FUNCTION: sortc - sorts all columns of a Matrix using shellsort
//  Input:  Matrix A
//  Return:  Matrix object 
//  Usage:  sort(Matrix)
//  Errors:  
//  Dependencies:  none
//  Notes: from Sedgewick, 1992, p. 109 with modifications (0 indexing)

   friend Matrix sortc (const Matrix& A);


/**********************************************************************/
//  FUNCTION: Cholesky - Cholesky Decomposition of a Sym. Pos-Def. Matrix
//  Input:  Matrix A
//  Return:  Matrix object (the L Matrix s.t. L*L'=A)
//  Usage:  cholesky (Matrix A)
//  Errors:  0019-0020
//  Dependencies:  none


   friend Matrix cholesky (const Matrix & A);

/**********************************************************************/
//  FUNCTION: Chol_solve - Solves Ax=b for x via backsubstitution using
//  Cholesky Decomposition  (NOTE: function is overloaded)
//  A must be symmetric and positive definite
//  Input:  Matrix A, Matrix b
//  Return:  Matrix x  (solution to Ax=b)
//  Usage:  chol_solve( Matrix A, Matrix b)
//  Errors:  0021
//  Dependencies:  cholesky()

   friend Matrix chol_solve (const Matrix & A, const Matrix & b);

/**********************************************************************/
//  FUNCTION: Chol_solve - Solves Ax=b for x via backsubstitution using
//  Cholesky Decomposition. This function takes in the lower
//  triangular L as input and does not depend upon cholesky()
//  A must be symmetric and positive definite
//  NOTE: function is overloaded
//  Input:  Matrix A, Matrix b, Matrix L
//  Return:  Matrix x  (solution to Ax=b)
//  Usage:  chol_solve( Matrix A, Matrix b, Matrix& L)
//  Errors:  0022
//  Dependencies:  none

   friend Matrix chol_solve (const Matrix & A, const Matrix & b,
    	    const Matrix & L);

/**********************************************************************/
//  FUNCTION: Invpd - Calculates the inverse of a Sym. Pos. Def. Matrix
//  (NOTE: function is overloaded)
//  Input:  Matrix A 
//  Return:  Matrix A^(-1)
//  Usage:  invpd(Matrix A)
//  Errors:  0023-0024
//  Dependencies:  none

   friend Matrix invpd (const Matrix & A);

/**********************************************************************/
//  FUNCTION: Invpd - Calculates the inverse of a Sym. Pos. Def. Matrix
//  (NOTE: function is overloaded)
//  Input:  Matrix A, Matrix L
//  Return:  Matrix A^(-1)
//  Usage:  invpd(Matrix A, Matrix L)
//  Errors: 
//  Dependencies:  none

   friend Matrix invpd (const Matrix & A, const Matrix & L);

/**********************************************************************/
//  FUNCTION: Lu_decomp - Calculates the LU Decomposition of a square
//  Matrix
//  Input:  Matrix A, Matrix L, Matrix U, Matrix perm_vec
//  Return:  integer (0 signifies success)
//  Usage:  lu_decomp(Matrix& A, Matrix& L, Matrix& U, Matrix& perm_vec)
//  Errors:  0025-0026
//  Dependencies:  fabs()

   friend int lu_decomp(const Matrix& A, Matrix& L, Matrix& U, 
           Matrix& perm_vec);  

/**********************************************************************/
//  FUNCTION: Lu_solve - Solve Ax=b for x via forward and
//  backsubstitution using the LU Decomp of Matrix A
//  NOTE: This function is overloaded
//  Input:  Matrix A, Matrix b 
//  Return:  Matrix x
//  Usage:  lu_solve(Matrix& A, Matrix& b)
//  Errors:  0027-0030
//  Dependencies:  fabs()

   friend Matrix lu_solve(const Matrix& A, const Matrix& b);

/**********************************************************************/
//  FUNCTION: Lu_solve - Solve Ax=b for x via forward and
//  backsubstitution using the LU Decomp of Matrix A
//  NOTE: This function is overloaded
//  Input:  Matrix A, Matrix b, Matrix L, Matrix U, Matrix p 
//  Return:  Matrix x
//  Usage:  lu_solve(Matrix& A, Matrix& b, Matrix& L, Matrix& U,
//    		Matrix& p)
//  Errors:  0031-0035
//  Dependencies:  fabs()

   friend Matrix lu_solve(const Matrix& A, const Matrix& b, 
    const Matrix& L, const Matrix& U, const Matrix& p);

/**********************************************************************/
//  FUNCTION: Row_interchange - Interchanges the rows of A with those
//  in vector p and returns the modified Matrix.
//  Input:  Matrix A, Matrix p 
//  Return:  Matrix object
//  Usage:  row_interchange(Matrix& A, Matrix& p)
//  Errors:  0036-0037
//  Dependencies:  

   friend Matrix row_interchange(const Matrix& A, const Matrix& pp);

/**********************************************************************/
//  FUNCTION: Inv - Calculate the Inverse of a square Matrix A via
//  LU decomposition
//  Input:  Matrix A 
//  Return:  Matrix object
//  Usage:  inv(Matrix& A)
//  Errors:  0038-0039
//  Dependencies:  row_interchange(), det()

  friend Matrix inv(const Matrix& A);

/**********************************************************************/
//  FUNCTION: Det - Calculates the determinant of Matrix A via
//  LU Decomposition
//  Input:  Matrix A 
//  Return:  double determinant
//  Usage:  det( Matrix& A)
//  Errors:  0040
//  Dependencies:  

   friend double det(const Matrix& A);

/**********************************************************************/
//  FUNCTION: Cbind - Column bind 2 matrices
//  Input:  Matrix A, Matrix B 
//  Return:  Matrix object
//  Usage:  cbind(Matrix& A, Matrix& B)
//  Errors:  0041
//  Dependencies:  

   friend Matrix cbind (const Matrix & A, const Matrix & B);

/**********************************************************************/
//  FUNCTION: Rbind - Row bind 2 matrices
//  Input:  Matrix A, Matrix B 
//  Return:  Matrix object
//  Usage:  rbind(Matrix& A, Matrix& B)
//  Errors:  0042
//  Dependencies:  

  friend Matrix rbind (const Matrix & A, const Matrix & B);

//  ACCESSOR: Rows - Returns the number of rows in a Matrix
//  ACCESSOR: Cols - Returns the number of columns in a Matrix
//  ACCESSOR: Size - Returns the size of a Matrix (size = rows x cols)
/**********************************************************************/
//  ACCESSOR: Rows - Returns the number of rows in a Matrix
//  Input:  Matrix A
//  Return:  integer
//  Usage:  rows(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend inline int rows (const Matrix & A) {
       return (A.rowsize);
   }

/**********************************************************************/
//  ACCESSOR: Cols - Returns the number of columns in a Matrix
//  Input:  Matrix A
//  Return:  integer
//  Usage:  cols(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend inline int cols (const Matrix & A) {
      return (A.colsize);
   }

/**********************************************************************/
//  ACCESSOR: Size - Returns the size of a Matrix (size = rows x cols)
//  Input:  Matrix A
//  Return:  integer
//  Usage:  size(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend inline int size (const Matrix & A) {
     	return (A.size);
   }

/**********************************************************************/
//  FUNCTION: Sumc - Calculate the sum of each column of a Matrix
//  Input:  Matrix A
//  Return:  Matrix object (row vector)
//  Usage:  sumc(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend Matrix sumc (const Matrix & A);

/**********************************************************************/
//  FUNCTION: Prodc - Calculate the product of each column of a Matrix
//  Input:  Matrix A
//  Return:  Matrix object (row vector)
//  Usage:  prodc(Matrix& A)
//  Errors:  
//  Dependencies:  

  friend Matrix prodc (const Matrix& A);

/**********************************************************************/
//  FUNCTION: Meanc - Calculate the mean of each column of a Matrix
//  Input:  Matrix A
//  Return:  Matrix object (row vector)
//  Usage:  meanc(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend Matrix meanc (const Matrix & A);

/**********************************************************************/
//  FUNCTION: Varc - Calculate the variance of each Matrix column
//  Input:  Matrix A
//  Return:  Matrix object (row vector)
//  Usage:  varc(Matrix& A)
//  Errors:  
//  Dependencies:  meanc()

   friend Matrix varc (const Matrix & A);

/**********************************************************************/
//  FUNCTION: Stdc - Calculate the std deviation of each Matrix column
//  Input:  Matrix A
//  Return:  Matrix object (row vector)
//  Usage:  stdc(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend Matrix stdc (const Matrix & A);

/**********************************************************************/
//  FUNCTION: Sqrt - Calculate the sqrt of each element of a Matrix
//  Input:  Matrix A
//  Return:  Matrix object
//  Usage:  sqrt(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend Matrix sqrt (const Matrix & A);

/**********************************************************************/
//  FUNCTION: Fabs - Calculate the absolute value of each Matrix element
//  Input:  Matrix A
//  Return:  Matrix object
//  Usage:  fabs(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend Matrix fabs (const Matrix & A);

/**********************************************************************/
//  FUNCTION: Exp - Calculate the value of e^x for each  individual
//  Matrix element
//  Input:  Matrix A
//  Return:  Matrix object
//  Usage:  exp(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend Matrix exp(const Matrix& A);

/**********************************************************************/
//  FUNCTION: Log - Calculate the natural log of each Matrix element
//  Input:  Matrix A
//  Return:  Matrix object
//  Usage:  log(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend Matrix log(const Matrix& A);

/**********************************************************************/
//  FUNCTION: Log10 - Calculate the Base 10 Log of each Matrix element
//  Input:  Matrix A
//  Return:  Matrix object
//  Usage:  log10(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend Matrix log10(const Matrix& A);

/**********************************************************************/
//  FUNCTION: Pow - Raise each Matrix element to a specified power
//  Input:  Matrix A, double e
//  Return:  Matrix object
//  Usage:  pow(Matrix& A. double& e)
//  Errors:  
//  Dependencies:  

   friend Matrix pow(const Matrix& A, const double& e);

/**********************************************************************/
//  FUNCTION: Max - Calculates the maximum element in a Matrix
//  Input:  Matrix A
//  Return:  double max
//  Usage:  max(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend double max (const Matrix & A);

/**********************************************************************/
//  FUNCTION: Min - Calculates the minimum element in a Matrix
//  Input:  Matrix A
//  Return:  double min
//  Usage:  min(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend double min (const Matrix & A);

/**********************************************************************/
//  FUNCTION: Maxc - Calculates the maximum of each Matrix column
//  Input:  Matrix A
//  Return:  Matrix object (row vector)
//  Usage:  maxc(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend Matrix maxc (const Matrix & A);

/**********************************************************************/
//  FUNCTION: Minc - Calculates the minimum of each Matrix column
//  Input:  Matrix A
//  Return:  Matrix object (row vector)
//  Usage:  minc(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend Matrix minc (const Matrix & A);

/**********************************************************************/
//  FUNCTION: Maxindc - Finds the index of the max of each Matrix column
//  Input:  Matrix A
//  Return:  Matrix object (row vector)
//  Usage:  maxindc(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend Matrix maxindc(const Matrix& A);
  
/**********************************************************************/
//  FUNCTION: Minindc - Finds the index of the min of each Matrix column
//  Input:  Matrix A
//  Return:  Matrix object (row vector)
//  Usage:  minindc(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend Matrix minindc(const Matrix& A);

/**********************************************************************/
//  FUNCTION: Order - Calculates the order of each element in a Matrix
//  Input:  Matrix A (must be column vector)
//  Return:  Matrix object (column vector)
//  Usage:  order(Matrix& A)
//  Errors:  0043
//  Dependencies:  sumc()

   friend Matrix order(const Matrix& A);

/**********************************************************************/
//  FUNCTION: Selif - Selects all the rows of Matrix A for which
//  binary column vector e has an element equal to 1
//  Input:  Matrix A, Matrix e (binary data)
//  Return:  Matrix object
//  Usage:  selif(Matrix& A, MAtrix& e)
//  Errors:  0044-0046
//  Dependencies:  

   friend Matrix selif(const Matrix& A, const Matrix& e);

/**********************************************************************/
//  FUNCTION: Unique - Finds unique elements in a Matrix
//  Input:  Matrix A
//  Return:  Matrix object (column vector)
//  Usage:  unique(Matrix& A)
//  Errors:  
//  Dependencies:  

   friend Matrix unique(const Matrix& A);

/**********************************************************************/
//  FUNCTION: Vecr - Turn Matrix into Column vector by stacking rows
//  NOTE: Vecr is much faster than Vecc
//  Input:  Matrix A
//  Return:  Matrix object (column vector)
//  Usage:  vecr(Matrix& A)
//  Errors:  
//  Dependencies:

   inline friend Matrix vecr(const Matrix& A) {
	   Matrix temp = Matrix(A.data, A.size, 1);
	   return temp;
   }

/**********************************************************************/
//  FUNCTION: Vecc - Turn Matrix into Column vector by stacking columns
//  NOTE: Vecr is much faster than Vecc
//  Input:  Matrix A
//  Return:  Matrix object (column vector)
//  Usage:  vecc(Matrix& A)
//  Errors:  
//  Dependencies:

   friend Matrix vecc(const Matrix& A);

/**********************************************************************/
//  FUNCTION: Reshape - Reshapes a row major order Matrix or Vector
//  Input:  Matrix A, int rows, int columns
//  Return:  Matrix object (same size, but different # of rows/columns
//  Usage:  reshape(Matrix& A, int rows, int cols)
//  Errors:  0047
//  Dependencies:

   friend Matrix reshape(const Matrix& A, const int r, const int c);

/**********************************************************************/
//  FUNCTION: Vech - Make vector out of unique elements of a symmetric
//  Matrix
//  Input:  Matrix A
//  Return:  Matrix object
//  Usage:  vech(Matrix& A)
//  Errors:  0048
//  Dependencies:

   friend Matrix vech(const Matrix& A);

/**********************************************************************/
//  FUNCTION: Xpnd - Get symmetric Matrix B back from A = vech(B)
//  Input:  Matrix A
//  Return:  Matrix object (symmetric)
//  Usage:  xpnd(Matrix& A)
//  Errors:  0049
//  Dependencies: fmod()

//   friend Matrix xpnd(const Matrix& A);

/**********************************************************************/
//  FUNCTION: Diag - get the diagonal of a Matrix
//  Input:  Matrix A
//  Return:  Matrix object (column vector)
//  Usage:  diag(Matrix& A)
//  Errors:  0050
//  Dependencies:

   friend Matrix diag(const Matrix& A);

/**********************************************************************/
//  FUNCTION: Gaxpy - Fast calculation of A*B + C
//  Input:  Matrix A, Matrix B, Matrix C
//  Return:  Matrix object (result of calculation)
//  Usage:  gaxpy(Matrix& A, Matrix& B, Matrix& C)
//  Errors:  0051-0054
//  Dependencies:

   friend Matrix gaxpy(const Matrix& A, const Matrix& B, const Matrix& C);

/**********************************************************************/
//  FUNCTION: Crossprod - Fast calculation of A'A
//  Input:  Matrix A
//  Return:  Matrix object (result of calculation)
//  Usage:  crossprod(Matrix& A)
//  Errors:  
//  Dependencies:

   friend Matrix crossprod(const Matrix& A);

   friend Matrix crossprod2(const Matrix& A);


/************************    OPERATORS     ****************************/
/**********************************************************************/
//  OPERATOR: Addition
//  NOTE: This operator is overloaded
//  Input:  Matrix A, Matrix B
//  Return:  Matrix object (result of calculation)
//  Usage:  A+B
//  Errors:  0055
//  Dependencies:

   friend Matrix operator + (const Matrix & A, const Matrix & B);

/**********************************************************************/
//  OPERATOR: Addition
//  NOTE: This operator is overloaded
//  Input:  Matrix A, double b
//  Return:  Matrix object (result of calculation)
//  Usage:  A+b
//  Errors:  
//  Dependencies:

   friend Matrix operator + (const Matrix & A, const double &b);

/**********************************************************************/
//  OPERATOR: Addition
//  NOTE: This operator is overloaded
//  Input:  double a, Matrix B
//  Return:  Matrix object (result of calculation)
//  Usage:  a+B
//  Errors:  
//  Dependencies:

   friend Matrix operator + (const double &a, const Matrix & B);

/**********************************************************************/
//  OPERATOR: Subtraction
//  NOTE: This operator is overloaded
//  Input:  Matrix A, Matrix B
//  Return:  Matrix object (result of calculation)
//  Usage:  A-B
//  Errors:  0056
//  Dependencies:

   friend Matrix operator - (const Matrix & A, const Matrix & B);

/**********************************************************************/
//  OPERATOR: Subtraction
//  NOTE: This operator is overloaded
//  Input:  Matrix A, double b
//  Return:  Matrix object (result of calculation)
//  Usage:  A-b
//  Errors:  
//  Dependencies:

   friend Matrix operator - (const Matrix & A, const double &b);

/**********************************************************************/
//  OPERATOR: Subtraction
//  NOTE: This operator is overloaded
//  Input:  double a, Matrix B
//  Return:  Matrix object (result of calculation)
//  Usage:  a-B
//  Errors:  
//  Dependencies:

   friend Matrix operator - (const double &a, const Matrix & B);

/**********************************************************************/
//  OPERATOR: Multiplication
//  NOTE: This operator is overloaded
//  Input:  Matrix A, Matrix B
//  Return:  Matrix object (result of calculation)
//  Usage:  A*B
//  Errors:  0057
//  Dependencies:

  friend Matrix operator * (const Matrix & A, const Matrix & B);

/**********************************************************************/
//  OPERATOR: Multiplication
//  NOTE: This operator is overloaded
//  Input:  Matrix A, double b
//  Return:  Matrix object (result of calculation)
//  Usage:  A*b
//  Errors:  
//  Dependencies:

  friend Matrix operator * (const Matrix & A, const double & b);

/**********************************************************************/
//  OPERATOR: Multiplication
//  NOTE: This operator is overloaded
//  Input:  double a, Matrix B
//  Return:  Matrix object (result of calculation)
//  Usage:  a*B
//  Errors:  
//  Dependencies:

  friend Matrix operator * (const double & a, const Matrix & B);

/**********************************************************************/
//  OPERATOR: Division
//  NOTE: This operator is overloaded
//  Input:  Matrix A, Matrix B
//  Return:  Matrix object (result of calculation)
//  Usage:  A/B
//  Errors:  
//  Dependencies:

   friend Matrix operator / (const Matrix& A, const Matrix& B);

/**********************************************************************/
//  OPERATOR: Division
//  NOTE: This operator is overloaded
//  Input:  Matrix A, double b
//  Return:  Matrix object (result of calculation)
//  Usage:  A/b
//  Errors:  0058
//  Dependencies:

   friend Matrix operator / (const Matrix& A, const double& b);


/**********************************************************************/
//  OPERATOR: Division
//  NOTE: This operator is overloaded
//  Input:  double a, Matrix B
//  Return:  Matrix object (result of calculation)
//  Usage:  a/B
//  Errors:  
//  Dependencies:

   friend Matrix operator / (const double& a, const Matrix& B);

/**********************************************************************/
//  OPERATOR: Kronecker multiplication
//  Input:  Matrix A, Matrix B
//  Return:  Matrix object (result of calculation)
//  Usage:  A%B
//  Errors:  
//  Dependencies:

   friend Matrix operator % (const Matrix& A, const Matrix& B);

/**********************************************************************/
//  OPERATOR: Equality
//  Input:  Matrix A, Matrix B
//  Return:  Integer (False=0, True=1)
//  Usage:  A==B
//  Errors:  
//  Dependencies:

   friend int operator == (const Matrix& A, const Matrix & B);

/**********************************************************************/
//  OPERATOR: Inequality
//  Input:  Matrix A, Matrix B
//  Return:  Integer (False=0, True=1)
//  Usage:  A!=B
//  Errors:  
//  Dependencies:

   friend int operator != (const Matrix& A, const Matrix & B);

/**********************************************************************/
//  OPERATOR: Element-by-element Greater Than
//  NOTE: This operator is overloaded
//  Input:  Matrix A, Matrix B
//  Return:  Matrix of ints (False=0, True=1)
//  Usage:  A>>B
//  Errors:  0059-0060
//  Dependencies:

   friend Matrix operator >> (const Matrix& A, const Matrix& B);

/**********************************************************************/
//  OPERATOR: Element-by-scalar Greater Than
//  NOTE: This operator is overloaded
//  Input:  Matrix A, double b
//  Return:  Matrix of ints (False=0, True=1)
//  Usage:  A>>b
//  Errors:  
//  Dependencies:

   friend Matrix operator >> (const Matrix& A, const double& b);

/**********************************************************************/
//  OPERATOR: Element-by-element Less Than
//  NOTE: This operator is overloaded
//  Input:  Matrix A, Matrix B
//  Return:  Matrix of ints (False=0, True=1)
//  Usage:  A>>B
//  Errors:  0061-0062
//  Dependencies:

   friend Matrix operator << (const Matrix& A, const Matrix& B);

/**********************************************************************/
//  OPERATOR: Element-by-scalar Less Than
//  NOTE: This operator is overloaded
//  Input:  Matrix A, double b
//  Return:  Matrix of ints (False=0, True=1)
//  Usage:  A<<b
//  Errors:  
//  Dependencies:

   friend Matrix operator << (const Matrix& A, const double& b);


/**********************************************************************/
//  OPERATOR: Element-by-element Equality
//  NOTE: This operator is overloaded
//  Input:  Matrix A, Matrix B
//  Return:  Matrix of ints (False=0, True=1)
//  Usage:  A^=B
//  Errors:  0063-0064
//  Dependencies:

   friend Matrix operator ^= (const Matrix& A, const Matrix& B);

/**********************************************************************/
//  OPERATOR: Element-by-scalar Equality
//  NOTE: This operator is overloaded
//  Input:  Matrix A, double b
//  Return:  Matrix of ints (False=0, True=1)
//  Usage:  A^=b
//  Errors:  
//  Dependencies:

   friend Matrix operator ^= (const Matrix& A, const double& b);

};

} // namespace declaration

#endif /* SCYTHE_DOUBLE_MATRIX_H */
