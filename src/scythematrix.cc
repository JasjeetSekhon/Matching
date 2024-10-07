/* Edited by Jasjeet S. Sekhon <jasjeet_sekhon@berkeley.edu> */
/* https://www.jsekhon.com                            */
/*                                                           */
/* April 26, 2013                                            */
/* get rid of friend-injection for ones, zero, seqa          */
/* January 9, 2012                                           */
/* May 9, 2010: Solaris compatibility issues                 */
/* October 24, 2006                                          */
/* May 25, 2006                                              */
/* __NATE__ additions by Nate Begeman (Apple)                */

//
// Memeber function definitions for the Scythe_Double_Matrix.h
// header file.  These functions make up the Matrix class as used
// in the Scythe project.  
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

#ifndef SCYTHE_DOUBLE_MATRIX_CC
#define SCYTHE_DOUBLE_MATRIX_CC

#include "scythematrix.h"

/*! \class Matrix
 *  \brief The Matrix Class that is the cornerstone of the Scythe 
 *	Statistical Library
 
 *  The Matrix class contains functions that modify the Matrix object 
 *	in a variety of ways.
 */
// Avoid NameSpace Pollution
namespace SCYTHE {
    using namespace std;
    
/********************     BASIC FUNCTIONS     *************************/
/**********************************************************************/
    
  /*!
  * \brief CONSTRUCTOR:  Matrix  - Creates Matrix obj. with specificed 
   * rows and columns, sets each element = 0
   * \param rows An integer that reflects the number of rows.
   * \param cols An integer that reflects the numbef of columns.
   */
  // #ifdef __NATE__
  Matrix::Matrix (const int& rows, const int& cols) {
    if (rows < 1 || cols < 1) {
      Rf_error("Improper row or column dimension in Matrix constructor");
    }
    rowsize = rows;    // assign Matrix rowsize
    colsize = cols;    // assign Matrix colsize 
    size = rows * cols;    // assign Matrix size
    //data = new double[size];
    data = (double *) malloc(size * sizeof(double));
    memset(data, 0, size * sizeof(double));
  }


/*!
 * \brief CONSTRUCTOR:  Matrix - Creates Matrix object filled with data 
 * from array
 * \param inputarray An array containing the data for the Matrix in 
 * row major order.
 * \param rows An integer reflecting the number of rows.
 * \param cols An integer that reflects the numbef of columns.
 */
    Matrix::Matrix (const double *inputarray, const int& rows, const int& cols)
    {
	if (rows < 1 || cols < 1) {
	  Rf_error("Improper row or column dimension in Matrix constructor");
	}
	rowsize = rows;    // assign Matrix rowsize
	colsize = cols;    // assign Matrix colsize
	size = rows * cols;    // assign Matrix size
	// data = new double[size];
	data = (double *) malloc(size * sizeof(double));
	
	/*
	for (int i = 0; i < size; ++i) {
	    data[i] = inputarray[i];
	}
	*/
	memcpy(data, inputarray, size*sizeof(double));
    }
    
/*! 
 * \brief CONSTRUCTOR:  Matrix - Creates Matrix object from old Matrix
 * \param old_Matrix A Matrix object that contains the data from 
 * another Matrix.
 */
    Matrix::Matrix (const Matrix & old_Matrix)
    {
	rowsize = old_Matrix.rowsize;  // assign Matrix rowsize
	colsize = old_Matrix.colsize;  // assign Matrix colsize
	size = old_Matrix.size;  // assign Matrix size
	//data = new double[size];
	data = (double *) malloc(size * sizeof(double));
	
	/*
	for (int i = 0; i < size; ++i) {
	    data[i] = old_Matrix.data[i];
	}
	*/
	memcpy(data, old_Matrix.data, size*sizeof(double));
    }
    
/*!
 * \brief OPERATOR:  Matrix operator '=' - Allows for the copying of a 
 * Matrix
 * \param B A Matrix from which the data will be copied.
 * \return A Matrix object.
 */
    
    Matrix & 
    Matrix::operator = (const Matrix & B)
    {
	rowsize = B.rowsize;
	colsize = B.colsize;
	size = B.size;
	
	// delete[]data;
	free(data);
	//data = new double[size];
	data = (double *) malloc(size * sizeof(double));
	memcpy(data, B.data, size * sizeof(double));
	
	return *this;
    }
    
    
    
    
/*!
 * \brief OPERATOR: Matrix operator () - Retrieves all Matrix elements 
 * in row \a i.  Please Note: Indexing starts at 0.
 * \param i a constant integer referring to the number of the 
 * row to be retrived.
 * \param a a constant all_elements 
 *  \return Matrix all elements in row \a i.
 */
    
    Matrix 
    Matrix::operator () (const int& i, const all_elements& a)
    {
	if (i >= rowsize || i < 0) {
	  Rf_error("Index out of range in () operator");
	}
	
	int newrowsize = 1;
	int newcolsize = colsize;
	Matrix newdata(newrowsize, newcolsize);

	/*
	for (int j = 0; j < newcolsize; ++j) {
	    newdata.data[j] = data[j + i*colsize];
	}
	*/
	memcpy(newdata.data, data+i*colsize, newcolsize*sizeof(double));
	return newdata;
    }
    
    
    
/*! 
 * \brief OPERATOR: Matrix operator () - Retrieves all Matrix elements 
 * in column \a j. Please Note: Indexing starts at 0.
 * \param _ a constant of type all_elements.
 * \param j a constant integer referring to the number of the 
 * column to be retrived.
 *  \return Matrix all elements in column \a j.
 */
    
    Matrix 
    Matrix::operator () (const all_elements& a, const int& j)
    {
	if (j >= colsize || j < 0) {
	  Rf_error("Index out of range in () operator");
	}
	
	int newrowsize = rowsize;
	int newcolsize = 1;
	Matrix newdata(newrowsize, newcolsize);
	
	for (int i = 0; i < newrowsize; ++i) {
	    newdata.data[i] = data[j + i*colsize];
	}
	
	return newdata;
    }
    
    
    
/*!
 * \brief OPERATOR: Matrix operator () - Retrieves Matrix elements 
 * in row \a i.  Please Note: Indexing starts at 0.
 * \param i a constant integer referring to the number of the 
 * row to be retrived.
 * \param J a constant reference to the Matrix from which the 
 * data will be extracted.
 *  \return Matrix elements in row \a i.
 */
    
    Matrix
    Matrix::operator () (const int& i, const Matrix& J)
    {
	
	if (i >= rowsize || i < 0) {
	  Rf_error("Index out of range in () operator");
	}
	
	if (J.colsize != 1 && J.rowsize != 1) {
	  Rf_error("Either rows or cols of J != 1 in () operator");
	}
	
	int newrowsize = 1;
	int newcolsize = J.size;
	Matrix newdata(newrowsize, newcolsize);

	/*
	for (int j = 0; j < newcolsize; ++j) {
	    int index = static_cast < int >(J.data[j]);
	    if (index >= colsize || index < 0) {
	        Rf_error("Index out of range in () operator");
	    }
	    index = index + i * colsize;
	    newdata.data[j] = data[index];
	}
	*/
	memcpy(newdata.data, data+i*colsize, newcolsize*sizeof(double));
	
	return newdata;
    }
    
/*! 
 * \brief OPERATOR: Matrix operator () - Retrieves Matrix elements 
 * in column \a j. Please Note: Indexing starts at 0.
 * \param I a constant reference to the Matrix from which the data 
 * will be extracted.
 * \param j a constant integer referring to the number of the 
 * column to be retrived.
 *  \return Matrix elements in column \a j.
 */
    Matrix
    Matrix::operator () (const Matrix& I, const int& j)
    {
	
	if (j >= colsize || j < 0) {
	  Rf_error("Index out of range in () operator");
	}
	
	if (I.colsize != 1 && I.rowsize != 1) {
	  Rf_error("Either rows or cols of I != 1 in () operator");
	}
	
	int newrowsize = I.size;
	int newcolsize = 1;
	Matrix newdata(newrowsize, newcolsize);
	
	for (int i = 0; i < newrowsize; ++i) {
	    int index = static_cast < int >(I.data[i]);
	    if (index >= rowsize || index < 0) {
	      Rf_error("Index out of range in () operator");
	    }
	    index = j + index * colsize;
	    newdata.data[i] = data[index];
	}

	return newdata;
    }
    
    
    
/*! 
 * \brief OPERATOR: Matrix operator () - Extracts submatrix 
 * from existing matrix
 
 * OPERATOR: Matrix operator () - Extracts submatrix from existing matrix.  
 * Get elements \a i,\a j from a Matrix where \a i in \a I and \a j 
 * in \a J.  Indexing starts at 0.
 * \param I a constant reference to the Matrix from which the data 
 * will be extracted.
 * \param J a constant reference to another Matrix from which the 
 * data will be extracted.
 * \return a new Matrix created from the selected data from the 
 * previous two.
 */
    Matrix 
    Matrix::operator () (const Matrix& I, const Matrix& J){
	if (I.colsize != 1 && I.rowsize != 1) {
	  Rf_error("Either Rows or Cols of I != 1 in () operator");
	}
	if (J.colsize != 1 && J.rowsize != 1) {
	  Rf_error("Either rows or cols of J != 1 in () operator");
	}
	if (I.size > rowsize){
	  Rf_error("size(I) > rowsize of Matrix in Matrix operator ()");
	}
	if (J.size > colsize){
	  Rf_error("size(J) > colsize of Matrix in Matrix operator ()");
	}
	
	int place = 0;
	int indexi, indexj;
	Matrix newdata(I.size,  J.size);

	for (int i = 0; i < I.size; i++) {
	    for (int j = 0; j < J.size; j++) {
		indexi = static_cast < int > (I.data[i]);
		indexj = static_cast < int > (J.data[j]);
		if (indexi >= rowsize || indexi < 0) {
		  Rf_error("Row index out of range in () operator");
		}
		if (indexj >= colsize || indexj < 0) {
		  Rf_error("Column index out of range in () operator");
		}
		newdata.data[place] = data[indexi * colsize + indexj];
		place++;
	    }
	}
	
	return newdata;
    }
    
    
    
/*!
 * \brief Prints the Matrix to the screen
 * \param width a constant integer reflecting the screen width 
 * (in characters) for each element of the Matrix.  This is used 
 * in the C++ \e setw() function.
 * \param prec a constant integer reflecting the decimal precision 
 * of each element in the Matrix.  This is used in the C++ \e 
 * setprecision() function.
 * \return void
 */
/*
    void
    Matrix::print (const int width, const int prec)
    {
	
	int count = 0;
	
	for (int i = 0; i < rowsize; ++i) {
	    if (i > 0) {
		cout << endl;
	    }
	    for (int j = 0; j < colsize; ++j) {
		cout << setw (width) << setprecision (prec) 
		     << data[i * colsize + j] << " ";
		++count;
	    }
	}
	cout << endl << endl;
    }
*/
    
/********************   MORE ADVANCED FUNCTIONS  **********************/
/**********************************************************************/
    
//  FUNCTION: c  - concatenates a sequence of doubles into a Matrix
/*!
 * \brief Concatenates a sequence of doubles into a Matrix.
 * 
 * Concatenates a sequence of doubles into a Matrix.
 * \param a first double to be concatenated.
 * \param b second double to be concatenated.
 * \param ... other doubles (number determined by user (up to 26)).
 * \return A Matrix object, the column vector formed by concatenating the 
 * input doubles.
 */
    Matrix c (const double& a, const double& b){
      Matrix newdata(2, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c){
      Matrix newdata(3, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d){ 
      Matrix newdata(4, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e){ 
      Matrix newdata(5, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f){ 
      Matrix newdata(6, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g){ 
      Matrix newdata(7, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h){ 
      Matrix newdata(8, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i){ 
      Matrix newdata(9, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j){ 
      Matrix newdata(10, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k){ 
      Matrix newdata(11, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l){ 
      Matrix newdata(12, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m){ 
      Matrix newdata(13, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n){ 
      Matrix newdata(14, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	newdata.data[13] = n;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o){ 
      Matrix newdata(15, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	newdata.data[13] = n;
	newdata.data[14] = o;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p){ 
      Matrix newdata(16, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	newdata.data[13] = n;
	newdata.data[14] = o;
	newdata.data[15] = p;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q){ 
      Matrix newdata(17, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	newdata.data[13] = n;
	newdata.data[14] = o;
	newdata.data[15] = p;
	newdata.data[16] = q;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r){ 
      Matrix newdata(18, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	newdata.data[13] = n;
	newdata.data[14] = o;
	newdata.data[15] = p;
	newdata.data[16] = q;
	newdata.data[17] = r;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s){ 
      Matrix newdata(19, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	newdata.data[13] = n;
	newdata.data[14] = o;
	newdata.data[15] = p;
	newdata.data[16] = q;
	newdata.data[17] = r;
	newdata.data[18] = s;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t){ 
      Matrix newdata(20, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	newdata.data[13] = n;
	newdata.data[14] = o;
	newdata.data[15] = p;
	newdata.data[16] = q;
	newdata.data[17] = r;
	newdata.data[18] = s;
	newdata.data[19] = t;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t, const double& u){ 
      Matrix newdata(21, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	newdata.data[13] = n;
	newdata.data[14] = o;
	newdata.data[15] = p;
	newdata.data[16] = q;
	newdata.data[17] = r;
	newdata.data[18] = s;
	newdata.data[19] = t;
	newdata.data[20] = u;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t, const double& u,
	      const double& v){ 
      Matrix newdata(22, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	newdata.data[13] = n;
	newdata.data[14] = o;
	newdata.data[15] = p;
	newdata.data[16] = q;
	newdata.data[17] = r;
	newdata.data[18] = s;
	newdata.data[19] = t;
	newdata.data[20] = u;
	newdata.data[21] = v;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t, const double& u,
	      const double& v, const double& w){ 
      Matrix newdata(23, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	newdata.data[13] = n;
	newdata.data[14] = o;
	newdata.data[15] = p;
	newdata.data[16] = q;
	newdata.data[17] = r;
	newdata.data[18] = s;
	newdata.data[19] = t;
	newdata.data[20] = u;
	newdata.data[21] = v;
	newdata.data[22] = w;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t, const double& u,
	      const double& v, const double& w, const double& x){ 
      Matrix newdata(24, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	newdata.data[13] = n;
	newdata.data[14] = o;
	newdata.data[15] = p;
	newdata.data[16] = q;
	newdata.data[17] = r;
	newdata.data[18] = s;
	newdata.data[19] = t;
	newdata.data[20] = u;
	newdata.data[21] = v;
	newdata.data[22] = w;
	newdata.data[23] = x;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t, const double& u,
	      const double& v, const double& w, const double& x,
	      const double& y){ 
      Matrix newdata(25, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	newdata.data[13] = n;
	newdata.data[14] = o;
	newdata.data[15] = p;
	newdata.data[16] = q;
	newdata.data[17] = r;
	newdata.data[18] = s;
	newdata.data[19] = t;
	newdata.data[20] = u;
	newdata.data[21] = v;
	newdata.data[22] = w;
	newdata.data[23] = x;
	newdata.data[24] = y;
	
	return newdata;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t, const double& u,
	      const double& v, const double& w, const double& x,
	      const double& y, const double& z){ 
      Matrix newdata(26, 1);
	newdata.data[0] = a;
	newdata.data[1] = b;
	newdata.data[2] = c;
	newdata.data[3] = d;
	newdata.data[4] = e;
	newdata.data[5] = f;
	newdata.data[6] = g;
	newdata.data[7] = h;
	newdata.data[8] = i;
	newdata.data[9] = j;
	newdata.data[10] = k;
	newdata.data[11] = l;
	newdata.data[12] = m;
	newdata.data[13] = n;
	newdata.data[14] = o;
	newdata.data[15] = p;
	newdata.data[16] = q;
	newdata.data[17] = r;
	newdata.data[18] = s;
	newdata.data[19] = t;
	newdata.data[20] = u;
	newdata.data[21] = v;
	newdata.data[22] = w;
	newdata.data[23] = x;
	newdata.data[24] = y;
	newdata.data[25] = z;
	
	return newdata;
    }
    
    
    
//  FUNCTION: Transpose  - computes the transpose of a Matrix
/*!
 * \brief Computes the transpose of the given Matrix.
 * 
 * Computes the transpose of the given Matrix.
 * \param old_matrix a constant reference to a Matrix object.  
 * This Matrix will be transposed.
 * \return A Matrix object, the transpose of the input Matrix object.
 */
// #ifdef __NATE__
Matrix t (const Matrix & old_matrix)
{
  int newrowsize = old_matrix.colsize;
  int newcolsize = old_matrix.rowsize;
  Matrix temp(newrowsize, newcolsize);
  for (int i = 0; i < newcolsize; ++i) {
    for (int j = 0; j < newrowsize; ++j) {
      temp.data[i + newcolsize * j] = old_matrix.data[j + newrowsize * i];
    }
  }
  return temp;
}
    
/*!
 * \brief Creates a Matrix of Ones
 *
 * Creates a Matrix filled with Ones, given a specified size.
 * \param rows a constant int reflecting the number of rows in the Matrix.
 * \param cols a constant int reflecting the number of columns in the Matrix.
 * \return a new Matrix filled with 1's.
 */
  Matrix Matrix::ones (const int& rows, const int& cols)
 {
   if (rows < 1 || cols < 1) {
     Rf_error("improper row or column dimension in ones()");
   }
   Matrix newdata(rows, cols);
   int size = rows * cols;
      for (int i = 0; i < size; ++i) {
	newdata.data[i] = 1.0;
      }
      return newdata;
 }

    
/*!
 * \brief Creates a Matrix of Zeros
 *
 * Creates a Matrix filled with Zeros, given a specified size.
 * \param rows a constant int reflecting the number of rows in the Matrix.
 * \param cols a constant int reflecting the number of columns in the Matrix.
 * \return a new Matrix filled with 1's.
 */
// #ifdef __NATE__
  Matrix Matrix::zeros (const int& rows, const int& cols)
{
  if (rows < 1 || cols < 1) {
    Rf_error("Error 0018: improper row or column dimension in ones()");
  }
  Matrix temp(rows, cols); // ctor zeros data
  return temp;
} // end of zeros
    
//  FUNCTION: Eye - creates an Identity Matrix of size k x k
/*! 
 * \brief Creates an Identity Matrix
 *
 * Creates an Identity Matrix of size \a k \a x \a k.
 * \param k a constant integer reflecting the length and width 
 * of the identity matrix.
 * \return the Identity Matrix
 */
    Matrix eye (const int& k)
    {
      Matrix newdata(k, k);
	double hold;
	for (int i = 0; i < k; ++i) {
	    for (int j = 0; j < k; ++j) {
		if (i == j)
		    hold = 1.0;
		else
		    hold = 0.0;
		newdata.data[k * i + j] = hold;
	    }
	}
	return newdata;
    }
    
/*! 
 * \brief Creates a Vector-additive sequence Matrix
 *
 * Creates a Vector-additive sequence Matrix of (\a size x 1)
 * \param start a constant double reflecting the start value of 
 * the first element in the vector.
 * \param incr a double constant reflecting the incremental step 
 * value between each matrix element.
 * \param size a constant integer reflecting the size of the vector.
 * \return a new Matrix (vector).
 */
Matrix Matrix::seqa (const double& start, const double& incr, const int& size)
    {
      Matrix newdata(size, 1);
	double val = start;
	for (int i = 0; i < size; ++i) {
	    newdata.data[i] = val;
	    val += incr;
	}
	return newdata;
    }
    
    
/*! 
 * \brief Sorts all elements of a Matrix (not column by column) using 
 * shellsort
 *
 * Sorts all elements of a Matrix (not column by column) using shellsort
 * \param A the Matrix to be sorted.
 * \return a new Matrix the same size as the original in 
 * which all elements have been sorted.
 */
    Matrix sort(const Matrix& A){
	int i, j, h;
	double v;
	
	Matrix newdata(A.rowsize, A.colsize);
	for (i = 0; i<A.size; ++i)
	    newdata.data[i] = A.data[i];
	for (h = 1; h <= A.size/9; h = 3*h+1);
	for (; h > 0; h /= 3)
	    for (i = h+1; i <= A.size; i += 1){
		v = newdata.data[i-1]; 
		j = i;
		while (j>h && newdata.data[j-h-1] > v)
		{newdata.data[j-1] = newdata.data[j-h-1]; j -= h;}
		newdata.data[j-1] = v;
	    }
	return newdata;     
    }
    
    
    
/*! 
 * \brief Sorts all columns of a Matrix using shellsort
 *
 * Sorts all columns of a Matrix using shellsort
 * \param A the Matrix to be sorted.
 * \return a new Matrix the same size as the original in 
 * which all elements have been sorted.
 */
    Matrix sortc(const Matrix& A){
	int i, j, h;
	double v;
	
	Matrix newdata(A.rowsize, A.colsize);
	for (i = 0; i<A.size; ++i)
	    newdata.data[i] = A.data[i];
	
	for (int colindex=0; colindex<A.colsize; ++colindex){
	    
	    for (h = 1; h <= A.rowsize/9; h = 3*h+1);
	    for (; h > 0; h /= 3)
		for (i = h+1; i <= A.rowsize; i += 1){
		    v = newdata.data[(i-1)*A.colsize + colindex]; 
		    j = i;
		    while (j>h && newdata.data[(j-h-1)*A.colsize + colindex] > v){
			newdata.data[(j-1)*A.colsize + colindex] = 
			    newdata.data[(j-h-1)*A.colsize + colindex]; 
			j -= h;
		    }
		    newdata.data[(j-1)*A.colsize + colindex] = v;
		}
	    
	}
	return newdata;     
    }
    
    
//  4/29/2001 (KQ)
/*! 
 * \brief Interchanges the rows of \a A with those in vector \a p
 *
 * Interchanges the rows of \a A with those in vector \a p and 
 * returns the modified Matrix. Useful for putting \a A into the form 
 * of its permuted LU factorization.
 * \param A a constant reference to the Matrix \a A.
 * \param pp a constant reference to the Matrix (vector) \a p from 
 * which the interchange row will come from.
 * \return the modified Matrix \a A.
 */
    Matrix row_interchange(const Matrix& A, const Matrix& pp){
	Matrix PA = A;
	Matrix p = pp;
	if (p.colsize != 1){
	  Rf_error("Vector p not a column vector in row_interchange()");
	}
	if ( (p.rowsize +1) != A.rowsize){
	  Rf_error("Matrices A and p not of consistent sizes in row_interchange()");
	}
	
	for (int i=0; i<(A.rowsize-1); ++i){
	    //swap A(i,.) and A(p[i],.)
	    int swap_row = static_cast<int>(p.data[i]);
	    for (int j=0; j<A.colsize; ++j){
		double temp = PA.data[j+A.colsize*i];
		PA.data[j+A.colsize*i] = PA.data[j+A.colsize*swap_row];
		PA.data[j+A.colsize*swap_row] = temp;
	    }
	}
	
	return PA;
    }
    
    
//  4/29/2001 (KQ)
/*! 
 * \brief Calculate the Inverse of a square Matrix \a A
 *
 * Calculate the Inverse of a square Matrix \a A via LU decomposition.
 * \param AA a constant reference to the Matrix \a A to be inverted.
 * \return the inverted Matrix \a A.
 */
    Matrix inv (const Matrix & AA)
    {
	if (AA.rowsize != AA.colsize){
	  Rf_error("Matrix A not square in SCYTHE::inv()");
	}
	
	Matrix b = Matrix (AA.rowsize, 1); 
	Matrix Ainv = Matrix(AA.rowsize, AA.rowsize);    
	
	// step 1 compute the LU factorization 
	// taken from lu_decomp()
	Matrix A = AA;
	Matrix L, U, perm_vec;
	
	if (A.rowsize == 1){
	  L = Matrix::ones(1,1);
	  U = AA;
	  perm_vec = Matrix(1,1);
	} else {
	    L = U = Matrix(A.rowsize, A.rowsize);
	    perm_vec = Matrix(A.rowsize-1,1);
	    
	    for (int k=0; k<A.rowsize-1; ++k){
		int pivot = k;
		// find pivot
		for (int i=k; i<A.rowsize; ++i){
		    if ( std::fabs(A(pivot,k)) < std::fabs(A(i,k)) ) pivot = i;  
		}
		
		if(A(pivot,k) == 0.0){
		  Rf_error("Matrix A is singular in SCYTHE::inv()");
		}
		
		// permute 
		if (k != pivot){
		    for (int i=0; i<A.rowsize; ++i){
			double temp = A(pivot,i);
			A(pivot,i) = A(k,i);
			A(k,i) = temp;
		    }
		}
		perm_vec[k] = pivot;
		
		for (int i = k+1; i<A.rowsize; ++i){
		    A(i,k) = A(i,k)/A(k,k);
		    for (int j = k+1; j <A.rowsize; ++j){
			A(i,j) = A(i,j) - A(i,k)*A(k,j);
		    }
		}
	    }
	    
	    L = A; 
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=i; j<A.rowsize; ++j){
		    U(i,j) = A(i,j);
		    L(i,j) = 0.0;
		    L(i,i) = 1.0;
		}
	    }
	}
	
	// step 2 repeated solving of A*hold = b
	for (int j = 0; j < A.rowsize; ++j) {
	    b.data[j] = 1.0;
	    //Matrix hold = lu_solve(A, b, L, U, p);
	    
	    // step 2.1 solve L*y = Pb via forward substitution
	    Matrix bb = row_interchange(b, perm_vec);
	    Matrix y = Matrix(A.rowsize,1);
	    for (int i=0; i<A.rowsize; ++i){
		double sum = 0.0;
		for (int j=0; j<i; ++j)
		    sum += L.data[j+A.colsize*i]*y.data[j];
		y.data[i] = (bb.data[i] - sum)/L.data[i+A.colsize*i];
	    }
	    
	    // step 2.2 solve U*x = y via backsubstitution
	    Matrix x = Matrix(A.rowsize,1);
	    for (int i=A.rowsize-1; i>=0; --i){
		double sum = 0.0;
		for (int j=i+1; j<A.rowsize; ++j)
		    sum += U.data[j+A.colsize*i] * x.data[j];
		x.data[i] = (y.data[i] - sum)/U.data[i+A.colsize*i];
	    }
	    
	    // step 3 reset b and put the solution in Ainv
	    b.data[j] = 0.0;
	    for (int k=0; k<A.rowsize; ++k)
		Ainv(k,j) = x[k];
	}
	
	return Ainv;
    }
    
    
//  NOTE: LU decomposition algorithm is based on  Algorithm 3.4.1 
//      of Golub and Van Loan 3rd edition, 1996. 
/*! 
 * \brief Calculate the determinant of a square Matrix \a A
 *
 * Calculate the determinant of a square Matrix \a A via LU decomposition.
 * \param AA a constant reference to the Matrix.
 * \return the determinant of \a A.
 */
    double det(const Matrix& AA){
	
	Matrix A = AA;
	
	if (A.rowsize != A.colsize){
	  Rf_error("Matrix A not square in SCYTHE::det()");
	}
	
	if(A.rowsize == 1)
	    return A(0,0);
	
	Matrix L = Matrix(A.rowsize, A.rowsize);
	Matrix U = L;
	double sign = 1.0;
	
	for (int k=0; k<A.rowsize-1; ++k){
	    int pivot = k;
	    // find pivot
	    for (int i=k; i<A.rowsize; ++i){
		if ( A(pivot,k) < std::fabs(A(i,k)) ) pivot = i;  
	    }
	    
	    if(A(pivot,k) == 0.0){
		return 0.0;
	    }
	    
	    // permute 
	    if (k != pivot){
		sign = -1*sign;
		for (int i=k; i<A.rowsize; ++i){
		    double temp = A(pivot,i);
		    A(pivot,i) = A(k,i);
		    A(k,i) = temp;
		}
	    }
	    
	    for (int i = k+1; i<A.rowsize; ++i){
		A(i,k) = A(i,k)/A(k,k);
		for (int j = k+1; j <A.rowsize; ++j){
		    A(i,j) = A(i,j) - A(i,k)*A(k,j);
		}
	    }
	}
	
	double det = 1.0;
	for (int i = 0; i<A.rowsize; ++i)
	    det = det*A(i,i);
	
	return sign*det;
    }
    
//! Column bind 2 matrices
/*!
 * Column bind 2 matrices,\a A and \a B.
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the new Matrix \a A
 */
    Matrix cbind (const Matrix & A, const Matrix & B)
    {
	if (A.rowsize != B.rowsize) {
	  Rf_error("Matrices A and B do not have some number of rows in SCYTHE::cbind()");
	}
	
	int totalcols = A.colsize + B.colsize;
	//double *newdata = new double[A.rowsize * totalcols];
	Matrix newdata(A.rowsize, totalcols);
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		newdata.data[i * totalcols + j] = A.data[i * A.colsize + j];
	    }
	    for (int k = 0; k < B.colsize; ++k) {
		newdata.data[i * totalcols + k + A.colsize] = B.data[i * B.colsize + k];
	    }
	}
	
	return newdata;
    }
    
//! Row bind 2 matrices
/*!
 * Row bind 2 matrices,\a A and \a B.
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the new Matrix \a A.
 */
    Matrix rbind (const Matrix & A, const Matrix & B)
    {
	if (A.colsize != B.colsize) {
	  Rf_error("Matrices A and B do not have some number of cols in SCYTHE::rbind()");
	}
	
	int totalrows = A.rowsize + B.rowsize;
	//double *newdata = new double[totalrows * A.colsize];
	Matrix newdata(totalrows, A.colsize);
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j)  {
		newdata.data[i * A.colsize + j] = A.data[i * A.colsize + j];
	    }
	}
	for (int k = 0; k < B.rowsize; ++k) {
	    for (int j = 0; j < B.colsize; ++j) {
		newdata.data[k * B.colsize + (A.rowsize * A.colsize) + j] =
		    B.data[k * B.colsize + j];
	    }
	}
	
	return newdata;
    }
    
//! Calculate the sum of each column of a Matrix
/*!
 * Calculate the sum of each column of a Matrix.  
 * Returns the row vector of sums.
 * \param A a constant reference to a Matrix \a A.
 * \return the vector of sums for each corresponding column.
 */
// #ifdef __NATE__
Matrix sumc (const Matrix & A)
{
  Matrix temp = Matrix::zeros(1, A.colsize);
  
  for (int i = 0; i < A.rowsize; ++i) {
    double *rowptr = A.data + (A.colsize * i);
    for (int j = 0; j < A.colsize; ++j) {
      temp.data[j] += rowptr[j];
    }
  }
  return temp;
}

    
//! Calculate the product of each column of a Matrix
/*!
 * Calculate the product of each column of a Matrix.  
 * Returns the row vector of products.
 * \param A a constant reference to a Matrix \a A.
 * \return the vector of products for each corresponding column.
 */
    Matrix prodc (const Matrix & A)
    {
      Matrix newdata(1, A.colsize);
	
	for (int i = 0; i < A.colsize; ++i)
	    newdata.data[i] = 1.0;
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		newdata.data[j] = newdata.data[j] * A.data[A.colsize * i + j];
	    }
	}
	
	return newdata;
    }
    
//! Calculate the mean of each column of a Matrix
/*!
 * Calculate the mean of each column of a Matrix.  
 * Returns the row vector of means.
 * \param A a constant reference to a Matrix \a A.
 * \return the vector of means for each corresponding column.
 */
    Matrix meanc (const Matrix & A)
    {
      Matrix newdata(1, A.colsize);
	for (int i = 0; i < A.colsize; ++i)
	    newdata.data[i] = 0.0;
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		newdata.data[j] += A.data[A.colsize * i + j];
	    }
	}
	for (int i = 0; i < A.colsize; ++i) {
	    newdata.data[i] = (1.0 / A.rowsize) * newdata[i];
	}
	
	return newdata;
    }
    
//! Calculate the variances of each column of a Matrix
/*!
 * Calculate the variances of each column of a Matrix.  
 * Returns the row vector of variances.
 * \param A a constant reference to a Matrix \a A.
 * \return the vector of variances for each corresponding column.
 */
    Matrix varc (const Matrix & A)
    {
	Matrix mu = meanc (A);
	Matrix newdata(1, A.colsize);
	for (int i = 0; i < A.colsize; ++i)
	    newdata.data[i] = 0.0;
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		newdata.data[j] += std::pow (mu.data[j] - A.data[A.colsize * i + j], 2) *
		    (1.0 / (A.rowsize-1));
	    }
	}
	
	return newdata;
    }
    
//! Calculate the standard deviation of each column of a Matrix
/*!
 * Calculate the standard deviation of each column of a Matrix.  
 * Returns the row vector of standard deviations.
 * \param A a constant reference to a Matrix \aA.
 * \return the vector of standard deviations for each corresponding column.
 */
    Matrix stdc (const Matrix & A)
    {
	Matrix mu = meanc (A);
	Matrix newdata(1, A.colsize);
	for (int i = 0; i < A.colsize; ++i)
	    newdata.data[i] = 0.0;
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		newdata.data[j] += std::pow (mu.data[j] - A.data[A.colsize * i + j], 2) *
		    (1.0 / (A.rowsize-1));
	    }
	}
	
	for (int i = 0; i < A.colsize; ++i)
	    newdata.data[i] = std::sqrt (newdata[i]);
	
	return newdata;
    }
    
//! Calculate the square root of each element of a Matrix
/*!
 * Calculate the square root of each element of a Matrix.
 * \param A a constant reference to a Matrix \a A.
 * \return the Matrix of square roots.
 */
    Matrix sqrt (const Matrix & A)
    {
      Matrix newdata(A.rowsize, A.colsize);
	
	for (int i = 0; i < A.size; ++i)
	    newdata.data[i] = std::sqrt (A.data[i]);
	
	return newdata;
    }
    
//! Calculate the absolute value of each element of a Matrix
/*!
 * Calculate the absolute value of each element of a Matrix.
 * \param A a constant reference to a Matrix \a A.
 * \return the Matrix of absolute values.
 */
    Matrix fabs (const Matrix & A)
    {
      Matrix newdata(A.rowsize, A.colsize);
	
	for (int i = 0; i < A.size; ++i)
	    newdata.data[i] = std::fabs (A.data[i]);
	
	return newdata;
    }
    
    
//! Calculate the value of \a e^x for each individual Matrix element
/*!
 * Calculate the value of \a e^x for each individual Matrix element.
 * \param A a constant reference to a Matrix \a A.
 * \return the new Matrix of exponentials.
 */
    Matrix exp(const Matrix& A){
      Matrix newdata(A.rowsize, A.colsize);
	
	for (int i = 0; i < A.size; ++i)
	    newdata.data[i] = std::exp (A.data[i]);
	
	return newdata;
    }
    
//! Calculate the natural log of each individual Matrix element
/*!
 * Calculate the natural log of each individual Matrix element.
 * \param A a constant reference to a Matrix \a A.
 * \return the new Matrix of natural logs.
 */
    Matrix log(const Matrix& A){
      Matrix newdata(A.rowsize, A.colsize);
	
	for (int i = 0; i < A.size; ++i)
	    newdata.data[i] = std::log (A.data[i]);
	
	return newdata;
    }
    
//! Calculate the Base 10 Log of each Matrix element
/*!
 * Calculate the Base 10 Log of each Matrix element.
 * \param A a constant reference to a Matrix \a A.
 * \return the new Matrix of Base 10 logs.
 */
    Matrix log10(const Matrix& A){
      Matrix newdata(A.rowsize, A.colsize);

	for (int i = 0; i < A.size; ++i)
	    newdata.data[i] = std::log10(A.data[i]);
	
	return newdata;
    }
    
//! Raises each Matrix element to a specified power
/*!
 * Raises each Matrix element to a specified power.
 * \param A a constant reference to a Matrix \a A.
 * \return the new modified Matrix.
 */
    Matrix pow(const Matrix& A, const double& e){
      Matrix newdata(A.rowsize, A.colsize);
	
	for (int i = 0; i < A.size; ++i)
	    newdata.data[i] = std::pow(A.data[i], e);
	
	return newdata;
    }
    
    
//! Calculates the maximum element in a Matrix
/*!
 * Calculates the maximum element in a Matrix.
 * \param A a constant reference to a Matrix \a A.
 * \return the maximum element (a double).
 */
    double max (const Matrix & A)
    {
	double max = A.data[0];
	for (int i = 1; i < A.size; ++i) {
	    if (A.data[i] > max)
		max = A.data[i];
	}
	return max;
    }
    
//! Calculates the minimum element in a Matrix
/*!
 * Calculates the minimum element in a Matrix.
 * \param A a constant reference to a Matrix \a A.
 * \return the minimum element (a double).
 */
    double min (const Matrix & A)
    {
	double min = A.data[0];
	for (int i = 1; i < A.size; ++i) {
	    if (A.data[i] < min)
		min = A.data[i];
	}
	return min;
    }
    
//! Calculates the maximum of each Matrix column
/*!
 * Calculates the maximum of each Matrix column.
 * \param A a constant reference to a Matrix \a A.
 * \return a Matrix (vector) of the maximum elements.
 */
    Matrix maxc (const Matrix & A)
    {
      Matrix newdata(1, A.colsize);
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		if (i == 0) {
		    newdata.data[j] = A.data[A.colsize * i + j];
		} else if (A.data[A.colsize * i + j] > newdata.data[j]) {
		    newdata.data[j] = A.data[A.colsize * i + j];
		}
	    }
	}
	
	return newdata;
    }
    
//! Calculates the minimum of each Matrix column
/*!
 * Calculates the minimum of each Matrix column.
 * \param A a constant reference to a Matrix \a A.
 * \return a Matrix (vector) of the minimum elements.
 */
    Matrix minc (const Matrix & A)
    {
      Matrix newdata(1, A.colsize);
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		if (i == 0) {
		    newdata.data[j] = A.data[A.colsize * i + j];
		} else if (A.data[A.colsize * i + j] < newdata.data[j]) {
		    newdata.data[j] = A.data[A.colsize * i + j];
		}
	    }
	}
	return newdata;
    }
    
    
//! Finds the index of the maximum of each Matrix column
/*!
 * Finds the index of the maximum of each Matrix column.
 * \param A a constant reference to a Matrix \a A.
 * \return a Matrix (vector) of the index of each maximum element.
 */
    Matrix maxindc(const Matrix& A){
      Matrix newdata(1, A.colsize);
	Matrix maxvec = Matrix(1,A.colsize);
	
	for (int i = 0; i < A.rowsize; ++i){
	    for (int j = 0; j < A.colsize; ++j){
		if (i == 0){
		    maxvec[j] = A.data[A.colsize * i + j]; 
		    newdata.data[j] = 0;
		} else if (A.data[A.colsize * i + j] > maxvec[j]){
		    maxvec[j] = A.data[A.colsize * i + j];
		    newdata.data[j] = i;
		}
	    }
	}
	return newdata;
    }
    
    
//! Finds the index of the minimum of each Matrix column
/*!
 * Finds the index of the minimum of each Matrix column.
 * \param A a constant reference to a Matrix \a A.
 * \return a Matrix (vector) of the index of each minimum element.
 */
    Matrix minindc(const Matrix& A){
        Matrix newdata(1, A.colsize);
	Matrix minvec = Matrix(1,A.colsize);
	
	for (int i = 0; i < A.rowsize; ++i){
	    for (int j = 0; j < A.colsize; ++j){
		if (i == 0){
		    minvec[j] = A.data[A.colsize * i + j]; 
		    newdata.data[j] = 0;
		} else if (A.data[A.colsize * i + j] < minvec[j]){
		    minvec[j] = A.data[A.colsize * i + j];
		    newdata.data[j] = i;
		}
	    }
	}
	return newdata;
    }
    
    
//! Calculates the order of each element in a Matrix
/*!
 * Calculates the order of each element in a Matrix.
 * \param A a constant reference to a Matrix A.
 * \return a Matrix (vector) in which the \e i'th element 
 * gives the order position of the \e i'th element of \a A.
 */
    Matrix order(const Matrix& A){
	if (A.colsize != 1){
	  Rf_error("Matrix A not a column vector in SCYTHE::order()");
	}
	Matrix newdata(A.rowsize, 1);
	
	for (int i=0; i<A.rowsize; ++i){
	    newdata.data[i] = sumc(A << A.data[i])[0];    
	}
	
	return newdata;
    }
    
    
//! Selects all the rows of Matrix \a A for which binary column vector \a e has an element equal to 1
/*!
 * Selects all the rows of Matrix \a A for which binary column 
 * vector \a e has an element equal to 1.
 * \param A a constant reference to a Matrix \a A (n x k).
 * \param e a constant reference to a vector \a e (n x 1).
 * \return a Matrix of all rows of \a A for \a e equal to one.
 */
    Matrix selif(const Matrix& A, const Matrix& e){
	
	// check to see if rowsize matches  
	if (A.rowsize != e.rowsize){
	  Rf_error("Matrices not conformable in SCYTHE::selif()");
	}
	
	// check to see if e is a column vector
	if (e.colsize > 1){
	  Rf_error("Not a column vector in SCYTHE::selif()");
	}
	
	// loop to check if e contains binary data, and count number
	// of output rows 
	int N = 0;
	for (int i=0; i<e.rowsize; ++i){
	    if (e.data[i] != 0 && e.data[i] != 1){
	      Rf_error("Vector contains non binary data in SCYTHE::selif()");
	    }
	    if (e.data[i] == 1) {
		N += 1;
	    }
	}
	
	// declare and form output Matrix
	Matrix newdata(N, A.colsize);
	int count = 0;
	for (int i=0; i<A.rowsize; ++i){
	    if (e.data[i] == 1){
		for (int j=0; j<A.colsize; ++j){
		    newdata.data[count] = A.data[A.colsize*i + j]; 
		    ++count;
		}
	    }
	}
	return newdata;
    }
    
    
//! Finds unique elements in a Matrix
/*!
 * Finds unique elements in a Matrix.
 * \param A a constant reference to a Matrix \a A.
 * \return a Matrix of all unique elements.
 */
    Matrix unique(const Matrix& A){
      // double *newdata = new double[A.size];
      double *newdata = (double *) malloc( (A.size)*sizeof(double));
	
	newdata[0] = A.data[0];
	int count = 1;
	for (int i=1; i<A.size; ++i){
	    int uniq = 1;
	    for (int j=0; j<count; ++j){
		if (newdata[j] == A.data[i]){
		    uniq = 0;
		    break;
		}
	    }
	    if (uniq==1){
		newdata[count] = A.data[i];
		++count;
	    }
	}
	
	Matrix temp = Matrix(newdata, count, 1);
	// delete[] newdata;
	free(newdata);
	return temp;
    }
    
//! Turn Matrix into Column vector by stacking columns
/*!
 * Turn Matrix into Column vector by stacking columns. 
 * NOTE: \e Vecr() is much faster than \e Vecc().
 * \param A a constant reference to a Matrix \a A.
 * \return a Column vector.
 */
    Matrix vecc(const Matrix& A){
	// first transposes the input Matrix's data 
	int newrowsize = A.colsize;
	int newcolsize = A.rowsize;
	Matrix newdata(A.size, 1);
	for (int i = 0; i < newcolsize; ++i) {
	    for (int j = 0; j < newrowsize; ++j) {
		newdata.data[i + newcolsize * j] = A.data[j + newrowsize * i];
	    }
	}
	// then takes this transposed data and vecrs
	return newdata;
    }
    
    
//! Reshapes a row major order Matrix or Vector
/*!
 * Reshapes a row major order Matrix or Vector.
 * \param A a constant reference to a Matrix \a A.
 * \param r a constant integer reflecting the new number of rows.
 * \param c a constant integer reflecting the new number of columns.
 * \return a new Matrix with the specified dimensions.
 */
    Matrix reshape(const Matrix& A, const int r, const int c){
	if (A.size != r*c){
	  Rf_error("Input dimensions to SCYTHE::reshape() not consistent with size of input Matrix");
	}
	Matrix temp = Matrix(A.data, r, c);
	return temp;
    }
    
//! Make vector out of unique elements of a symmetric Matrix
/*!
 * Make vector out of unique elements of a symmetric Matrix.  
 * NOTE: DOES NOT CHECK FOR SYMMETRY!!!
 * \param A a constant reference to a Matrix \a A.
 * \return a Column vector.
 */
    Matrix vech(const Matrix& A){
	if (A.rowsize != A.colsize){
	  Rf_error("Input Matrix not square in SCYTHE::vech()");
	}
	
	int newsize = static_cast<int>(0.5*(A.size - A.rowsize) + A.rowsize);
	Matrix newdata(newsize, 1);
	int count = 0;
	for (int i=0; i<A.rowsize; ++i){
	    for(int j=i; j<A.colsize; ++j){
		newdata.data[count] = A.data[i*A.colsize + j];
		++count;
	    }
	}
	
	return newdata;
    }
    
    
//! Get symmetric Matrix B back from \a A = \a vech(B)
/*!
 * Get symmetric Matrix \a B back from \a A = \a vech(B)
 * \param A a constant reference to a Matrix \a A.
 * \return a Symmetric Matrix.
 */
/* 
    Matrix xpnd(const Matrix& A){
	double newrowsize_d = -.5 + .5*std::sqrt(1+8*A.size);
	if (fmod(newrowsize_d,1.0) != 0.0){
	    Rf_error("Not possible to make square Matrix out of input Matrix to SCYTHE::xpnd()");
	}
	int newrowsize = static_cast<int>(newrowsize_d);
	Matrix newdata(newrowsize, newrowsize);
	int count = 0;
	for (int i=0; i<newrowsize; ++i){
	    for (int j=i; j<newrowsize; ++j){
		newdata.data[i*newrowsize +j] = newdata.data[j*newrowsize + i] = A.data[count];
		++count;
	    }
	}
	return newdata;
    }
*/
    
    
//! Get the diagonal of a Matrix
/*!
 * Get the diagonal of a Matrix.
 * \param A a constant reference to a Matrix \a A.
 * \return a Matrix (vector) containing the diagonal of \a A.
 */
    Matrix diag(const Matrix& A){
	if (A.rowsize != A.colsize){
	  Rf_error("Matrix is not square in SCYTHE::diag()");
	}
	Matrix newdata(A.rowsize, 1);
	for (int i=0; i<A.rowsize; ++i){
	    newdata.data[i] = A.data[i*A.colsize + i];
	}
	return newdata;
    }
    
    
//! Fast calculation of \a A * \a B + \a C
/*!
 * Fast calculation of \a A * \a B + \a C
 * \param A a constant reference to a Matrix \a A.
 * \return the final Matrix (the value of \a A * \a B + \a C).
 */
    Matrix gaxpy(const Matrix& A, const Matrix& B, const Matrix& C){
	// Case 1: A is 1 x 1 and B is n x k
	if (A.rowsize == 1 && A.colsize == 1) {
	    if (B.rowsize == C.rowsize && B.colsize == C.colsize){
	      Matrix prod(B.rowsize, B.colsize);
		for (int i = 0; i < B.size; ++i) {
		    prod.data[i] = A.data[0] * B.data[i] + C.data[i];
		}
		return prod;
	    } else {
	      Rf_error("A*B and C not conformable in SCYTHE::gaxpy()");
	    }
	} else if (B.rowsize == 1 && B.colsize == 1) {
	    // Case 2: A is n x k and B is 1 x 1
	    if (A.rowsize == C.rowsize && A.colsize == C.colsize){
	      Matrix prod(A.rowsize, A.colsize);
		for (int i = 0; i < A.size; ++i) {
		    prod.data[i] = A.data[i] * B.data[0] + C.data[i];
		}
		return prod;
	    } else {
	      Rf_error("A*B and C not conformable in SCYTHE::gaxpy()");
	    }
	} else if (A.colsize != B.rowsize) {
	    // Case 3: A is n x k and B is m x j (m !=j)
	  Rf_error("Matrices not conformable for multiplication in SCYTHE::gaxpy()");
	} else if (A.rowsize == C.rowsize && B.colsize == C.colsize){
	    // Case 4: A is n x k and B is k x j
	  Matrix newdata(A.rowsize, B.colsize);
	    for (int i = 0; i < A.rowsize; ++i){
		for (int j = 0; j < B.colsize; ++j){
		    newdata.data[i*B.colsize + j] = C.data[i*B.colsize + j];
		    for (int k = 0; k < B.rowsize; ++k){
			newdata.data[i * B.colsize + j] += A.data[i * A.colsize + k] *
			    B.data[k * B.colsize + j];
		    }
		}
	    }
	    
	    return newdata;
	} else {
	  Rf_error("A*B and C not conformable in SCYTHE::gaxpy()");
	}
    }
    
    
//! Fast calculation of \a A'A
/*!
 * Fast calculation of \a A'A
 * \param A a constant reference to a Matrix \a A.
 * \return the final Matrix (the value of \a A'A). 
 */
// original
    Matrix crossprod(const Matrix& A){
      Matrix newdata(A.colsize, A.colsize);

	for (int i = 0; i < A.colsize; ++i){
	    for (int j = i; j < A.colsize; ++j){
		newdata.data[i*A.colsize + j] = 0.0;
		for (int k = 0; k < A.rowsize; ++k){
		    newdata.data[i * A.colsize + j] += A.data[k * A.colsize + i] *
			A.data[k * A.colsize + j];
		    newdata.data[j*A.colsize + i] = newdata.data[i*A.colsize + j];
		}
	    }
	}
	
	return newdata;
    }
    
    // better loop ordering
    Matrix crossprod2(const Matrix& A){
      Matrix newdata(A.colsize, A.colsize);
	const int nr = A.rowsize;
	const int nc = A.colsize;
	
	for (int k = 0; k < nr; ++k){
	    for (int i = 0; i < nc; ++i){
		for (int j = i; j < nc; ++j){
		    newdata.data[i * nc + j] += A.data[k * nc + i] *
			A.data[k * nc + j];
		    newdata.data[j*nc + i] = newdata[i*nc + j];
		}
	    }
	}
	
	/*
	  for (int i = 0; i < A.colsize; ++i){
	  for (int j = i; j < A.colsize; ++j){
	  newdata[i*A.colsize + j] = 0.0;
	  for (int k = 0; k < A.rowsize; ++k){
	  newdata[i * A.colsize + j] += A.data[k * A.colsize + i] *
	  A.data[k * A.colsize + j];
	  newdata[j*A.colsize + i] = newdata[i*A.colsize + j];
	  }
	  }
	  }
	*/
	
	return newdata;

    }
    
    
/************************    OPERATORS     ****************************/
    
// ADDITION OPERATORS IN ALL THEIR MANY FLAVORS
    
/*!
 * \brief OPERATOR: Addition  (Matrix + Matrix).
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the sum of the two Matrices.
 */
    Matrix operator + (const Matrix & A, const Matrix & B)
    {
	if (A.rowsize == 1 && A.colsize == 1) {
	    // Case 1: A is 1 x 1 and B is n x k
	  Matrix sum(B.rowsize, B.colsize);
	    for (int i = 0; i < B.size; ++i) {
		sum.data[i] = A.data[0] + B.data[i];
	    }
	    return sum;
	} else if (B.rowsize == 1 && B.colsize == 1) {
	    // Case 2: A is n x k and B is 1 x 1
	  Matrix sum(A.rowsize, A.colsize);
	    for (int i = 0; i < A.size; ++i) {
		sum.data[i] = A.data[i] + B.data[0];
	    }
	    return sum;
	} else if (A.rowsize != B.rowsize || A.colsize != B.colsize) {
	    // Case 3: A is n x k and B is m x j (n != m or k != m)
	  Rf_error("Matrices not conformable for addition");
	} else {
	    // Case 4: A is n x k and B is also n x k
  	  Matrix sum(A.rowsize, A.colsize);
	    for (int i = 0; i < A.size; ++i) {
		sum.data[i] = A.data[i] + B.data[i];
	    }
	    return sum;
	}
    }
    
//  OPERATOR: Addition
//  NOTE: This operator is overloaded
//  Matrix + scalar
/*!
  \overload Matrix operator + (const Matrix & A, const double &b)
*/
    Matrix operator + (const Matrix & A, const double &b)
    {
      Matrix sum(A.rowsize, A.colsize);
	for (int i = 0; i < A.size; ++i) {
	    sum.data[i] = A.data[i] + b;
	}
	return sum;
    }
    
//  OPERATOR: Addition
//  NOTE: This operator is overloaded
//  scalar + Matrix
/*!
  \overload Matrix operator + (const double &a, const Matrix & B)
*/
    Matrix operator + (const double &a, const Matrix & B)
    {
      Matrix sum(B.rowsize, B.colsize);
	for (int i = 0; i < B.size; ++i) {
	    sum.data[i] = a + B.data[i];
	}
	return sum;
    }
    
// SUBTRACTION OPERATORS IN ALL THEIR MANY FLAVORS
    
/*!
 * \brief OPERATOR: Subtraction  (Matrix - Matrix).
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the difference of the two Matrices.
 */
    Matrix operator - (const Matrix & A, const Matrix & B)
    {
	// Case 1: A is 1 x 1 and B is n x k
	if (A.rowsize == 1 && A.colsize == 1) {
	  Matrix sum(B.rowsize, B.colsize);
	    for (int i = 0; i < B.size; ++i) {
		sum.data[i] = A.data[0] - B.data[i];
	    }
	    return sum;
	} else if (B.rowsize == 1 && B.colsize == 1) {
	    // Case 2: A is n x k and B is 1 x 1
	  Matrix sum(A.rowsize, A.colsize);
	    for (int i = 0; i < A.size; ++i) {
		sum.data[i] = A.data[i] - B.data[0];
	    }
	    return sum;
	} else if (A.rowsize != B.rowsize || A.colsize != B.colsize) {
	    // Case 3: A is n x k and B is m x j (n != m or k != m)
	  Rf_error("Matrices not conformable for subtraction");
	} else {
	    // Case 4: A is n x k and B is also n x k
 	  Matrix sum(A.rowsize, A.colsize);
	    for (int i = 0; i < A.size; ++i) {
		sum.data[i] = A.data[i] - B.data[i];
	    }
	    return sum;
	}
    }
    
//  OPERATOR: Subtraction
//  NOTE: This operator is overloaded
//  Matrix - scalar
/*!
  \overload Matrix operator - (const Matrix & A, const double &b)
*/
    Matrix operator - (const Matrix & A, const double &b)
    {
        Matrix sum(A.rowsize, A.colsize);
	for (int i = 0; i < A.size; ++i) {
	    sum.data[i] = A.data[i] - b;
	}
	return sum;
    }
    
//  OPERATOR: Subtraction
//  NOTE: This operator is overloaded
//  scalar - Matrix
/*!
  \overload Matrix operator - (const double &a, const Matrix & B)
*/
    Matrix operator - (const double &a, const Matrix & B)
    {
      Matrix sum(B.rowsize, B.colsize);
	for (int i = 0; i < B.size; ++i) {
	    sum.data[i] = a - B.data[i];
	}
	return sum;
    }
    
// MULTIPLICATION OPERATORS IN ALL THEIR MANY FLAVORS
    
/*!
 * \brief OPERATOR: Multiplication  (Matrix * Matrix).
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the product of the two Matrices.
 */
    Matrix operator *(const Matrix & A, const Matrix & B)
    {
	// Case 1: A is 1 x 1 and B is n x k
	if (A.rowsize == 1 && A.colsize == 1) {
	  Matrix prod(B.rowsize, B.colsize);
	    for (int i = 0; i < B.size; ++i) {
		prod.data[i] = A.data[0] * B.data[i];
	    }
	    return prod;
	} else if (B.rowsize == 1 && B.colsize == 1) {
	    // Case 2: A is n x k and B is 1 x 1
	  Matrix prod(A.rowsize, A.colsize);
	    for (int i = 0; i < A.size; ++i) {
		prod.data[i] = A.data[i] * B.data[0];
	    }
	    return prod;
	} else if (A.colsize != B.rowsize) {
	    // Case 3: A is n x k and B is m x j (m !=j)
	  Rf_error("Matrices not conformable for multiplication");
	} else {
	    // Case 4: A is n x k and B is k x j
	  Matrix newdata(A.rowsize, B.colsize);
	    for (int i = 0; i < A.rowsize; ++i) {
		for (int j = 0; j < B.colsize; ++j) {
		    newdata.data[i*B.colsize + j] = 0.0;
		    for (int k = 0; k < B.rowsize; ++k) {
			newdata.data[i * B.colsize + j] += A.data[i * A.colsize + k] *
			    B.data[k * B.colsize + j];
		    }
		}
	    }
	    
	    return newdata;
	}
    }
    
//  OPERATOR: Multiplication
//  NOTE: This operator is overloaded
// Matrix * scalar
/*!
  \overload Matrix operator *(const Matrix & A, const double &b)
*/
    Matrix operator *(const Matrix & A, const double &b)
    {
      Matrix prod(A.rowsize, A.colsize);
	for (int i = 0; i < A.size; ++i) {
	    prod.data[i] = A.data[i] * b;
	}
	return prod;
    }
    
//  OPERATOR: Multiplication
//  NOTE: This operator is overloaded
// scalar * Matrix
/*!
  \overload Matrix operator *(const double &a, const Matrix & B)
*/
    Matrix operator *(const double &a, const Matrix & B)
    {
      Matrix prod(B.rowsize, B.colsize);
	for (int i = 0; i < B.size; ++i) {
	    prod.data[i] = a * B.data[i];
	}
	return prod;
    }
    
// ELEMENT BY ELEMENT DIVISION
/*!
 * \brief OPERATOR: Division  (Matrix / Matrix).
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the dividend of each element of the two Matrices.
 */
    Matrix operator / (const Matrix & A, const Matrix & B)
    {
	// Case 1: A is 1 x 1 and B is n x k
	if (A.rowsize == 1 && A.colsize == 1) {
 	  Matrix quot(B.rowsize, B.colsize);
	    for (int i = 0; i < B.size; ++i) {
		quot.data[i] = A.data[0] / B.data[i];
	    }
	    return quot;
	} else if (B.rowsize == 1 && B.colsize == 1) {
	    // Case 2: A is n x k and B is 1 x 1
  	  Matrix quot(A.rowsize, A.colsize);
	    for (int i = 0; i < A.size; ++i) {
		quot.data[i] = A.data[i] / B.data[0];
	    }
	    return quot;
	} else if (A.rowsize != B.rowsize || A.colsize != B.colsize) {
	    // Case 3: A is n x k and B is m x j (n != m or k != m)
	  Rf_error("Matrices not conformable for division");
	} else {
	    // Case 4: A is n x k and B is also n x k
	  Matrix quot(A.rowsize, A.colsize);
	    for (int i = 0; i < A.size; ++i) {
		quot.data[i] = A.data[i] / B.data[i];
	    }
	    return quot;
	}
    }
    
//  OPERATOR: Division
//  NOTE: This operator is overloaded
//  Matrix / scalar
/*!
  \overload Matrix operator / (const Matrix & A, const double &b)
*/
    Matrix operator / (const Matrix & A, const double &b)
    {
	
      Matrix quot(A.rowsize, A.colsize);
	for (int i = 0; i < A.size; ++i) {
	    quot.data[i] = A.data[i] / b;
	}
	return quot;
    }
    
//  OPERATOR: Division
//  NOTE: This operator is overloaded
//  scalar / Matrix
/*!
  \overload Matrix operator / (const double &a, const Matrix & B)
*/
    Matrix operator / (const double &a, const Matrix & B)
    {
	
      Matrix quot(B.rowsize, B.colsize);
	for (int i = 0; i < B.size; ++i) {
	    quot.data[i] = a / B.data[i];
	}
	return quot;
    }
    
    
//!  OPERATOR: Kronecker multiplication
/*!
 * OPERATOR: Kronecker Multiplication  (Matrix % Matrix).
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the result of the Kronecker Multiplication of the 2 Matrices.
 */
    Matrix operator % (const Matrix& A, const Matrix& B){
      Matrix newdata(A.rowsize*B.rowsize, A.colsize*B.colsize);
	int count = 0;
	for (int i=0; i<A.rowsize; ++i){
	    for (int j=0; j<B.rowsize; ++j){
		for (int k=0; k<A.colsize; ++k){
		    for (int m=0; m<B.colsize; ++m){
			newdata.data[count] = A.data[i*A.colsize + k] * B.data[j*B.colsize + m];
			++count;
		    }
		}
	    }
	}
	return newdata;
    }
    
    
//!  OPERATOR: Equality
/*!
 * OPERATOR: Equality.
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return an integer, 1 if the Matrices are the same, 0 otherwise.
 */
    int operator == (const Matrix& A, const Matrix & B){
	if (A.rowsize != B.rowsize || A.colsize != B.colsize) return 0;
	for(int i=0; i<A.size; ++i){
	    if (A.data[i] != B.data[i]) return 0;
	}
	return 1;
    }
    
//  OPERATOR: Inequality
/*!
 * OPERATOR: Inequality
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return an integer, 1 if the Matrices are different, 0 if they are identical.
 */
    int operator != (const Matrix& A, const Matrix & B){
	return !(A==B);
    }
    
//!  OPERATOR: Element-by-element Greater Than
/*!
 * OPERATOR: Element-by-element Greater Than.
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return A Matrix of 1's and 0's, 1's if the element is 
 * greater than the other element, 0 otherwise.
 */
    Matrix operator >> (const Matrix& A, const Matrix& B){
	if (A.rowsize != B.rowsize && A.colsize != B.colsize && (B.size > 1)){
	  Rf_error("Matrices not conformable for >> operator");
	}
	
	if (A.rowsize == B.rowsize && A.colsize == B.colsize){
	  Matrix newdata(A.rowsize, A.colsize);
	    for (int i = 0; i<A.size; ++i){
		newdata.data[i] = A.data[i] > B.data[i];
	    }
	    return newdata;
	}
	
	if (A.rowsize == B.rowsize && B.colsize == 1){
    	    Matrix newdata(A.rowsize, A.colsize);
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=0; j<A.colsize; ++j){
		    newdata.data[i*A.colsize + j] = A.data[i*A.colsize + j] > B.data[i];
		}
	    }
	    return newdata;
	}
	
	if (A.colsize == B.colsize && B.rowsize == 1){
    	    Matrix newdata(A.rowsize, A.colsize);
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=0; j<A.colsize; ++j){
		    newdata.data[i*A.colsize + j] = A.data[i*A.colsize + j] > B.data[j];
		}
	    }
	    return newdata;
	}
	
	if (B.size == 1){
   	    Matrix newdata(A.rowsize, A.colsize);
	    for (int i=0; i<A.size; ++i){
		newdata.data[i] = A.data[i] > B.data[0];
	    }
	    return newdata;  
	} else {
	  Rf_error("Matrices not conformable for >> operator");
	}  
    }
    
//  OPERATOR: Element-by-element Greater Than
/*!
  \overload Matrix operator >> (const Matrix& A, const double& b)
*/
    Matrix operator >> (const Matrix& A, const double& b){
	Matrix newdata(A.rowsize, A.colsize);
	for (int i=0; i<A.size; ++i){
	    newdata.data[i] = A.data[i] > b;
	}
	return newdata;
    }
    
//!  OPERATOR: Element-by-scalar Less Than
/*!
 * OPERATOR: Element-by-element Less Than.
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return A Matrix of 1's and 0's, 1's if the element 
 * is less than the other element, 0 otherwise.
 */
    Matrix operator << (const Matrix& A, const Matrix& B){
	if (A.rowsize != B.rowsize && A.colsize != B.colsize && (B.size > 1)){
	  Rf_error("Matrices not conformable for << operator");
	}
	
	if (A.rowsize == B.rowsize && A.colsize == B.colsize){
  	  Matrix newdata(A.rowsize, A.colsize);
	    for (int i = 0; i<A.size; ++i){
		newdata.data[i] = A.data[i] < B.data[i];
	    }
	    return newdata;
	}
	
	if (A.rowsize == B.rowsize && B.colsize == 1){
	  Matrix newdata(A.rowsize, A.colsize);
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=0; j<A.colsize; ++j){
		    newdata.data[i*A.colsize + j] = A.data[i*A.colsize + j] < B.data[i];
		}
	    }
	    return newdata;
	}
	
	if (A.colsize == B.colsize && B.rowsize == 1){
	  Matrix newdata(A.rowsize, A.colsize);
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=0; j<A.colsize; ++j){
		    newdata.data[i*A.colsize + j] = A.data[i*A.colsize + j] < B.data[j];
		}
	    }
	    return newdata;
	}
	
	if (B.size == 1){
	  Matrix newdata(A.rowsize, A.colsize);
	    for (int i=0; i<A.size; ++i){
		newdata.data[i] = A.data[i] < B.data[0];
	    }
	    return newdata;  
	}
	
	else {
	  Rf_error("Matrices not conformable for << operator");
	}  
    }
    
//  OPERATOR: Element-by-scalar Less Than
/*!
  \overload Matrix operator << (const Matrix& A, const double& b)
*/
    Matrix operator << (const Matrix& A, const double& b)
    {
        Matrix newdata(A.rowsize, A.colsize);
	for (int i=0; i<A.size; ++i){
	    newdata.data[i] = A.data[i] < b;
	}
	return newdata;
    }
    
    
//!  OPERATOR: Element-by-element Equality
/*!
 * OPERATOR: Element-by-element Equality.
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return A Matrix of 1's and 0's, 1's if the elements are equal, 
 * 0 otherwise.
 */
    Matrix operator ^= (const Matrix& A, const Matrix& B){
	if (A.rowsize != B.rowsize && A.colsize != B.colsize && (B.size > 1)){
	  Rf_error("Matrices not conformable for ^= operator");
	}
	
	if (A.rowsize == B.rowsize && A.colsize == B.colsize){
	  Matrix newdata(A.rowsize, A.colsize);
	    for (int i = 0; i<A.size; ++i){
		newdata.data[i] = A.data[i] == B.data[i];
	    }
	    return newdata;
	}
	
	if (A.rowsize == B.rowsize && B.colsize == 1){
	  Matrix newdata(A.rowsize, A.colsize);
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=0; j<A.colsize; ++j){
		    newdata.data[i*A.colsize + j] = A.data[i*A.colsize + j] == B.data[i];
		}
	    }
	    return newdata;
	}
	
	if (A.colsize == B.colsize && B.rowsize == 1){
	  Matrix newdata(A.rowsize, A.colsize);
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=0; j<A.colsize; ++j){
		    newdata.data[i*A.colsize + j] = A.data[i*A.colsize + j] == B.data[j];
		}
	    }
	    return newdata;
	}
	
	if (B.size == 1){
	  Matrix newdata(A.rowsize, A.colsize);
	    for (int i=0; i<A.size; ++i){
		newdata.data[i] = A.data[i] == B.data[0];
	    }
	    return newdata;  
	}
	
	else {
	  Rf_error("Matrices not conformable for ^= operator");
	}  
    }
    
    
//  OPERATOR: Element-by-element Equality
/*!
  \overload Matrix operator ^= (const Matrix& A, const double& b)
*/
    Matrix operator ^= (const Matrix& A, const double& b)
    {
      Matrix newdata(A.rowsize, A.colsize);
	for (int i=0; i<A.size; ++i){
	    newdata.data[i] = A.data[i] == b;
	}
	return newdata;
    }
    
    
    
    
} // namespace dec

#endif /* SCYTHE_DOUBLE_MATRIX_CC */
