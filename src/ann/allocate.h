/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack_config.h"

/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

/*
 * File: allocate.h
 * -----------------
 * This file provides the interface for
 * a utilities library to run toms Neural Networks.
 */

#ifndef _allocate_h_
#define _allocate_h_

/*
 * Function: nrerror
 * Usage: void nrerror(error_text)
 * -------------------------------
 * This funcition will print error messages.
 */
void nrerror(const char* error_text);

/*
 * Function: vector
 * Usage: float *vector(nl,nh)
 * ---------------------------
 * This function will allocate memory for a vector.
 */
float *_vector(int nl,int nh);

/*
 * Function: i_vector
 * Usage: int *i_vector(nl,nh)
 * --------------------------
 * This function will allocate memory for a i_vector.
 */
int *i_vector(int nl,int nh);

/*
 * Function: d_vector
 * Usage: int *d_vector(nl,nh)
 * --------------------------
 * This function will allocate memory for a d_vector.
 */
double *d_vector(int nl,int nh);

/*
 * Function: matrix
 * Usage: float **matrix(nrl,nrh,ncl,nch)
 * --------------------------------------
 * This function will allocate memory for a matrix.
 */
float **matrix(int nrl,int nrh,int ncl,int nch);

/*
 * Function: d_matrix
 * Usage: double **d_matrix(nrl,nrh,ncl,nch)
 * ----------------------------------------
 * This function will allocate memory for a d_matrix.
 */
double **d_matrix(int nrl,int nrh,int ncl,int nch);

/*
 * Function: i_matrix
 * Usage: int **i_matrix(nrl,nrh,ncl,nch)
 * -------------------------------------
 * This function will allocate memory for a i_matrix.
 */
int **i_matrix(int nrl,int nrh,int ncl,int nch);

/*
 * Function: submatrix
 * Usage: float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
 * ---------------------------------------------------------------
 * This function will allocate memory for a submatrix.
 */
float **submatrix(float **a, int oldrl,int oldrh,int oldcl,int oldch,int newrl,
                  int newcl);

/*
 * Function: free_vector
 * Usage: void free_vector(v,nl,nh)
 * --------------------------------
 * This function will free the memory allocated for a vector.
 */
void free_vector(float *v,int nl,int nh);

/*
 * Function: free_i_vector
 * Usage: void free_i_vector(v,nl,nh)
 * ---------------------------------
 * This function will free the memory allocated for a i_vector.
 */
void free_i_vector(int *v,int nl,int nh);

/*
 * Function: free_d_vector
 * Usage: void free_d_vector(v,nl,nh)
 * ---------------------------------
 * This function will free the memory allocated for a d_vector.
 */
void free_d_vector(double *v,int nl,int nh);

/*
 * Function: free_matrix
 * Usage: void free_matrix(m,nrl,nrh,ncl,nch)
 * ------------------------------------------
 * This function will free the memory allocated for a matrix.
 */
void free_matrix(float **m,int nrl,int nrh,int ncl,int nch);

/*
 * Function: free_d_matrix
 * Usage: void free_d_matrix(m,nrl,nrh,ncl,nch)
 * -------------------------------------------
 * This function will free the memory allocated for a d_matrix.
 */
void free_d_matrix(double **m,int nrl,int nrh,int ncl,int nch);

/*
 * Function: free_i_matrix
 * Usage: void free_i_matrix(m,nrl,nrh,ncl,nch)
 * -------------------------------------------
 * This function will free the memory allocated for a i_matrix.
 */
void free_i_matrix(int **m,int nrl,int nrh,int ncl,int nch);

/*
 * Function: free_submatrix
 * Usage: void free_submatrix(b,nrl,nrh,ncl,nch)
 * ---------------------------------------------
 * This function will free the memory allocated for a submatrix.
 */
void free_submatrix(float **b,int nrl,int nrh,int ncl,int nch);

/*
 * Function: convert_matrix
 * Usage: float **convert_matrix(a,nrl,nrh,ncl,nch)
 * ------------------------------------------------
 * This function will allocated memory for a convert_matrix.
 */
float **convert_matrix(float *a,int nrl,int nrh,int ncl,int nch);

/*
 * Function: free_convert_matrix
 * Usage: void free_convert_matrix(b,nrl,nrh,ncl,nch)
 * --------------------------------------------------
 * This function will free the memory allocated for a convert_matrix.
 */
void free_convert_matrix(float **b,int nrl,int nrh,int ncl,int nch);

/* dtensor3.c */
/* This file is stored in dtensor3.c.
   First function:
   Its purpose is to dynamically allocate a 3D tensor of user specified
   dimensions.  The function returns a pointer to the pointer to
   the pointer to the first element in the tensor.
   Second function:
   Its purpose is to free the memory that was dynamically allocated
        with the function tensor3().
*/
/*
 * Function: dtensor3
 * Usage: double ***dtensor3(n1l, n1h, n2l, n2h, n3l, n3h)
 * ---------------------------------------------------------
 * Its purpose is to dynamically allocate a 3D tensor of user specified
 * dimensions.  The function returns a pointer to the pointer to
 * the pointer to the first element in the tensor.
 */
double ***dtensor3(int n1l, int n1h, int n2l, int n2h, int n3l, int n3h);

/*
 * Function: free_dtensor3
 * Usage: void free_dtensor3(t, n1l, n1h, n2l, n2h, n3l, n3h)
 * ------------------------------------------------------------
 * Its purpose is to free the memory that was dynamically allocated
 * with the function tensor3().
 */
void free_dtensor3(double ***t, int n1l, int n1h, int n2l, int n2h, int n3l,
                   int n3h);

/* lmatrix.c */
/* This file is stored in lmatrix.c.
   The first function can be used to allocate a matrix of long
   integer variables.
   The second function can be used to deallocate the matrix.
*/

/*
 * Function: lmatrix
 * Usage: long int **lmatrix(nrl,nrh,ncl,nch)
 * ------------------------------------------
 * This function will free the memory allocated for a lmatrix of long integers.
 */
long int **lmatrix(int nrl,int nrh,int ncl,int nch);

/*
 * Function: free_lmatrix
 * Usage: void free_lmatrix(m,nrl,nrh,ncl,nch)
 * -------------------------------------------
 * This function will free the memory allocated for a lmatrix.
 */
void free_lmatrix(long int **m,int nrl,int nrh,int ncl,int nch);
#endif

