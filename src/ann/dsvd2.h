#include "config.h"

/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

/*
 * File: dsvd2.h
 * -----------------
 * This file provides the header file for dsvd2.c
 */

#ifndef _dsvd2_h_
#define _dsvd2_h_

/*
 * Function: dsvd2 
 * Usage: void dsvd2(a,m,n,u,w,v)
 * ------------------------------
 * Its purpose is to perform an svd on input data.  It does this by
 * accepting the data, performing the svd in dsvdcmp1(), reordering
 * the svd in reorder_dsvd().  The function reorder_dsvd() reorders
 * the singular values in descending order, and reorders the
 * columns of the u and v matrices accordingly.
 *
 * To use this function one must link dsvdcmp1.o, reorder_dsvd.o, dindexx.o,
 * and nrutil.o.  The latter file is required because of the vector
 * generation in the first three functions.  The file dindexx.o performs
 * the reordering of the singular values.
 *
 * Inputs:
 * **a The matrix whose svd is to be computed. (double)
 * m   The number of rows in a. (int)
 * n   The number of columns in a. (int)
 *
 * Outputs:
 * **u The coeficient matrix of the svd. (double)
 * *w  The vector of singular values. (double)
 * **v The matrix of singular vectors. (double)
 */
void dsvd2(double **a, int m, int n, double **u,double *w,double **v);

/*
 * Function: dsvdcmp1
 * Usage: void dsvdcmp1(a,m,n,w,v)
 * -------------------------------
 * This file will perform the acutal svd, when called from
 * dsvd2.c, you can use reorder_dsvd.c if you want after
 * this function is called to order the singular values
 * in descending order.
 */
void dsvdcmp1(double **a, int m,int n,double *w,double **v);

/*
 * Function: reorder_dsvd
 * Usage: void reorder_dsvd(u,w,v,m,n)
 * -----------------------------------
 * Its purpose is to reorder the terms in an svd so that the values
 * in the singular values vector are in decreasing order, and the
 * u and v matrices are rearranged accordingly.
 */
void reorder_dsvd(double **u,double *w,double **v,int m,int n);

/*
 * Function: dindexx
 * Usage: void dindexx(n,arrin,indx)
 * ---------------------------------
 * This function will index an array arrin[1...n],
 * i.e. outputs the array indx[1...n] such that arrin[indx[j]]
 * is in ascending order for j=1,2,...,n.  The
 * input quantities,n and arrin are not changed.
 */
void dindexx(int n,double arrin[],int indx[]);

#endif


