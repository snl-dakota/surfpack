/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack_config.h"

/*
 * File: convert.h
 * ---------------
 * This file will convert data from std::vector< std::vector< double > > to a **matrix
 * and back.  It will also convert a DakotaRealArray to a *array and
 * back.
 */


#ifndef _convert_h
#define _convert_h

//#include "data_types.h"

/*
 * Function: DakotaToPtrMatrix
 * Usage: DakotaToPtrMatrix(matrix,drm)
 * ------------------------------------
 * This function will convert data from a std::vector< std::vector< double > >
 * to a **matrix.
 */
void DakotaToPtrMatrix(std::vector< std::vector< double > >& drm,double*** matrix);

/*
 * Function: PtrMatrixToDakota
 * Usage: PtrMatrixToDakota(drm,matrix)
 * ------------------------------------
 * This function will convert data from a **matrix to a
 * std::vector< std::vector< double > >.
 */
void PtrMatrixToDakota(double ***matrix, std::vector< std::vector< double > >& drm);

/*
 * Function: DakotaToPtrArray
 * Usage: DakotaToPtrArray(array,dra)
 * ----------------------------------
 * This function will convert data from a DakotaRealArray
 * to a *array.
 */
void DakotaToPtrArray(std::vector<double>& dra,double** array);

/*
 * Function: PtrArrayToDakota
 * Usage: PtrArrayToDakota(dra,array)
 * ----------------------------------
 * This function will convert data from a **array to a
 * DakotaRealArray.
 */
void PtrArrayToDakota(double **array, std::vector<double>& dra);


#endif
