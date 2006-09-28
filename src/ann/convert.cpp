/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

/*
 * File: convert.C
 * ---------------
 * This file will convert data from vector< vector< double > > to a **matrix
 * and back.  It will also convert a DakotaRealArray to a *array and 
 * back.
 */

#include <vector>
#include "vector_enhancements.h"
#include "surfpack_system_headers.h"
#include "convert.h"
#include "allocate.h"

using std::vector;

/*
 * Function: DakotaToPtrMatrix 
 * Usage: DakotaToPtrMatrix(matrix,drm)
 * ------------------------------------
 * This function will convert data from a vector< vector< double > >
 * to a **matrix.
 */
void DakotaToPtrMatrix(vector< vector< double > >& drm,double*** matrix)
{
   int j,k;
   int num_rows_,num_columns_;

   num_rows_ = num_rows(drm);
   num_columns_ = num_columns(drm);

   *matrix = d_matrix(0,num_rows_-1,0,num_columns_-1);
   for(j=0;j<num_rows_;j++) {
      for(k=0;k<num_columns_;k++) {
         (*matrix)[j][k] = drm[j][k];
      }
   }
}


/*
 * Function: PtrMatrixToDakota
 * Usage: PtrMatrixToDakota(drm,matrix)
 * ------------------------------------
 * This function will convert data from a **matrix to a
 * vector< vector< double > >.
 */
void PtrMatrixToDakota(double ***matrix, vector< vector< double > >& drm)
{
   int j,k;
   int num_rows_,num_columns_;

   num_rows_ = num_rows(drm);
   num_columns_ = num_columns(drm);

   for(j=0;j<num_rows_;j++) {
      for(k=0;k<num_columns_;k++) {
         drm[j][k] = (*matrix)[j][k];
      }
   }
   free_d_matrix(*matrix,0,num_rows_-1,0,num_columns_-1);
}
// Convert an Array

/*
 * Function: DakotaToPtrArray
 * Usage: DakotaToPtrArray(array,dra)
 * ----------------------------------
 * This function will convert data from a DakotaRealArray
 * to a *array.
 */
void DakotaToPtrArray(vector<double>& dra,double** array)
{
   int j;
   int num_columns_;

   num_columns_ = dra.size();

   *array = d_vector(0,num_columns_-1);
   for(j=0;j<num_columns_;j++) {
         (*array)[j] = dra[j];
   }
}

/*
 * Function: PtrArrayToDakota
 * Usage: PtrArrayToDakota(dra,array)
 * ----------------------------------
 * This function will convert data from a **array to a
 * DakotaRealArray.
 */
void PtrArrayToDakota(double **array, vector<double>& dra)
{
   int j;
   int num_columns_;

   num_columns_ = dra.size();

   for(j=0;j<num_columns_;j++) {
         dra[j] = (*array)[j];
   }
   free_d_vector(*array,0,num_columns_-1);
}
