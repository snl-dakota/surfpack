#include "config.h"

/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

/* This function will compute the arctanh for given
   value x
*/

#include <vector>
#include "vector_enhancements.h"

//#include "data_types.h"
#include "system_defs.h"
#include "utilities.h"
#include <cmath>
using namespace std;
/*
 * Function: my_atanh
 * Usage: double my_atanh(x)
 * -------------------------
 * This function will compute the arctanh(x).
 */
double my_atanh(double x)
{
  if (fabs(x) >= 1.0) {
    cerr << "atanh(x) ERROR: |x| must be < 1.0" << endl;
    exit(-1);
  }

  return( 0.5*log((1.0+x)/(1.0-x)) );
}

/*
 * Function: multiplyMM 
 * Usage: DakotaRealMatix = multiplyMM(vector< vector< double > > m1, vector< vector< double > > m2)
 * -----------------------------------------------------------------------------
 * This function will multiply two matrices (result = m1*m2). 
 */
vector< vector< double > > multiplyMM(vector< vector< double > >& m1, vector< vector< double > >& m2)
{
  int j,k,l;
  double sum;   
  int num_rows_ = num_rows(m1);
  int num_cols_ = num_columns(m2);
  int inner    = num_columns(m1);
  vector< vector< double > > result;
  reshape_2d(result,num_rows_,num_cols_);
   
  if (inner != num_rows(m2)) {
    cerr << "multiplyMM: Inner matrix dimensions must agree." << endl;
    exit(-1);
  }

  for (j=0;j<num_rows_;j++) {
    for (k=0;k<num_cols_;k++) {
      sum = 0.0;
      for (l=0;l<inner;l++)
	sum += m1[j][l] * m2[l][k];
      result[j][k] = sum;
    }
  }

  return(result);
}

/*
 * Function: FindMax 
 * Usage: result = FindMax(vector< vector< double > > m) 
 * -------------------------------------------
 * This function will return the maximum value for each
 * column of the matrix m in the form of a vector. 
 */
vector<double> FindMax(vector< vector< double > > m)
{
  int j,k,num_cols = num_columns(m);
  vector<double> result(num_cols);

  for(j=0;j<num_cols;j++)
    result[j] = m[0][j];

  for(j=0;j<num_rows(m);j++)
    for(k=0;k<num_cols;k++)
      if(m[j][k] > result[k])
	result[k] = m[j][k];

  return(result);
}

/*
 * Function: FindMin
 * Usage: result = FindMin(vector< vector< double > > m)
 * -------------------------------------------
 * This function will return the minimum value for each
 * column of the matrix m in the form of a vector.
 */
vector<double> FindMin(vector< vector< double > > m)
{
  int j,k,num_cols = num_columns(m);
  vector<double> result(num_cols);

  for(j=0;j<num_cols;j++)
    result[j] = m[0][j];

  for(j=0;j<num_rows(m);j++)
    for(k=0;k<num_cols;k++)
      if(m[j][k] < result[k])
	result[k] = m[j][k];

  return(result);
}
