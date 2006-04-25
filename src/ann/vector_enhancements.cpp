/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.

    Surfpack: A Software Library of Multidimensional Surface Fitting Methods

    Surfpack is distributed under the DAKOTA GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#include "surfpack_config.h"

#include <vector>
using namespace std;

void reshape_2d(vector< vector< double > >& matrix, int rows, int columns)
{
  matrix.resize(rows);
  for (int i = 0; i < rows; i++) {
    matrix[i].resize(columns);
  } 
}

int num_rows(std::vector< std::vector< double> >& matrix)
{
  return matrix.size();
}

int num_columns(std::vector< std::vector< double> >& matrix)
{
  return (matrix.size() == 0) ? 0 : matrix[0].size();
}

