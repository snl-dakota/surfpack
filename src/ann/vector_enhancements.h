/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.

    Surfpack: A Software Library of Multidimensional Surface Fitting Methods

    Surfpack is distributed under the DAKOTA GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#include "surfpack_config.h"

void reshape_2d(std::vector< std::vector< double > >& matrix, int rows, int columns);
int num_rows(std::vector< std::vector< double> >& matrix);
int num_columns(std::vector< std::vector< double> >& matrix);
