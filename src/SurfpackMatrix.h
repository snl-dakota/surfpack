/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */
#ifndef SURFPACK_MATRIX
#define SURFPACK_MATRIX

#include "surfpack_config.h"
#include "surfpack_system_headers.h"

/// Helper class for SurfpackMatrix.  Allows for the double subscripting of
/// a contiguous chunk of memory.
template< typename T >
class SurfpackVector
{
public:
  /// First argument gives address of first element of the row.
  /// Second argument gives number of elements in the row.
  SurfpackVector(T* ptr_to_first, unsigned n_elements);

  /// First argument gives address of first element of the row.
  /// Second argument gives number of elements in the row.
  /// Third and fourth arguments are necessary if object is used in conjunction
  /// with a SurfpackMatrix that is configured for passing to Fortran APIs.
  /// Row
  SurfpackVector(T* ptr_to_first, unsigned n_elements,
    unsigned n_rows, bool for_fortran = true);

  /// Default constructor.  Data members default to 0/false 
  SurfpackVector();

  /// Copy constructor.  Performs deep copy.
  SurfpackVector(const SurfpackVector< T >& other);

  /// Assignment operator.  Performs deep copy.
  SurfpackVector<T>& operator=(const SurfpackVector< T >& other);

  /// Sets the address of the first element of this row
  SurfpackVector<T>& operator()(T* ptr_to_first);

  /// Returns list of row elements, separated by spaces, NOT newline-terminated
  const std::string asString();

  /// Returns the element of the row indexed by the argument
  T& operator[](unsigned vector_index);
private:
  /// First element of the row
  T* ptrToFirst;

  /// Number of elements in the row
  unsigned nElements;

  /// Number of rows in the Matrix for which this vector is a row.
  /// Only needed if Matrix is configured for use with Fortran APIs
  unsigned nRows;

  /// True if the Matrix this row is a part of is configured for use with
  /// Fortran.  If so, the elements of the row are not contiguous, but are
  /// spaced every nRows elements apart in the contiguous block that the
  bool forFortran;
};
  
/// A SurfpackMatrix uses an STL vector-- where all elements reside in a 
/// contiguous block of memory-- to represent a two-dimensional matrix.
/// Overloading of operator[], together with the use of helper class
/// SurfpackVector, allows for the familiar, convenient, intuitive double-
/// subscripting.  A[i][j] is the element in the i-th row and j-th column of A.
/// In Surfpack, the data in such matrices are frequently passed through to
/// Fortran APIs (for use with BLAS and LAPACK routines).  Fortran expects the
/// data to be organized with elements adjacent to other elements in the column.
/// In C, on the other hand, doubly-subscripted arrays are organized so that
/// elements of the same row are adjacent to each other in memory.  The client
/// of SurfpackMatrix/Vector must specify at the time of object construction
/// whether the data are to be organized after the manner of C or Fortran, via
/// the constructor paramater "for_fortran."  After that, the client may trust
/// that A[i][j] means row-i and column-j, without worrying about the underlying
/// organization of the data.
template< typename T >
class SurfpackMatrix
{
public:
  /// Creates a n_rows x n_cols matrix
  SurfpackMatrix(unsigned n_rows, unsigned n_cols, bool for_fortran = true);
  
  /// Defaults to a 1x1 matrix for use with Fortran
  SurfpackMatrix(bool for_fortran = true);

  /// Copy constructor.  Makes deep copy.
  SurfpackMatrix(const SurfpackMatrix<T>& other);

  /// Resize the matrix.
  void reshape(unsigned n_rows, unsigned n_cols);

  /// Elements separated by spaces; rows delimited by newlines.
  const std::string asString();

  /// Returns the elements of a matrix as a list, in the order they appear in
  /// memory.
  const std::string asArrayString();

  /// For selecting a matrix row.  Should only be used as the first in a pair
  /// of subscripts-- never alone.
  SurfpackVector<T>& operator[](unsigned row_index);

  /// Assignment operator.  Performs deep copy.
  SurfpackMatrix<T>& operator=(const SurfpackMatrix<T>& other);

  /// Assignment operator.  Performs deep copy.
  SurfpackMatrix<T>& operator=(const std::vector< std::vector<T> >);
  
  /// Sum of the elements of the matrix.  Data type T should support operator+.
  T sum();

  /// Number of rows in the matrix.
  unsigned getNRows();

  /// Number of columns in the matrix.
  unsigned getNCols();
private:
  /// True iff matrix can be safely passed to a Fortran API.
  bool forFortran;

  /// Number of rows in the matrix.
  unsigned nRows;

  /// Number of columns in the matrix.
  unsigned nCols;

  /// Contigous block of memory stores elements of the matrix.
  std::vector< T > rawData;

  /// Row of matrix most recently references by operator[]
  SurfpackVector< T > oneRow;
};

// Definitions must be included in same file as declarations
#include "SurfpackMatrix.tpl"

#endif
