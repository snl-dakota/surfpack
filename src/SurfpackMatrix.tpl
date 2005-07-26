// See SurfpackMatrix.h for class and method details

template< typename T >
SurfpackVector< T >::SurfpackVector(T* ptr_to_first, unsigned n_elements)
  : ptrToFirst(ptr_to_first), nElements(n_elements), 
  nRows(0), forFortran(false)
{

}

template< typename T >
SurfpackVector< T >::SurfpackVector(T* ptr_to_first, unsigned n_elements,
  unsigned n_rows, bool for_fortran)
  : ptrToFirst(ptr_to_first), nElements(n_elements), 
  nRows(n_rows), forFortran(for_fortran)
{

}

template< typename T >
SurfpackVector< T >::SurfpackVector()
  : ptrToFirst(0), nElements(0), nRows(0), forFortran(false)
{

}

template< typename T >
SurfpackVector< T >::SurfpackVector(const SurfpackVector< T >& other)
  : ptrToFirst(other.ptrToFirst), nElements(other.nElements), 
  nRows(other.nRows), forFortran(other.forFortran)
{

}

template< typename T >
SurfpackVector<T>& SurfpackVector<T>::operator()(T* ptr_to_first)
{
  ptrToFirst = ptr_to_first;
  return *this;
}

template< typename T >
const std::string SurfpackVector< T >::asString()
{
  std::ostringstream os;
  unsigned index;
  for (unsigned i =0; i < nElements; i++) {
    // In C, members of the row are contiguous.  In Fortran, if the matrix has
    // n rows, the elements of any row are spaced n elements apart.
    index = (forFortran) ? i*nRows : i;
    os << ptrToFirst[index] << " ";
  }
  return os.str();
}

template< typename T >
T& SurfpackVector< T >::operator[](unsigned index)
{
  // In C, members of the row are contiguous.  In Fortran, if the matrix has
  // n rows, the elements of any row are spaced n elements apart.
  return forFortran ? ptrToFirst[nRows*index] : ptrToFirst[index];
}

template< typename T >
SurfpackVector< T >& SurfpackVector<T>::operator=(const SurfpackVector<T>& other)
{
  ///\todo Check for assignment to self
  ptrToFirst = other.ptrToFirst;
  nElements = other.nElements;
  nRows = other.nRows;
  forFortran = other.forFortran;
  return *this;
}

template< typename T >
SurfpackMatrix< T >::SurfpackMatrix(unsigned n_rows, unsigned n_cols,
  bool for_fortran)
  : forFortran(for_fortran), nRows(n_rows), nCols(n_cols), 
  oneRow(0,n_cols,n_rows, for_fortran)
{
  rawData.resize(nRows*nCols);
}

template< typename T >
SurfpackMatrix< T >::SurfpackMatrix(bool for_fortran)
  : forFortran(for_fortran), nRows(1), nCols(1), 
  oneRow(0,nCols,nRows,for_fortran)
{
  rawData.resize(1);
}

template< typename T >
SurfpackMatrix< T>::SurfpackMatrix(const SurfpackMatrix<T>& other)
  : forFortran(other.forFortran), nRows(other.nRows), nCols(other.nCols),
  rawData(other.rawData), oneRow(other.oneRow)
{

}


template< typename T >
void SurfpackMatrix< T >::reshape(unsigned n_rows, unsigned n_cols)
{
  nRows = n_rows;
  nCols = n_cols;
  rawData.resize(n_rows*n_cols);
  oneRow = SurfpackVector< T >(0,n_cols,n_rows,forFortran);
}

template< typename T >
SurfpackVector< T >& SurfpackMatrix< T >::operator[](unsigned row_index)
{
  // In Fortran, element i is the first element of the i-th row.
  // But in C, element i is the 1st-row element at column i, and the first-
  // column element of row i is actually at i*nCols. 
  return forFortran ?
    oneRow(&rawData[row_index]) :
    oneRow(&rawData[row_index*nCols]);
}

template< typename T >
SurfpackMatrix< T >& SurfpackMatrix< T >::operator=(const SurfpackMatrix<T>& other)
{
 ///\todo check for assignment to self
 forFortran = other.forFortran;
 nRows = other.nRows;
 nCols = other.nCols;
 rawData = other.rawData;
 oneRow = other.oneRow;
 return *this;
}

template< typename T >
const std::string SurfpackMatrix< T >::asString()
{
  std::ostringstream os;
  unsigned index;
  for (unsigned i = 0; i < nRows; i++) {
    // In C, members of the row are contiguous.  In Fortran, if the matrix has
    // n rows, the elements of any row are spaced n elements apart.
    index = forFortran ? i : i * nCols;
    os << oneRow(&rawData[index]).asString() << std::endl;
  }
  return os.str();
}

template< typename T >
const std::string SurfpackMatrix< T >::asArrayString()
{
  std::ostringstream os;
  for (unsigned i = 0; i < rawData.size(); i++) {
    os << rawData[i] << " "; 
  }
  return os.str();
}

template< typename T >
unsigned SurfpackMatrix< T >::getNRows()
{
  return nRows;
}

template< typename T >
unsigned SurfpackMatrix< T >::getNCols()
{
  return nCols;
}
