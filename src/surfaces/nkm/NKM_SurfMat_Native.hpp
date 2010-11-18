#ifndef __SURFMAT_HPP__
#define __SURFMAT_HPP__
#include <vector>
//#include <iostream>
#include <string>
#include <sstream>

namespace nkm {

///allows a #define statement to toggle between intializing and not initializing to zero, use this to allocate data BUT NOT data2D (since we ALWAYS need to set that ourselves anyway); 
inline void* data_alloc(size_t nelem, size_t elemsize){
#ifdef __SURFMAT_ZERO_MEM__
  return calloc(nelem,elemsize);
#else
  return malloc(nelem*elemsize);
#endif  
};

/************************************************************************/
/************************************************************************/
/************************************************************************/
/**************** The SurfMat Class template Starts Here ****************/
/************************************************************************/
/************************************************************************/
/************************************************************************/

template<typename T>
class SurfMat
{
public:
  //default constructor, does not allocate memory
  SurfMat() {
    tol  =0;  //of type T, so one type conversion every construction
    NRows=NCols=0;
    data  =NULL;
    data2D=NULL;
    return;
  };
  
  //typical matrix constructor, defaults to 1 column (a vector) if the number of columns is not specified
  SurfMat(int nrows_in, int ncols_in=1);
  
  //copy constructor, make a deep copy
  SurfMat(const SurfMat<T>& other);
  
  //copy from vector constructor
  SurfMat(const std::vector<T> &vecT, int nrows_in=-99999, int ncols_in=1);

  //destructor 
  ~SurfMat() {
    clear();
    return;
  };
  
  //deallocate memory
  void clear();
  
  //change this matrix to a vector with nrows_new rows, it does not attempt to preserve values in memory which, if a change is needed, makes it faster than reshape2 or resize2: this function doesn't do anything if you request a new size the same as the current size
  inline void newSize(int nrows_new){
    if((NRows!=nrows_new)||(NCols!=1))
      newSize2(nrows_new,1);
    return;};
  
  //change this matrix to a matrix with nrows_new rows and ncols_new columns, it does not attempt to preserve values in memory which, if a change is needed, makes it faster than reshape2 or resize2: this function doesn't do anything if you request a new size the same as the current size
  inline void newSize(int nrows_new, int ncols_new){
    if((NRows!=nrows_new)||(NCols!=ncols_new)) 
      newSize2(nrows_new,ncols_new);
    return;};
  
  //enlarge or shrink the matrix/Vector while keeping continguous in memory elements contiguous in memory: this function doesn't do anything if you request a new size the same as the current size
  inline void reshape(int nrows_new){
    if((NRows!=nrows_new)||(NCols!=1))
      reshape2(nrows_new,1);
    return;};
  
  //enlarge or shrink the matrix while keeping contigous in memory elements contiguous in memory: this function doesn't do anything if you request a new size the same as the current size
  inline void reshape(int nrows_new, int ncols_new){
    if((NRows!=nrows_new)||(NCols!=ncols_new)) 
      reshape2(nrows_new,ncols_new);
    return;};
  
  //enlarge or shrink the matrix/Vector while adding zeros after or chopping everything off the end of the first row, this function doesn't do anything if you request a new size the same as the current size
  inline void resize(int nrows_new){
    if((NRows!=nrows_new)||(NCols!=1))
      resize2(nrows_new,1);
    return;};
  
  //enlarge or shrink the matrix while adding zeros after the last row and/or last column and/or chopping off the rows and/or columns after the newly requested last row and/or column: this function doesn't do anything if you request a new size the same as the current size
  inline void resize(int nrows_new, int ncols_new){
    if((NRows!=nrows_new)||(NCols!=ncols_new)) 
      resize2(nrows_new,ncols_new);
    return;
  };
  
  ///set all of the matrix's elements to zero
  //void zero(){int nelems=getNElems(); for(int k=0; k<nelems; k++) data[k]=0; return;};
  inline void zero(){int nelems=getNElems(); for(int k=0; k<nelems; k++) data[k]=0; return;};

  ///set all of the matrix's elements to zero, except for the diagonal (i=j) which is set to 1, works for rectangular matrices
  inline void identity(){
    zero(); 
    int n=(NRows<NCols)?NRows:NCols; 
    for(int i=0; i<n; i++) data2D[i][i]=1;
    return;
  }
  
  ///make a deep copy
  SurfMat<T>& copy(const SurfMat<T>& other);


  // Assignment operator.  Performs deep copy.
  inline SurfMat<T>& operator=(const SurfMat<T>& other){return (copy(other));};
  
  //return the number of rows
  inline int getNRows() const {return NRows;};
  
  //return the number of columns
  inline int getNCols() const {return NCols;};
  
  //return the number of elements
  inline int getNElems() const {return (NRows*NCols);};

  inline const T getTol() const {return tol;};
  inline void putTol(T tol_in){tol=tol_in; return;};
  
  
  //vector style retrieve of a single element from the matrix by value
  inline const T& operator()(int k) const {
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=k)&&(k<(NRows*NCols)));
#endif
    return (const_cast<T&> (data[k]));
  };
  
  //vector style retrieve of a single element from the matrix by reference
  inline T& operator()(int k) {
    //T& operator()(int k) {
#ifdef __SURFMAT_ERR_CHECK__
    if(!((0<=k)&&(k<NRows*NCols))) {
      printf("T& SurfMat::operator()(int k): k=%d NRows=%d NCols=%d",
	     k,NRows,NCols);
      printf("\n");
      assert((0<=k)&&(k<NRows*NCols));
    }
#endif
    return data[k];};
  
  //matrix style retrieve of a single element from the matrix by value
  inline const T& operator()(int i, int j) const {
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=i)&&(i<NRows)&&(0<=j)&&(j<NCols));
#endif
    return const_cast<T&> (data2D[j][i]);};
  
  //matrix style retrieve of a single element from the matrix by reference
  inline T& operator()(int i, int j) {
    //T& operator()(int i, int j) {
#ifdef __SURFMAT_ERR_CHECK__
    if(!((0<=i)&&(i<NRows)&&(0<=j)&&(j<NCols))){
      printf("i=%d NRows=%d j=%d NCols=%d",i,NRows,j,NCols);
      printf("\n"); //so can breakpoint on this line in debugger
      assert((0<=i)&&(i<NRows)&&(0<=j)&&(j<NCols));
    }
#endif
    return data2D[j][i];
  };

  
  //vector style retrieve of pointer to element, for passing to BLAS & LAPACK convenience
  inline T* ptr(int k){
    //T* ptr(int k){
#ifdef __SURFMAT_ERR_CHECK__
    if(!((0<=k)&&(k<NRows*NCols))) {
      printf("ERROR: SurfMat T* ptr(k): need 0<k=%d, k<NRows*NCols=%d, NRows=%d, NCols=%d\n",k,NRows*NCols,NRows,NCols);
      assert((0<=k)&&(k<NRows*NCols));
    }
#endif
    return (data+k);
  };

  inline const T* ptr(int k) const {
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=k)&&(k<NRows*NCols));
#endif
    return (data+k);
  };

  //matrix style retrieve of pointer to element, for passing to BLAS & LAPACK convenience
  inline T* ptr(int i, int j){
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=i)&&(i<NRows)&&(0<=j)&&(j<NCols));
#endif
    return &(data2D[j][i]);
  };

  inline const T* ptr(int i, int j) const {
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=i)&&(i<NRows)&&(0<=j)&&(j<NCols));
#endif
    return &(data2D[j][i]);
  };

  ///returns the smallest value contained in this matrix, wraps T& minElem(T& val ,int& loc)
  inline T minElem() const;
  
  ///returns the location of the smallest value contained in this matrix, if there is a tie it returns the location of the first one encountered, wraps T& minElem(T& val ,int& loc)
  inline int iMinElem() const;
  
  ///returns the smallest value contained in this matrix and its location in loc, if there is a tie it records the location of the first smallest value encountered
  inline void minElem(T& val, int& loc) const;
  
  ///returns the largest value contained in this matrix, wraps T& maxElem(T& val ,int& loc)
  inline T maxElem() const;
  
  ///returns the location of the largest value contained in this matrix, if there is a tie it returns the location of the first one encountered, wraps T& maxElem(T& val ,int& loc)
  inline int iMaxElem() const;

  ///returns the largest value contained in this matrix and its location in loc, if there is a tie it records the location of the first largest value encountered
  inline void maxElem(T& val, int& loc) const;
  
  ///returns the smallest and largest values contained in this matrix
  inline void minMaxElem(T& smallest, T& largest) const {
    int nelem=NRows*NCols;
#ifdef __SURFMAT_ERR_CHECK__
    assert(nelem>0);
#endif
    smallest=largest=data[0];
    for(int k=1; k<nelem; k++){
      //note I did not use an else if for the largest so that both if statements could be evaluated simultaneous by different processing units
      if(data[k]<smallest) smallest=data[k];
      if(data[k]>largest ) largest =data[k];
    }    
    return;
  };
  
  ///replace the values in row "irow" of this matrix with the values stored in vecT_row; this function will NOT expand the size of this matrix if you specify a row index larger than NRows-1.
  inline void putRows(const std::vector<T>& vecT_row, int irow) {
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=irow)&&(irow<NRows)&&(vecT_row.size()==NCols));
#endif
    for(int j=0; j<NCols; ++j) data2D[j][irow]=vecT_row[j];
    return;    
  };

  ///replace the values in row "irow" of this matrix with the entries listed in string s; it returns 0 if the string and row have the same number of entries and 1 otherwise
  inline int putRows(const std::string& s, int irow) {
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=irow)&&(irow<NRows));
#endif
    std::istringstream is(s);
    int j;
    //T temp;
    for(j=0; j<NCols; ++j) {
      if(is.eof())
	break;
      //is >> temp; data2D[j][irow]=temp;
      //printf("irow=%d NRows=%d j=%d NCols=%d temp=%g\n",irow,NRows,j,NCols,temp);
      //printf("data2D[j=%d][i=%d]=%g\n",j,irow,data2D[j][irow]);
      is >> data2D[j][irow];
    }
    return(!((j==NCols)&&(is.eof())));
  };    

  ///replace the values in row "irow" of this matrix with the values stored in "row"; this function will NOT expand the size of this matrix if you specify a row index larger than NRows-1.
  inline void putRows(const SurfMat<T>& row, int irow) {
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=irow)&&(irow<NRows)&&
	   (row.getNRows()==1)&&(row.getNCols()==NCols));
#endif
    for(int j=0; j<NCols; j++) data2D[j][irow]=row.data[j];
    return;    
  };
  
  ///replace the values in multiple rows (which rows they are is specified in irows) of this matrix with the values stored in "rows"; this function will NOT expand the size of this matrix if you specify a row index larger than NRows-1.
  inline void putRows(const SurfMat<T>& rows, SurfMat<int> irows) {
    int nrows_put=irows.getNElems();
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=irows.minElem())&&(irows.maxElem()<NRows)&&
	   (rows.getNRows()==nrows_put)&&(rows.getNCols()==NCols));
#endif
    for(int j=0; j<NCols; j++)
      for(int k=0; k<nrows_put; k++)
	data2D[j][irows(k)]=rows.data2D[j][k];
    return;
  };
  
  ///replace the values in column "jcol" of this matrix with the values stored in vector "vecT_col"; this function will NOT expand the size of this matrix if you specify a column index larger than NCols-1.
  inline void putCols(const std::vector<T>& vecT_col, int jcol) {
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=jcol)&&(jcol<NCols)&&(vecT_col.size()==NRows));
#endif
    for(int i=0; i<NRows; i++)
      data2D[jcol][i]=vecT_col[i];
    return;
  };

  ///replace the values in column "jcol" of this matrix with the entries listed in string s; it returns 0 if the string and column have the same number of entries and 1 otherwise this function will NOT expand the size of this matrix if you specify a column index larger than NCols-1.
  inline void putCols(const std::string& s, int jcol) {
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=jcol)&&(jcol<NCols));
#endif
    std::istringstream is(s);
    int i;
    for(i=0; i<NRows; ++i) {
      if(is.eof())
	break;
      is >> data2D[jcol][i];
    }    
    return(!((i==NRows)&&(is.eof()))); 
  };

  ///replace the values in column "jcol" of this matrix with the values stored in "col"; this function will NOT expand the size of this matrix if you specify a column index larger than NCols-1.
  inline void putCols(const SurfMat<T>& col, int jcol) {
#ifdef __SURFMAT_ERR_CHECK__
    if(!((0<=jcol)&&(jcol<NCols)&&
	 (col.getNRows()==NRows)&&(col.getNCols()==1))){
      //assert((0<=jcol)&&(jcol<NCols)&&
      //(col.getNRows()==NRows)&&(col.getNCols()==1));
      assert(0<=jcol);
      assert(jcol<NCols);
      assert(col.getNRows()==NRows);
      assert(col.getNCols()==1);
    }
#endif
    for(int i=0; i<NRows; i++)
      data2D[jcol][i]=col.data[i];
    return;
  };

  ///replace the values in multiple columns (which columns they are is specified in jcols) of this matrix with the values stored in "cols"; this function will NOT expand the size of this matrix if you specify a column index larger than NCols-1.
  inline void putCols(const SurfMat<T>& cols, SurfMat<int> jcols) {
    int ncols_put=jcols.getNElems();
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=jcols.minElem())&&(jcols.maxElem()<NCols)&&
	   (cols.getNRows()==NRows)&&(cols.getNCols()==ncols_put));
#endif
    for(int k=0; k<ncols_put; k++)
      for(int i=0; i<NRows; i++)
	data2D[jcols(k)][i]=cols.data2D[k][i];
    return;
  };
  
  ///get one row (index stored in irow) of this matrix and return as a "row vector" 
  inline SurfMat<T>& getRows(SurfMat<T>& result, int irow) const {
#ifdef __SURFMAT_ERR_CHECK__
    if(!((0<=irow)&&(irow<NRows))) {
      printf("irow=%d NRows=%d\n",irow,NRows); fflush(stdout);
      assert((0<=irow)&&(irow<NRows));
    }
#endif
    result.newSize(1,NCols);
    result.tol=tol;
    for(int j=0; j<NCols; j++) result(j)=data2D[j][irow];
    return result;
  };
  
  ///get multiple rows (indices stored in matrix irows) of this matrix and return as a new matrix
  inline SurfMat<T>& getRows(SurfMat<T>& result, SurfMat<int>& irows) const {
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=irows.minElem())&&(irows.maxElem()<NRows));
#endif
    int nrows_res=irows.getNElems();
    result.newSize(nrows_res,NCols);
    result.tol=tol;
    if(nrows_res>0) 
      for(int j=0; j<NCols; j++)
	for(int i=0; i<nrows_res; i++)
	  result(i,j)=data2D[j][irows(i)];	
    
    return result;
  };
  
  ///get one column (index stored in jcol) of this matrix and return as a vector
  inline SurfMat<T>& getCols(SurfMat<T>& result, int jcol) const {
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=jcol)&&(jcol<NCols));
#endif
    result.newSize(NRows);
    result.tol=tol;
    for (int i=0; i<NRows; i++) result(i)=data2D[jcol][i];
    return result;
  };
  
  ///get multiple columns (indices stored in matrix jcols) of this matrix and return as a new matrix
  inline SurfMat<T>& getCols(SurfMat<T>& result, SurfMat<int>& jcols) const {
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=jcols.minElem())&&(jcols.maxElem()<NCols));
#endif
    int ncols_res=jcols.getNElems();
    result.newSize(NRows,ncols_res); 
    result.tol=tol;
    if(ncols_res>0)
      for(int j=0; j<ncols_res; j++)
	for(int i=0; i<NRows; i++)
	  result(i,j)=data2D[jcols(j)][i];	
    
    return result;
  };
  
  ///returns a copy of the matrix excluding 1 row whose index is stored in irow, this isn't an inline because you will need to copy all but 1 row
  SurfMat<T>& excludeRows(SurfMat<T>& result, int irow);
  
  ///returns a copy of the matrix that excludes all rows whose indices are stored in the matrix irows, this isn't an inline because for the typical use you will need to copy most of the matrix
  SurfMat<T>& excludeRows(SurfMat<T>& result, SurfMat<int>& irows);
  
  ///returns a copy of the matrix excluding 1 column whose index is stored in jcol, this isn't an inline because you will need to copy all but 1 column
  SurfMat<T>& excludeCols(SurfMat<T>& result, int jcol);
  
  ///return a copy of the matrix that excludes all columns whose indices are stored in the matrix jcols, this isn't an inline because for the typical use case you will need to copy most of the matrix
  SurfMat<T>& excludeCols(SurfMat<T>& result, SurfMat<int>& jcols);
  
  ///used in the generic matrix quicksort; returns -1 if a-b<-tol, 0 if -tol<=a-b<=tol, +1 if tol<a-b
  inline int compareElemAElemB(int ia, int ib) {
    T diff=data[ia]-data[ib];
    return(-(diff<-tol)+(tol<diff));
  }
  
  ///used in the generic matrix quicksort; swap the ia-th and ib-th elements of the "vector"
  inline void swapElems(int ia, int ib) {
    T swap=data[ia];
    data[ia]=data[ib];
    data[ib]=swap;
    return;
  };
  
  ///ascending sort the elements of the matrix as if it were a vector, could easily add an integer argument (+1 or -1) that indicates whether you want to sort into ascending or descending order, would need to add a field to the matrix though
  inline void sortElems() {
    qsortElems(0,getNElems()-1);
    return;
  };
  

  void qsortElems(int istart, int istop);
  void qsortRows(int istart, int istop);
  void qsortCols(int istart, int istop);

  ///performs an ascending unique sort of the elements, eliminates duplicates, and reshapes it into a (shruken, if appropriate) vector
  inline void uniqueElems() {
    int nelems=getNElems();
    if(nelems>0) {
      sortElems();    
      int i,k;
      i=1; 
      while(i<nelems) {
	if(data[i]==data[i-1]) {
	  for(k=i+1; k<nelems; k++)
	    data[k-1]=data[k];
	  nelems--;}
	else i++;
      }
      reshape(nelems);
    }
    return;
  };
  
  ///used in the generic matrix quicksort; compares 2 rows (a and b) starting with column 0; it returns -1 if a-b<-tol, or +1 if tol<a-b, if -tol<=a-b<=tol it moves onto column 1, then 2 and so on, iff every column of the 2 rows are equal then it returns 0, could add a column order field to the matrix to add the capability to change the order in which the rows are sorted and whether or not you wanted to do an ascending or descending sort
  inline int compareRowARowB(int irowa, int irowb) {
    T diff;
    int j=0;
    do{
      diff=data2D[j][irowa]-data2D[j][irowb];
      j++;
    }while((-tol<=diff)&&(diff<=tol)&&(j<NCols));
    return (-(diff<-tol)+(tol<diff));  //could multiply outermost () by +1/-1 for ascending/descending
  };
  
  ///used in the generic matrix quicksort; swap rows a and b of the matrix
  inline void swapRows(int irowa, int irowb) {
    T swap;
    for(int j=0; j<NCols; j++) {
      swap=data2D[j][irowa];
      data2D[j][irowa]=data2D[j][irowb];
      data2D[j][irowb]=swap;
    }
    return;
  };
  
  ///ascending sort the rows of the matrix (comparing from left to right), you could easily add a MtxInt argument with 2 rows the first indicating order of columns being compared, the second being +1 or -1 to indicate to sort that row into ascending or descending order, would need to add a field to the matrix though
  inline void sortRows() {
    qsortRows(0,NRows-1);
    return;
  };
  
  ///performs an ascending sort of the rows, eliminates duplicates, and shrinks the matrix if it is appropriate to do so
  inline void uniqueRows() {
    int nrows=NRows;
    if(nrows>0) {
      sortRows();
      int i,j,k;
      i=1;
      while(i<nrows) {
	if(compareRowARowB(i,i-1)) {
	  for(j=0;j<NCols;j++) 
	    for(k=i+1;k<nrows;k++)
	      data2D[j][k-1]=data2D[j][k];
	  nrows--;}
	else i++;
      }
      resize(nrows,NCols); //not an error, want non-contiguous memory chopping
    }
    return;
  };
  
  ///used in the generic matrix quicksort; compares 2 columns starting with row 0; it returns -1  if a-b<-tol, or +1 if tol<a-b, if -tol<=a-b<=tol it moves onto the next row, iff every row of the 2 columns have -tol<=a-b<=tol then it returns 0, could add a row order field to the matrix to add the capability to change the order in which the columns are sorted and whether or not you wanted to do an ascending or descending sort
  inline int compareColAColB(int jcola, int jcolb) {
    T diff;
    int i=0;
    do{
      diff=data2D[jcola][i]-data2D[jcolb][i];
      i++;
    }while((-tol<=diff)&&(diff<=tol)&&(i<NRows));
    return (-(diff<-tol)+(tol<diff));
  };

  ///used in the generic matrix quicksort; swap columns a and b of the matrix
  inline void swapCols(int jcola, int jcolb) {
    T swap;
    for(int i=0; i<NRows; i++) {
      swap=data2D[jcola][i];
      data2D[jcola][i]=data2D[jcolb][i];
      data2D[jcolb][i]=swap;
    }
    return;
  };
  
  ///ascending sort the coluns of the matrix (comparing from row zero onward), you could easily add a MtxInt argument with 2 columns the first indicating order of columns being compared (some rows can be omitted if you don't want to include them as sorting criteria), the second column being +1 or -1 to indicate to sort that row into ascending or descending order, would need to add a field to the matrix though
  inline void sortCols() {
    qsortCols(0,NCols-1);
    return;
  };
  
  ///performs an ascending sort of the columns, eliminates duplicates, and shrinks the matrix if it is appropriate to do so
  inline void uniqueCols() {
    int ncols=NCols;
    if(ncols>0) {
      sortCols();
      int i,j,k;
      j=1;
      while(j<ncols) {
	if(compareColAColB(j,j-1)) {
	  for(k=j+1;k<ncols;k++) 
	    for(i=0;i<NRows;i++)
	      data2D[k-1][i]=data2D[k][i];
	  ncols--;}
	else j++;
      }
      reshape(NRows,ncols); //not an error, want contiguous memory chopping, reshape will call reshape2, resize would call resize2 which would call reshape2  
    }
    return;
  };

private:
  ///newSize2 should not be called by anything except newSize
  void newSize2(int nrows_new, int ncols_new);
  //reshape2 should not be called by anything but reshape and resize2
  void reshape2(int nrows_new, int ncols_new);
  //resize2 should not be called by anything but resize
  void resize2(int nrows_new, int ncols_new);
  int NRows;
  int NCols;
  T  *data;
  T **data2D;
  T   tol; //an inequaltiy tolerance for equality checking should be 0 for integers, WARNING TO WHOEVER FINISHES IMPLEMENTING THIS, you will have to decide what to do with tol when you clear(), newsize(), reshape(), or resize() the matrix, you will likely need to modify the contructors
  
  //this shouldn't be necessary but do it for good measure
  friend void mtxqsort(int istart, int istop,
		       int (* compare_a_b)(int ia, int ib),
		       void (* swap_a_b)(int ia, int ib));

};



//typical matrix constructor
template< typename T >
SurfMat<T>::SurfMat(int nrows_in, int ncols_in)
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((nrows_in>=0)&&(ncols_in>=0));
#endif
  tol=0;
  NRows=NCols=0;
  data=NULL;
  data2D=NULL;
  if((nrows_in>0)&&(ncols_in>0)) {
    data  =(T*) data_alloc(nrows_in*ncols_in,sizeof(T));
    data2D=(T**)malloc(ncols_in*sizeof(T *));
#ifdef __SURFMAT_ERR_CHECK__
    assert(data && data2D);
#endif
    NRows=nrows_in;
    NCols=ncols_in;
    *data2D=data;
    for(int j=1; j<NCols; j++)
      data2D[j]=data2D[j-1]+NRows;    
  }
  return;
}

//copy constructor, make a deep copy
template< typename T >
SurfMat<T>::SurfMat(const SurfMat<T>& other){
  //NRows=NCols=0;
  data  =NULL;
  data2D=NULL;
  tol=other.tol;
  NRows=other.NRows;
  NCols=other.NCols;
  if(!((NRows>0)&&(NCols>0))) {
    return;
    //printf("stop me\n");
    //assert((NRows>0)&&(NCols>0)); //for debug only
  }

  int k, nelem=NRows*NCols;
  data  =(T*) malloc(nelem*sizeof(T)); //we're going to fill it up completely so don't bother with data_alloc (which gives an option to calloc i.e. zero memory)
  data2D=(T**)malloc(NCols*sizeof(T*));
#ifdef __SURFMAT_ERR_CHECK__
  assert(data && data2D);
#endif
  *data2D=data;
  for(k=1; k<NCols; k++)
    data2D[k]=data2D[k-1]+NRows; 

  //for(k=0; k<nelem; k++)
  //data[k]=other.data[k];

  const T* od=other.data;
  for(k=0; k<nelem; k++)
    data[k]=od[k];

  return;
}

//copy from vector constructor
template< typename T >
SurfMat<T>::SurfMat(const std::vector<T>& vecT, int nrows_in, int ncols_in)
{
  tol=0;
  NRows=NCols=0;
  data=NULL;
  data2D=NULL;
  
  int size_vecT=vecT.size();

  if(nrows_in==-99999) 
    nrows_in=size_vecT;

#ifdef __SURFMAT_ERR_CHECK__
  assert((nrows_in>=0)&&(ncols_in>=0));
#endif

  if((nrows_in>0)&&(ncols_in>0)) {
    data  =(T*) data_alloc(nrows_in*ncols_in,sizeof(T));
    data2D=(T**)malloc(ncols_in*sizeof(T *));
#ifdef __SURFMAT_ERR_CHECK__
    assert(data && data2D);
#endif
    NRows=nrows_in;
    NCols=ncols_in;
    *data2D=data;
    for(int j=1; j<NCols; j++)
      data2D[j]=data2D[j-1]+NRows;    
  }  

  int n_to_copy=NRows*NCols;
  if(n_to_copy>size_vecT)
    n_to_copy=size_vecT;
  
  for(int k=0; k<n_to_copy; ++k)
    data[k]=vecT[k];

  return;
}


//deallocate the memory
template< typename T >
void SurfMat<T>::clear()
{
#ifdef __SURFMAT_ERR_CHECK__
  if(!(((NRows==0)&&(NCols==0)&&(data==NULL)&&(data2D==NULL))||
       (NRows && NCols && data && data2D))) {
    printf("NRows=%d NCols=%d data=%d data2D=%d COMMENTED OUT",NRows,NCols,data,data2D);
    assert(((NRows==0)&&(NCols==0)&&(data==NULL)&&(data2D==NULL))||
	   (NRows && NCols && data && data2D));
  }
#endif
  if(NRows) {
    free(data2D); data2D=NULL;
    free(data  ); data  =NULL;
    NRows=NCols=0;
  }
  return;
}

//enlarge or shrink the matrix without preserving it's contents, this is here for speed (resize2 and reshape2 are slower because they involve copying data, inherently done by realloc when it can't grab extra space at the end of the array) 
template< typename T >
void SurfMat<T>::newSize2(int nrows_new, int ncols_new)
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((nrows_new>=0)&&(ncols_new>=0));
#endif
  int nelem_new=nrows_new*ncols_new;
  if(nelem_new==0.0) {
    clear();
    return;
  }
  
  //if NRows*NCols==nelem_new then we only need to change data2D ptrs into data
  int nelem=NRows*NCols;
  //if(nelem!=nelem_new) {
  if(nelem<nelem_new) {
    free(data);
    data=(T*) data_alloc(nelem_new,sizeof(T));
  }
  else if(nelem>nelem_new) { 
    //decreasing the amount of memory allocated _may_be_ faster than deallocating the array and allocating a new one (couldn't find anything that said one way or another so just time it)
    data=(T*) realloc(data,nelem_new*sizeof(T));
#ifdef __SURFMAT_ZERO_MEM__
    for(int j=0; j<nelem_new; j++) data[j]=0;
#endif
  }
  
  NRows=nrows_new;
  NCols=ncols_new;
    
  //change data2D pointers into data
  free(data2D); 
  data2D=(T**) malloc(NCols*sizeof(T*)); 
#ifdef __SURFMAT_ERR_CHECK__
  assert(data2D);
#endif
  *data2D=data;
  for(int j=1; j<NCols; j++)
    data2D[j]=data2D[j-1]+NRows;
  
  return;
}

//enlarge or shrink the matrix while keeping contigous in memory elements contiguous in memory: reshape2() should not be called by anything but reshape() and resize2()
template< typename T >
void SurfMat<T>::reshape2(int nrows_new, int ncols_new)
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((nrows_new>=0)&&(ncols_new>=0));
#endif
  int nelem_new=nrows_new*ncols_new;
  if(nelem_new==0) {
    clear();
    return;
  }
  
  //if NRows*NCols==nelem_new then we only need to change data2D ptrs into data
  int nelem=NRows*NCols;
  if(nelem!=nelem_new) {
    data=(T*) realloc(data,nelem_new*sizeof(T));
#ifdef __SURFMAT_ZERO_MEM__
    for(int j=nelem; j<nelem_new; j++) data[j]=0;
#endif      
  }
  
  NRows=nrows_new;
  NCols=ncols_new;
  
  //change data2D pointers into data
  free(data2D); 
  data2D=(T**) malloc(NCols*sizeof(T*)); 
#ifdef __SURFMAT_ERR_CHECK__
  assert(data2D);
#endif
  *data2D=data;
  for(int j=1; j<NCols; j++)
    data2D[j]=data2D[j-1]+NRows;
  
  return;
}

//enlarge or shrink the matrix while adding "zeros" after the last row and/or last column and/or chopping off the rows and/or columns after the newly requested last row and/or column, behavior "should be" the same as Teuchos' Serial Dense Matrix's reshape function, this function should only be called by resize()
template< typename T >
void SurfMat<T>::resize2(int nrows_new, int ncols_new) {
#ifdef __SURFMAT_ERR_CHECK__
  assert((nrows_new>=0)&&(ncols_new>=0));
#endif
  int nelem_new=nrows_new*ncols_new;
  if(nelem_new==0) {
    clear();
    //printf("resize2 calling clear\n"); fflush(stdout);
    return;
  }
  
  if(NRows==nrows_new){
    //only need to chop off or extend the last few columns (which are contiguous) so reshape2 can do this
    reshape2(nrows_new,ncols_new);
    //printf("resize2 calling reshape2\n"); fflush(stdout);
    return;
  }

  //printf("resize2 doing the resizing itself\n"); fflush(stdout);
    
  int i, j;
  T **data2D_temp=(T **) malloc(ncols_new*sizeof(T*));
  T  *data_temp=(T *) data_alloc(nrows_new*ncols_new,sizeof(T));
  
  *data2D_temp=data_temp;
  for(j=1; j<ncols_new; j++)
    data2D_temp[j]=data2D_temp[j-1]+nrows_new;
  
  int nrows_smaller=(NRows<nrows_new)?NRows:nrows_new;
  int ncols_smaller=(NCols<ncols_new)?NCols:ncols_new;
  NRows=nrows_new;
  NCols=ncols_new;
  
  for(j=0; j<ncols_smaller; j++)
    for(i=0; i<nrows_smaller; i++)
      data2D_temp[j][i]=data2D[j][i];

  free(data2D); data2D=data2D_temp;
  free(data  ); data  =data_temp;
  
  return;
}

///Perform a deep copy
template< typename T >
SurfMat<T>& SurfMat<T>::copy(const SurfMat<T>& other){
  
  int k, nelem=NRows*NCols, nelem_o=other.NRows*other.NCols;
  
  if((NRows!=other.NRows)||(NCols!=other.NCols)) {
    if(NCols!=other.NCols) {
      free(data2D);
      data2D=(T**)malloc(other.NCols*sizeof(T*));
    }
    NRows=other.NRows;
    NCols=other.NCols;
    
    if(nelem!=nelem_o) {
      free(data);
      data=(T* ) malloc(nelem_o*sizeof(T));  // don't give this the option to fill with zeros because we are about to copy other to fill in all elements of data
    }

#ifdef __SURFMAT_ERR_CHECK__
    assert(data && data2D);
#endif
    *data2D=data;
    for(k=1; k<NCols; k++)
      data2D[k]=data2D[k-1]+NRows;
    //for(int j=1; j<NCols; j++)
    //data2D[j]=data2D[j-1]+NRows;    
  }

  tol=other.tol;
  
  const T *od=other.data;
  for(k=0; k<nelem_o; k++) 
    data[k]=od[k];

  //for(int ij=0; ij<nelem_o; ij++) 
  //data[ij]=other.data[ij];
  
  return *this;
}


///returns a copy of the matrix excluding 1 row whose index is stored in irow, this isn't an inline because you will need to copy all but 1 row
template< typename T >
SurfMat<T>& SurfMat<T>::excludeRows(SurfMat<T>& result, int irow) {
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=irow)&&(irow<NRows));
#endif
  if(NRows==1)
    result.clear();
  else{
    result.newSize(NRows-1,NCols);
    result.tol=tol;
    int isrc, ikeep, j;
    for(j=0; j<NCols; j++) {      
      for(isrc=ikeep=0; isrc<irow; isrc++, ikeep++)
	result.data2D[j][ikeep]=data2D[j][isrc];
      isrc=irow+1;
      for(;isrc<NRows;isrc++,ikeep++)
	result.data2D[j][ikeep]=data2D[j][isrc];
    }
  }
  return result;
}

///returns a copy of the matrix excluding all rows whose indices are stored in matrix irows, this isn't an inline because for the typical use you will need to copy most of the matrix
template< typename T >
SurfMat<T>& SurfMat<T>::excludeRows(SurfMat<T>& result, SurfMat<int>& irows) {
  int j;
  int nexclude=irows.getNElems();
  if(nexclude<1) {
    //the list of rows to exclude is empty so copy over the whole matrix
    result.newSize(NRows,NCols);
    result.tol=tol;
    for(j=0;j<NCols;j++)
      for(int i=0;i<NRows;i++)
	result.data2D[j][i]=data2D[j][i];
  }
  else{
    irows.uniqueElems(); //sort the rows to exclude into ascending order and eliminate duplicate listings
    nexclude=irows.getNElems();
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=irows(0))&&(irows(nexclude-1)<NRows));
#endif
    if(nexclude==NRows) {
      //the user wants us to eliminate _all_ rows
      result.clear();}
    else{
      //the user wants us to eliminate some but not all rows
      result.newSize(NRows-nexclude,NCols);
      result.tol=tol;
      int iexclude, ikeep, isrc;
      for(j=0;j<NCols;j++) {
	iexclude=ikeep=isrc=0;
	while(isrc<NRows) {
	  if(iexclude<nexclude) {
	    for(;isrc<irows(iexclude);isrc++, ikeep++)
	      result.data2D[j][ikeep]=data2D[j][isrc];
	    //at this point isrc=irows(iexclude)
	    iexclude++;
	    isrc++;}
	  else{
	    for(;isrc<NRows;isrc++, ikeep++)
	      result.data2D[j][ikeep]=data2D[j][isrc];
	    //at this point isrc=NRows and the while loop will terminate
	  }
	}//while loop terminates
      }//do the same thing with the next column  
    }
  }
  return result;
}

///returns a copy of the matrix excluding 1 column whose index is stored in jcol, this isn't an inline because you will need to copy all but 1 column
template< typename T >
SurfMat<T>& SurfMat<T>::excludeCols(SurfMat<T>& result, int jcol) {
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=jcol)&&(jcol<NCols));
#endif
  if(NCols==1)
    result.clear();
  else{
    result.newSize(NRows,NCols-1);
    result.tol=tol;
    int jsrc, jkeep, i;
    for(jsrc=jkeep=0; jsrc<jcol; jsrc++, jkeep++)
      for(i=0; i<NRows; i++)
	result.data2D[jkeep][i]=data2D[jsrc][i];
    jsrc++;
    for(; jsrc<NRows; jsrc++, jkeep++)
      for(i=0; i<NRows; i++)
	result.data2D[jkeep][i]=data2D[jsrc][i];
  }
  return result;
}

///return a copy of the matrix that excludes all columns whose indices are stored in the matrix jcols, this isn't an inline because for the typical use case you will need to copy most of the matrix
template< typename T >
SurfMat<T>& SurfMat<T>::excludeCols(SurfMat<T>& result, SurfMat<int>& jcols) {
  int nexclude=jcols.getNElems();
  int i;
  if(nexclude<1) {
    result.newSize(NRows,NCols);
    result.tol=tol;
    for(int j=0;j<NCols;j++)
      for(i=0;i<NRows;i++)
	result.data2D[j][i]=data2D[j][i];
  }
  else{
    jcols.uniqueElems();
    nexclude=jcols.getNElems();
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=jcols(0))&&(jcols(nexclude-1)<NCols));
#endif
    if(nexclude==NCols) {
      //the user wants us to eliminate _all_ columns
      result.clear(); }
    else {
      //the user wants us to eliminate some but not all columns
      result.newSize(NRows,NCols-nexclude);
      result.tol=tol;
      int jexclude, jkeep, jsrc;
      jexclude=jkeep=jsrc=0;
      while(jsrc<NCols) {
	if(jexclude<nexclude) {
	  for(;jsrc<jcols(jexclude);jsrc++,jkeep++)
	    for(i=0;i<NRows;i++)
	      result.data2D[jkeep][i]=data2D[jsrc][i];
	  //at this point jsrc=jrows(jexclude)
	  jexclude++;
	  jsrc++;
	}
	else{
	  for(;jsrc<NCols;jsrc++,jkeep++)
	    for(i=0;i<NRows;i++)
	      result.data2D[jkeep][i]=data2D[jsrc][i];
	  //at this point jsrc=NCols and the while loop will terminate
	}
      }//while loop terminates
    }
  }
  return result;
}

template< typename T >
void SurfMat<T>::qsortElems(int istart, int istop)
{
  int i,j,k;
  if( istart < istop) {
    k = (istart+istop)/2;
    swapElems(istart,k);
    i = istart+1;
    j = istop;
    while(i <= j){
      while((i <= istop) && (compareElemAElemB(i,istart)<=0))
	i++;
      while((j > istart) && (compareElemAElemB(istart,j)<0))
	j--;
      if( i < j)
	swapElems(i,j);
    }
    // swap two elements
    swapElems(istart,j);
    // recursively sort the lesser list
    qsortElems(istart,j-1);
    qsortElems(j+1,istop);
  }
  return;
}
  
template< typename T >
void SurfMat<T>::qsortRows(int istart, int istop)
{
  int i,j,k;
  if( istart < istop) {
    k = (istart+istop)/2;
    swapRows(istart,k);
    i = istart+1;
    j = istop;
    while(i <= j){
      while((i <= istop) && (compareRowARowB(i,istart)<=0))
	i++;
      while((j > istart) && (compareRowARowB(istart,j)<0))
	j--;
      if( i < j)
	swapRows(i,j);
    }
    // swap two elements
    swapRows(istart,j);
    // recursively sort the lesser list
    qsortRows(istart,j-1);
    qsortRows(j+1,istop);
  }
  return;
}

template< typename T >
void SurfMat<T>::qsortCols(int istart, int istop)
{
  int i,j,k;
  if( istart < istop) {
    k = (istart+istop)/2;
    swapCols(istart,k);
    i = istart+1;
    j = istop;
    while(i <= j){
      while((i <= istop) && (compareColAColB(i,istart)<=0))
	i++;
      while((j > istart) && (compareColAColB(istart,j)<0))
	j--;
      if( i < j)
	swapCols(i,j);
    }
    // swap two elements
    swapCols(istart,j);
    // recursively sort the lesser list
    qsortCols(istart,j-1);
    qsortCols(j+1,istop);
  }
  return;
}

///returns the smallest value contained in this matrix, wraps T& minElem(T& val ,int& loc)
template< typename T >
inline T SurfMat<T>::minElem() const {
  T val;
  int loc;
  minElem(val,loc);
  return val;
}

///returns the location of the smallest value contained in this matrix, if there is a tie it returns the location of the first one encountered, wraps T& minElem(T& val ,int& loc)
template< typename T >
inline int SurfMat<T>::iMinElem() const {
  T val;
  int loc;
  minElem(val,loc);
  return loc;
}

///returns the smallest value contained in this matrix and its location in loc, if there is a tie it records the location of the first smallest value encountered
template< typename T >
inline void SurfMat<T>::minElem(T& val, int& loc) const {
  int nelem=NRows*NCols; 
#ifdef __SURFMAT_ERR_CHECK__
  assert(nelem>0);
#endif
  val=data[0];
  loc=0;
  for(int k=1; k<nelem; k++)
    if(data[k]<val){
      val=data[k];
      loc=k;}
  return;
}


///returns the largest value contained in this matrix, wraps T& maxElem(T& val ,int& loc)
template< typename T >
inline T SurfMat<T>::maxElem() const {
  T val;
  int loc;
  maxElem(val,loc);
  return val;
}

///returns the location of the largest value contained in this matrix, if there is a tie it returns the location of the first one encountered, wraps T& maxElem(T& val ,int& loc)
template< typename T >
inline int SurfMat<T>::iMaxElem() const {
  T val;
  int loc;
  maxElem(val,loc);
  return loc;
}

///returns the largest value contained in this matrix and its location in loc, if there is a tie it records the location of the first largest value encountered
template< typename T >
inline void SurfMat<T>::maxElem(T& val, int& loc) const {
  int nelem=NRows*NCols; 
#ifdef __SURFMAT_ERR_CHECK__
  assert(nelem>0);
#endif
  val=data[0];
  loc=0;
  for(int k=1; k<nelem; k++)
    if(data[k]>val){
      val=data[k];
      loc=k;}
  return;
}



/* //incase we decide to make these not inline functions
   
///performs an ascending unique sort of the elements, eliminates duplicates, and reshapes it into a (shruken, if appropriate) vector, this isn't an inline because of the work you need to do to eliminate duplicates
void SurfMat<T>::uniqueElems() {
int nelems=getNElems();
if(nelems>0) {
sortElems();    
int i,k;
    i=1; 
    while(i<nelems) {
      if(data[i]==data[i-1]) {
	for(k=i+1; j<nelems; k++)
	  data[k-1]=data[k];
	nelems--;}
      else i++;
    }
    reshape(nelems);
  }
  return;
}

///performs an ascending sort of the rows, eliminates duplicates, and shrinks the matrix if it is appropriate to do so, this isn't an inline because of the work you need to do compare rows (which is an inline) and to eliminate duplicates
void SurfMat<T>::uniqueRows() {
  int nrows=NRows;
  if(nrows>0) {
    sortRows();
    int i,j,k;
    i=1;
    while(i<nrows) {
      if(compareRowARowB(i,i-1)==0) {
	for(j=0;j<NCols;j++) 
	  for(k=i+1;k<nrows;k++)
	    data2D[j][k-1]=data2D[j][k];
	nrows--;}
      else i++;
    }
    resize(nrows,NCols); //not an error, want non-contiguous memory chopping
  }
  return;
}

///performs an ascending sort of the columns, eliminates duplicates, and shrinks the matrix if it is appropriate to do so, this isn't an inline because of the work you need to do compare columns (which is an inline) and to eliminate duplicates
void SurfMat<T>::uniqueCols() {
  int ncols=NCols;
  if(ncols>0) {
    sortCols();
    int i,j,k;
    j=1;
    while(j<ncols) {
      if(compareColAColB(j,j-1)==0) {
	for(k=j+1;k<ncols;k++) 
	  for(i=0;i<NRows;i++)
	    data2D[k-1][i]=data2D[k][i];
	ncols--;}
      else j++;
    }
    reshape(NRows,ncols); //not an error, want contiguous memory chopping, reshape will call reshape2, resize would call resize2 which would call reshape2  
  }
  return;
}

*/


} // end namespace nkm


#endif
