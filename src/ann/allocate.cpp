#include "surfpack_config.h"

/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

/*
 * File: allocate.c
 * -----------------
 * This file allows for implementation of the utilites.h library interface.
 * Although this file was developed to make the running of the MAIN
 * file, these utility files are arranged in NO order. 
 */

#include "system_defs.h"
#include "allocate.h" 

using namespace std;

/* nrutil.c */
void nrerror(const char* error_text)
{
   //exit();

   fprintf(stderr,"Numerical Recipes run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n");
   exit(1);
}

float *_vector(int nl,int nh)
{
   float *v;

   v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
   v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
   if (!v) nrerror("allocation failure in vector()");
   return v-nl;
}

int *i_vector(int nl,int nh)
{
   int *v;

   v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
   if (!v) nrerror("allocation failure in i_vector()");
   return v-nl;
}

double *d_vector(int nl,int nh)
{
   double *v;

   v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
   if (!v) nrerror("allocation failure in d_vector()");
   return v-nl;
}

float **matrix(int nrl,int nrh,int ncl,int nch)
{
   int i;
   float **m;

   m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
   if (!m) nrerror("allocation failure 1 in matrix()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
      m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
      if (!m[i]) nrerror("allocation failure 2 in matrix()");
      m[i] -= ncl;
   }
   return m;
}

double **d_matrix(int nrl,int nrh,int ncl,int nch)
{
   int i;
   double **m;

   m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
   if (!m) nrerror("allocation failure 1 in d_matrix()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
      m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
      if (!m[i]) nrerror("allocation failure 2 in d_matrix()");
      m[i] -= ncl;
   }
   return m;
}

int **i_matrix(int nrl,int nrh,int ncl,int nch)
{
   int i,**m;

   m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
   if (!m) nrerror("allocation failure 1 in i_matrix()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
      m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
      if (!m[i]) nrerror("allocation failure 2 in i_matrix()");
      m[i] -= ncl;
   }
   return m;
}

float **submatrix(float **a,int oldrl,int oldrh,int oldcl,int oldch,int newrl,int newcl)
{
   int i,j;
   float **m;

   m=(float **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(float*));
   if (!m) nrerror("allocation failure in submatrix()");
   m -= newrl;

   for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

   return m;
}

void free_vector(float *v,int nl,int nh)
{
   free((char*) (v+nl));
}

void free_i_vector(int *v,int nl,int nh)
{

   free((char*) (v+nl));
}

void free_d_vector(double *v,int nl,int nh)
{
   free((char*) (v+nl));
}

void free_matrix(float **m,int nrl,int nrh,int ncl,int nch)
{
   int i;

   for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

void free_d_matrix(double **m,int nrl,int nrh,int ncl,int nch)
{
   int i;

   for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

void free_i_matrix(int **m,int nrl,int nrh,int ncl,int nch)
{
   int i;

   for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

void free_submatrix(float **b,int nrl,int nrh,int ncl,int nch)
{
   free((char*) (b+nrl));
}

float **convert_matrix(float *a,int nrl,int nrh,int ncl,int nch)
{
   int i,j,nrow,ncol;
   float **m;

   nrow=nrh-nrl+1;
   ncol=nch-ncl+1;
   m = (float **) malloc((unsigned) (nrow)*sizeof(float*));
   if (!m) nrerror("allocation failure in convert_matrix()");
   m -= nrl;
   for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;

   return m;
}

void free_convert_matrix(float **b,int nrl,int nrh,int ncl,int nch)
{
   free((char*) (b+nrl));
}


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

double ***dtensor3(int n1l,int n1h,int n2l,int n2h,int n3l,int n3h)
{
   int i, j;
   double ***t;

   /*  Allocate pointers to elements in first dimension.  */
   t = (double ***) malloc( (unsigned) (n1h-n1l+1)*sizeof(double**) );
   if( !t )
   {  printf("\n***  Dimension 1 alloc failure in tensor3().  ***\n");  }
   t -= n1l;

   /*  Allocate pointers to elements in second dimension. */
   for( i=n1l; i<=n1h; i++ )
   {  t[i] = (double **) malloc((unsigned)(n2h-n2l+1)*sizeof(double*));
      if( !t[i] )
      { printf("\n***  Dim 2 alloc failure in matrix().  ***\n"); }
      t[i] -= n2l;
   }

   /*  Allocate data storage in third dimension and pointers to it.  */
   for( i=n1l; i<=n1h; i++ )
   { for( j=n2l; j<=n2h; j++ )
     { t[i][j] = (double *) malloc((unsigned)(n3h-n3l+1)*sizeof(double));
       if( !t[i][j] )
       { printf("\n***  Dim 3 alloc failure in matrix().  ***\n"); }
       t[i][j] -= n3l;
     }
   }

   /*  Return pointer to pointer to array of pointers to stroage vectors.*/
   return t;
}

void free_dtensor3(double ***t,int n1l,int n1h,int n2l,int n2h,int n3l, int n3h)
{
        int i, j;

        for( i=n1h; i>=n1l; i-- )
   {  for( j=n2h; j>=n2l; j-- )
         {       free( (char*) (t[i][j]+n3l) );     }
      free( (char*) (t[i]+n2l) );
   }
        free( (char*) (t+n1l) );
}

/* lmatrix.c */
/* This file is stored in lmatrix.c.
   The first function can be used to allocate a matrix of long
   integer variables.
   The second function can be used to deallocate the matrix.
*/

long int **lmatrix(int nrl,int nrh,int ncl,int nch)
{
    int i;
   long int **m;

    m=(long int **)malloc((unsigned) (nrh-nrl+1)*sizeof(long int*));
    if (!m) nrerror("allocation failure 1 in lmatrix()");
    m -= nrl;

    for(i=nrl;i<=nrh;i++) {
        m[i]=(long int *)malloc((unsigned) (nch-ncl+1)*sizeof(long int));
        if (!m[i]) nrerror("allocation failure 2 in lmatrix()");
        m[i] -= ncl;
    }
    return m;
}

void free_lmatrix(long int **m,int nrl,int nrh,int ncl,int nch)
{
    int i;

    for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
    free((char*) (m+nrl));
}
