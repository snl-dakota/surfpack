#include "surfpack_config.h"

/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

/*	This file is stored in dsvd2.c.

	Its purpose is to perform an svd on input data.  It does this by 
	accepting the data, performing the svd in dsvdcmp1(), reordering
	the svd in reorder_dsvd().  The function reorder_dsvd() reorders
	the singular values in descending order, and reorders the 
	columns of the u and v matrices accordingly.

	To use this function one must link dsvdcmp1.o, reorder_dsvd.o, dindexx.o,
	and nrutil.o.  The latter file is required because of the vector
	generation in the first three functions.  The file dindexx.o performs
	the reordering of the singular values.

	Inputs:
	**a	The matrix whose svd is to be computed. (double)
	m	The number of rows in a. (int)
	n	The number of columns in a. (int)

	Outputs:
	**u	The coeficient matrix of the svd. (double)
	*w	The vector of singular values. (double)
	**v	The matrix of singular vectors. (double)
*/

#include "system_defs.h"
#include "allocate.h"
#include "dsvd2.h"
#ifdef HAVE_STD
#include <cmath>
#else
#include <math.h>
#endif

using namespace std;


void dsvd2(double **a, int m, int n, double **u,double *w,double **v)
{
	int i, j;

	/*  Move a into u.  */
	for( i=0; i<m; i++ )
	{	for( j=0; j<n; j++ )
		{	u[i][j] = a[i][j];	}
	}

	/*  Perform the svd.  */
	dsvdcmp1(u,m,n,w,v);

	/*  Reorder the svd.  */
	reorder_dsvd(u,w,v,m,n);
}

/***************  dsvdcmp1.c  *********************/
static double at,bt,ct;

#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static double maxarg1,maxarg2;

#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
    (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void dsvdcmp1(double **a, int m,int n,double *w,double **v)
{
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm=0.0, g=0.0, scale=0.0;
    double *rv1;

    if (m < n) nrerror("SVDCMP: You must augment A with extra zero rows");
    rv1 = d_vector(0,n-1);
    for ( i=0; i<n; i++)
    {
        l = i+1;
        rv1[i] = scale*g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for ( k=i; k<m; k++ ) scale += fabs(a[k][i]);
            if ( scale )
            {
                for ( k=i; k<m; k++)
                {
                    a[k][i] /= scale;
                    s += a[k][i]*a[k][i];
                }
                f = a[i][i];
                g = -SIGN( sqrt(s), f );
                h = f*g-s;
                a[i][i] = f-g;
                if (i != (n-1))
                {
                    for ( j=l; j<n; j++ )
                    {
                        for ( s=0.0, k=i; k<m; k++) s += a[k][i]*a[k][j];
                        f = s/h;
                        for ( k=i; k<m; k++ ) a[k][j] += f*a[k][i];
                    }
                }
                for ( k=i; k<m; k++ ) a[k][i] *= scale;
            }
        }
        w[i] = scale*g;
        g = s = scale = 0.0;
        if (i < m && i != (n-1))
        {
            for ( k=l; k<n; k++ ) scale += fabs(a[i][k]);
            if ( scale )
            {
                for ( k=l; k<n; k++ )
                {
                    a[i][k] /= scale;
                    s += a[i][k]*a[i][k];
                }
                f = a[i][l];
                g = -SIGN( sqrt(s), f );
                h = f*g-s;
                a[i][l] = f-g;
                for ( k=l; k<n; k++ ) rv1[k] = a[i][k]/h;
                if (i != (m-1))
                {
                    for ( j=l; j<m; j++ )
                    {
                        for ( s=0.0, k=l; k<n; k++ ) s += a[j][k]*a[i][k];
                        for ( k=l; k<n; k++ ) a[j][k] += s*rv1[k];
                    }
                }
                for ( k=l; k<n; k++ ) a[i][k] *= scale;
            }
        }
        anorm =MAX( anorm, (fabs(w[i])+fabs(rv1[i])) );
    }
    for ( i=(n-1); i>=0; i-- )
    {
        if (i < (n-1))
        {
            if (g)
            {
                for ( j=l; j<n; j++ )
                    v[j][i] = (a[i][j]/a[i][l])/g;
                for ( j=l; j<n; j++)
                {
                    for ( s=0.0, k=l; k<n; k++ ) s += a[i][k]*v[k][j];
                    for ( k=l; k<n; k++ ) v[k][j] += s*v[k][i];
                }
            }
            for ( j=l; j<n; j++ ) v[i][j]=v[j][i]=0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    for ( i=(n-1); i>=0; i-- )
    {
        l = i+1;
        g = w[i];
        if (i < (n-1))
            for ( j=l; j<n; j++ ) a[i][j]=0.0;
        if (g)
        {
            g = 1.0/g;
            if (i != (n-1))
            {
                for ( j=l; j<n; j++ )
                {
                    for ( s=0.0, k=l; k<m; k++ ) s += a[k][i]*a[k][j];
                    f=(s/a[i][i])*g;
                    for ( k=i; k<m; k++ ) a[k][j] += f*a[k][i];
                }
            }
            for ( j=i; j<m; j++ ) a[j][i] *= g;
        }
        else
        {
            for ( j=i; j<m; j++ ) a[j][i]=0.0;
        }
        ++a[i][i];
    }
    for ( k=(n-1); k>=0; k-- )
    {
	
        for ( its=1; its<=30; its++ )
        {
	//	printf("Entered 30 iteration limit\n");
            flag = 1;
            for ( l=k; l>=0; l-- )
            {
 	//	printf("nm: %d\n",nm);
                nm = l-1;
                if ( fabs(rv1[l])+anorm == anorm )
                {
                    flag = 0;
                    break;
                }
                if ( fabs(w[nm])+anorm == anorm ) break;
            }
            if (flag) // line changed by Mark Richards
            {
		printf("Entered the flag zone\n");
                nm = l-1;
 		printf("nm: %d\n",nm);
 		if (nm < 0) nm = 0;
                c = 0.0;
                s = 1.0;
                for ( i=l; i<=k; i++ )
                {
                    f = s*rv1[i];
                    if ( fabs(f)+anorm != anorm )
                    {
                        g = w[i];
                        h = PYTHAG(f,g);
                        w[i] = h;
                        h = 1.0/h;
                        c = g*h;
                        s = (-f*h);
                        for ( j=0; j<m; j++ )
                        {
                            y = a[j][nm];
                            z = a[j][i];
                            a[j][nm] = y*c+z*s;
                            a[j][i] = z*c-y*s;
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {
                if (z < 0.0)
                {
                    w[k] = -z;
                    for ( j=0; j<n; j++ ) v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its == 30) nrerror("No convergence in 30 SVDCMP iterations");
            x = w[l];
            nm = k-1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g = PYTHAG(f,1.0);
            f= ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c = s = 1.0;
            for ( j=l; j<=nm; j++ )
            {
                i = j+1;
                g = rv1[i];
                y = w[i];
                h = s*g;
                g = c*g;
                z = PYTHAG(f,h);
                rv1[j] = z;
                c = f/z;
                s = h/z;
                f = x*c+g*s;
                g = g*c-x*s;
                h = y*s;
                y = y*c;
                for ( jj=0; jj<n; jj++ )
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x*c+z*s;
                    v[jj][i] = z*c-x*s;
                }
                z = PYTHAG(f,h);
                w[j] = z;
                if (z) {
                    z = 1.0/z;
                    c = f*z;
                    s = h*z;
                }
                f = (c*g)+(s*y);
                x = (c*y)-(s*g);
                for ( jj=0; jj<m; jj++ )
                {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y*c+z*s;
                    a[jj][i] = z*c-y*s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }
    free_d_vector( rv1, 0, n-1 );
}

#undef SIGN
#undef MAX
#undef PYTHAG


/***********************  reorder_dsvd.c ******************/
/*      This file is stored in reorder_dsvd.c.

        Its purpose is to reorder the terms in an svd so that the values
        in the singular values vector are in decreasing order, and the
        u and v matrices are rearranged accordingly.
*/
void reorder_dsvd(double **u,double *w,double **v,int m,int n)
{
        int i, j, *idx;
        double *temp;

        /*  Allocate arrays.  */
        idx = i_vector( 0, n-1 );
        temp = d_vector( 0, n-1 );

        /*  Find the index order.  */
        dindexx(n,w-1,idx-1);

        /*  Reorder u, w, v.  */
        for( i=0; i<m; i++ )
        {       for( j=0; j<n; j++ )
                {       temp[j] = u[i][j];      }
                for( j=0; j<n; j++ )
                {       u[i][j] = temp[idx[n-1-j]-1];   }
        }
        for( j=0; j<n; j++ )
        {       temp[j] = w[j]; }
        for( j=0; j<n; j++ )
        {       w[j] = temp[idx[n-1-j]-1];      }
        for( i=0; i<n; i++ )
        {       for( j=0; j<n; j++ )
                {       temp[j] = v[i][j];      }
                for( j=0; j<n; j++ )
                {       v[i][j] = temp[idx[n-1-j]-1];   }
        }

        /*  Deallocate idx.  */
        free_i_vector( idx, 0, n-1 );
        free_d_vector( temp, 0, n-1 );
}

/******************** dindexx.c *********************/
/*      This file is stored in dindexx.c.
        It indexes an array arrin[1...n], i.e. outputs the array indx[1...n]
        such that arrin[indx[j]] is in ascending order for j=1,2,...,n.  The
        input quantities,n and arrin are not changed.
*/

void dindexx(int n,double arrin[],int indx[])
{
        int l, j, ir, indxt, i;
        double q;

        /*  Initialize the index array with consecutive integers.  */
        for( j=1; j<=n; j++ ) indx[j] = j;
        if( n== 1 ) return;

        /*  From here on, we just have heapsort, but with indirect indexing
        through indx in all references to arrin.  */
        l = (n>>1) + 1;
        ir = n;
        for( ; ; )
        {
                if( l>1 )
                {       q = arrin[(indxt=indx[--l])];   }
                else
                {       q = arrin[(indxt=indx[ir])];
                        indx[ir] = indx[1];
                        if( --ir == 1 )
                        {       indx[1] = indxt;
                                return;
                        }
                }
                i = l;
                j = l<<1;
                while( j<=ir )
                {       if( (j<ir) && (arrin[indx[j]]<arrin[indx[j+1]]) ) j++;
                        if( q<arrin[indx[j]] )
                        {       indx[i] = indx[j];
                                j += (i=j);
                        }
                        else j = ir + 1;
                }
                indx[i] = indxt;
        }
}

