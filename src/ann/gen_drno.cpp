/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack_config.h"

/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

/*	This file is stored in gen_drno.c.
	Its purpose is to generate random numbers from various random 
	number sources.  At each call the parameter type must be specified,
	and p must point to an array of parameters with the contents of
	the parameter list specified below.
	Definitions:
		type	This identifies the distribution that is the source 
			of the random number.  Currently given by the 
			following list with corresponding parameters.

			u	Uniform variate.
				p[0]	Lower limit.
				p[1]	Upper limit.

			n	Normal variate.
				p[0]	Mean.
				p[1]	Standard deviation.

	The user must include stdlib.h and should initiate random number 
	generation with a call to srand().

	type is a character variable, and p points to an array of double.
*/

#include "system_defs.h"

#define RCONST 32767.

double gen_drno(char c,double *p)
{
	int i;
	double x=0.0, sum=0.0;

    switch( c )
    {
	case 'u':
	{
		/*  Generate a uniformly distributed random variate. */
		sum = ((double)(rand()))/RCONST;
		x = sum * (p[1]-p[0]) + p[0];
		break;
	}
	case 'n':
	{
        	/*  Generate a normally distributed rand var with zero mean
        	and unit variance.  */
        	for( i=0; i<12; i++ )
        	{       sum += ((double)(rand()))/RCONST;        }
        	sum -= 6.0;
		x = p[1]*sum + p[0];
		break;
	}
	default:
	{	printf("\n************************************************** \
			\n**  gen_drno() called with inappropriate type.  *** \
			\n************************************************** \
			\n");
	}
    }

	return( x );
}
