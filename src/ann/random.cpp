/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack_config.h"

/*
 * File: random.C
 * -----------------
 * This file provides the implementation of functions used
 * to generate random real and integer values.  The function
 * tumble_idx will also randomize the indices in a given vector. 
 */

#include <vector>

#include "random.h" 
//#include "data_types.h"
#include "system_defs.h"
#ifdef HAVE_STD
#include <cstdlib>
#else
#include <stdlib.h>
#endif

using namespace std;
/*
 * Function: gen_dscrno
 * Usage: int gen_dscrno(low,high)
 * -------------------------------
 * This function will generate a discrete random
 * integer between "low" and "high".
 */

#define eps 1.0e-10
int gen_dscrno(int low,int high)
{
   double x;

   x = drand48();
   return(low+((int)(x*(((double)(high-low))+(1.0-eps)))));
}
#undef eps

/*
 * Function: tumble_idx 
 * Usage: tumble_idx(num,vector)
 * -----------------------------
 * This function will generate a discrete random
 * integer between "low" and "high". 
 */
void tumble_idx(int num,vector<int>& vector)
{
   int i,j,temp;

   for(i=0;i<num;i++) {
      j = gen_dscrno(0,num-1);
      temp = vector[j];
      vector[j] = vector[i];
      vector[i] = temp;
   }
}

/*
 * Function: RandomizeSeed
 * -----------------------
 * This function operates by setting the random number
 * seed to the given int. 
 */
//void srand48(long int seed);
void RandomizeSeed(long int seed)
{
    srand48(seed);
}


/*
 * Function: RandomizeTime
 * -----------------------
 * This function operates by setting the random number
 * seed to the current time.  The srand function is
 * provided by the <stdlib.h> library and requires an
 * integer argument.  The time function is provided
 * by <time.h>.
 */

void RandomizeTime(void)
{
    srand48((int) time(NULL));
}

/*
 * Function: RandomInteger
 * -----------------------
 * This function first obtains a random integer in
 * the range [0..RAND_MAX] by applying four steps:
 * (1) Generate a real number between 0 and 1.
 * (2) Scale it to the appropriate range size.
 * (3) Truncate the value to an integer.
 * (4) Translate it to the appropriate starting point.
 */

int RandomInteger(int low, int high)
{
    int k;
    double d;

    d = (double) rand() / ((double) RAND_MAX + 1);
//    d = (double) lrand48() / ((double) RAND_MAX + 1);
    k = (int) (d * (high - low + 1));
    return (low + k);
}
/*
 * Function: RandomReal
 * --------------------
 * The implementation of RandomReal is similar to that
 * of RandomInteger, without the truncation step.
 */

double RandomReal(double low, double high)
{
    double d;

    d = drand48();
//    d = (double) rand() / ((double) RAND_MAX + 1);
    return (low + d * (high - low));
}

