#include "surfpack_config.h"

/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

/*
 * File: random.h
 * -----------------
 * This file provides the implementation of functions used
 * to generate random real and integer values.  Some of the
 * following functions were taken from CS106 random library
 * offered by Roberts at Stanford. 
 * The function tumble_idx will also randomize the indices
 * in a given vector. 
 */

#ifndef _random_h
#define _random_h

//#include "data_types.h"

/*
 * Function: gen_dscrno
 * Usage: int gen_dscrno(low,high)
 * ------------------------------------
 * This function will generate a discrete random
 * integer between "low" and "high". 
 */
int gen_dscrno(int low,int high);

/*
 * Function: tumble_idx 
 * Usage: tumble_idx(num,vector)
 * -----------------------------
 * This function will generate a discrete random
 * integer between "low" and "high". 
 */
void tumble_idx(int num,std::vector<int>& vector);

/*
 * Constant: RAND_MAX
 * ------------------
 * Unfortunately, several libraries that supposedly conform to
 * the ANSI standard do not define RAND_MAX in <stdlib.h>.  To
 * reduce portability problems, this interface defines RAND_MAX
 * to be the largest positive integer if it is undefined.
 */

#ifndef RAND_MAX
#  define RAND_MAX ((int) ((unsigned) ~0 >> 1))
#endif

/*
 * Function: RandomizeSeed
 * -----------------------
 * This function operates by setting the random number
 * seed to the given int.
 */
void RandomizeSeed(long int seed);


/*
 * Function: RandomizeTime
 * -----------------------
 * This function operates by setting the random number
 * seed to the current time.  The srand function is
 * provided by the <stdlib.h> library and requires an
 * integer argument.  The time function is provided
 * by <time.h>.
 */
void RandomizeTime(void);

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
int RandomInteger(int low, int high);

/*
 * Function: RandomReal
 * --------------------
 * The implementation of RandomReal is similar to that
 * of RandomInteger, without the truncation step.
 */
double RandomReal(double low, double high);

#endif
