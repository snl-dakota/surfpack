/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef SURFPACK_SYSTEM_HEADERS_H
#define SURFPACK_SYSTEM_HEADERS_H

#include <algorithm>
#include <limits>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <climits>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip> 
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <sys/time.h>
#include <vector>

typedef std::vector<double> VecDbl;
typedef std::vector<double>::const_iterator VecDblIt;
typedef std::vector<unsigned> VecUns;
typedef std::vector<unsigned>::const_iterator VecUnsIt;
typedef std::vector< std::vector< unsigned > > VecVecUns;
typedef std::vector< std::vector< double > > VecVecDbl;
typedef std::vector< std::string > VecStr;
typedef std::set< unsigned > SetUns;

typedef std::pair< std::string, std::string > ModelParam;
typedef std::map< std::string, std::string> ParamMap;
typedef std::pair< std::string, ParamMap> Command;

#endif // SURFPACK_SYSTEM_HEADERS_H

