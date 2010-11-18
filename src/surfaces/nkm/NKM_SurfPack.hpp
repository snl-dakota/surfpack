#ifndef __SURFPACK_HPP__
#define __SURFPACK_HPP__

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif

#include "NKM_SurfMat.hpp"
//not sure which of these includes and usings are necessary achieved
#include <cstdlib>
#include <fstream>
#include <iomanip> 
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cstring>
#include <csignal>
#include <set>
#include <string>
#include <vector>

namespace nkm {

int if_close_enough(double a, double b);

int nchoosek(int n, int k);

int num_multi_dim_poly_coef(int nvarsr, int ndeg);

MtxInt& multi_dim_poly_power(MtxInt& poly, int nvarsr, int ndeg, int istart=0, int jstart=0, int iffirst=1);

MtxInt& main_effects_poly_power(MtxInt& poly, int nvarsr, int ndeg);

MtxDbl& evaluate_poly_basis(MtxDbl& g, const MtxInt& poly, const MtxDbl& xr);

MtxDbl& evaluate_poly_der_basis(MtxDbl& dg, const MtxInt& poly, const MtxInt& der, const MtxDbl& xr);

MtxDbl& evaluate_poly_der(MtxDbl& result, const MtxInt& poly, const MtxInt& der, const MtxDbl& coef, const MtxDbl& xr);


MtxDbl& gen_rot_mat(MtxDbl& Rot, const MtxDbl& EulAng, int nvarsr);

MtxDbl& gen_rand_rot_mat(MtxDbl& rot,int nvarsr);

MtxDbl& gen_rand_axis_bin_opt_samples_0to1(MtxDbl& xr, int nvarsr); 

//xr is for XReal, nvarsr is for Number of VARiables Real
inline MtxDbl& rotate_xr(MtxDbl& xr_rot, const MtxDbl& rot_or_eul_ang, const MtxDbl& xr) {
    int nvarsr=xr.getNCols();
    if((rot_or_eul_ang.getNRows()==nvarsr)&&
       (rot_or_eul_ang.getNCols()==nvarsr)) {
      //rot_or_eul_ang is a rotation matrix
      matrix_mult(xr_rot,xr,rot_or_eul_ang,0.0,1.0);
    }
    else if((rot_or_eul_ang.getNRows()==((nvarsr*(nvarsr-1))/2))&&
	    (rot_or_eul_ang.getNCols()==1)) {
      //rot_or_eul_ang is a vector containing the Euler Angles for the rotation matrix
      MtxDbl rot; 
      gen_rot_mat(rot,rot_or_eul_ang,nvarsr);
      matrix_mult(xr_rot,xr,rot,0.0,1.0);
    }
    else{
      printf("Error in rotate_xr(MtxDbl& xr_rot,MtxDbl& rot_or_eul_ang,MtxDbl& xr): rot_or_eul_ang has the wrong size!!!\n");
      assert(0);
    }
    return xr_rot;
}


template<typename T>
std::vector<T>& toVec(std::vector<T>& result, const std::string& s)
{
  std::istringstream is(s);
  result.clear();
  if (s == "") 
    return result;
  T temp;
  do {
    is >> temp;
    result.push_back(temp);
  } while (!is.eof());
  
  return result;
}


template<typename T>
std::string toString(const T arg)
{
  std::ostringstream os;
  os << arg;
  return os.str();
}

template<typename T>
std::string fromVec(const std::vector<T>& vec)
{
  std::ostringstream os;
  for (typename std::vector<T>::const_iterator itr = vec.begin();
	itr != vec.end(); ++itr) {
    if (itr != vec.begin()) os << " ";
    os << *itr;
  }
  return os.str();
}

namespace surfpack {
  ///should have variable # of sig figs (precision) control for output

  /// Precision of output for double precision numbers
  const unsigned output_precision = 3;

  /// Length of the field for double-precision number stream output
  const unsigned field_width = output_precision + 8;


  /// Thrown when an attempt to open a file for reading or writing fails
  class file_open_failure: public std::runtime_error
  {
  public:
    file_open_failure(const std::string& filename = "") 
      : std::runtime_error("File " + filename + " could not be opened.") {}
  };
    
  /// Thrown when end-of-file is reached unexpectedly, when an unrecognized or
  /// unacceptable file extension is encountered, or when a file contains
  /// unexpected or illegally formatted contents.
  class io_exception: public std::runtime_error
  {
  public:
    io_exception(const std::string& msg = "") : std::runtime_error(msg) {}
  };

  /// Throw an exception if end-of-file has been reached 
  void checkForEOF(std::istream& is);

  /// specified by parameter extension
  bool hasExtension(const std::string& filename, const std::string extension);
}

} // end namespace nkm

#endif
