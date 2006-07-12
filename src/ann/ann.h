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

// **************************
// ann.h - specification file
// **************************
#ifndef ann_h
#define ann_h

//#include "data_types.h"
#include <vector>

class ANNApprox
{
public:
    int   select_approximation_exemplars(void);
    int   normalize_data(std::vector< std::vector<double> > Inputs,std::vector< std::vector<double> > Outputs,double norm_bound);
    int   set_aside_test_exemplars(double TEST_PERCENT);
    int   build_approximation(double svdfactor,int Neurons);
    int   map(std::vector<double> ANNInput,std::vector<double> *ANNOutput);
    int   update_approximation(void);
    int   numExemplars;

// The following functions are add-ons to the ANN code found in the DAKOTA source listing.
// These functions are necessary for interfacing with Surfpack
  /// Write the surface in binary format
  void writeBinary(std::ostream& os); 

  /// Write the surface in text format
  void writeText(std::ostream& os); 

  /// Read the surface in binary format
  void readBinary(std::istream& is); 

  /// Read the surface in text format
  void readText(std::istream& is); 
private:

    int    numTestcases;
    int    numInputs;
    int    numOutputs;
    int    numNeurons;
    double normalize_factor;
    double svd_factor;

    std::vector< std::vector<double> >  trainInputs;
    std::vector< std::vector<double> >  trainOutputs;
    std::vector< std::vector<double> >  testInputs;
    std::vector< std::vector<double> >  testOutputs;
    std::vector< std::vector<double> >  hiddenWeights;
    std::vector<double>   hiddenOffsets;
    std::vector< std::vector<double> >  outputWeights;
    std::vector<double>   outputOffsets;

    std::vector<double>   minInputs;
    std::vector<double>   maxInputs;
    std::vector<double>   deltaInputs;
    std::vector<double>   minOutputs;
    std::vector<double>   maxOutputs;
    std::vector<double>   deltaOutputs;

};
#endif
