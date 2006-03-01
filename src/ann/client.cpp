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

// **********************
// client.C - client file
// **********************

#include <vector>
#include "vector_enhancements.h"
//#include "data_types.h"
#include "system_defs.h"
#include "ann.h"

using namespace std;
int main ()
{
   ANNApprox test_it;
   int num_in, num_out, num_neurons, num_exemplars;
   int    j,k,exemplar,input,output;
   double svdfactor, norm_bound;
   double buffer,percent;
   char in_file[20], data_file[20];
   ifstream input_file;
   vector< vector< double > >  InputPoints;
   vector< vector< double > >  OutputPoints;
   vector<double>  AnnInput;
   vector<double>  AnnOutput;


   cout.setf(ios::fixed, ios::floatfield);
   cout.sync_with_stdio();

   cout << "Enter the name of the input file:  ";
   cin.getline(in_file,20);
   cout << endl << "You entered:  " << in_file << endl;
   input_file.open(in_file);
   input_file >> num_exemplars; cout << "num_exemplars  = " << num_exemplars << endl;
   input_file >> num_in;        cout << "num_in         = " << num_in        << endl;
   input_file >> num_out;       cout << "num_out        = " << num_out       << endl;
   input_file >> num_neurons;   cout << "num_neurons    = " << num_neurons   << endl;
   input_file >> svdfactor;     cout << "svd_factor     = " << svdfactor     << endl;
   input_file >> norm_bound;    cout << "norm_bound     = " << norm_bound    << endl;
   input_file >> percent;       cout << "percent        = " << percent       << endl;
   input_file >> data_file;     cout << "data_file      = " << data_file     << endl;

   //**************************** From Input File **************
   input_file.close();
   input_file.open(data_file);
   reshape_2d(InputPoints,num_exemplars,num_in);
   reshape_2d(OutputPoints,num_exemplars,num_out);

   for(exemplar=0;exemplar<num_exemplars;exemplar++) {
      for(input=0;input<num_in;input++) {
         input_file >> buffer;
         InputPoints[exemplar][input] = 0.0;
         InputPoints[exemplar][input] = buffer;
         //cout << InputPoints[exemplar][input] << "   ";
      }
      for(output=0;output<num_out;output++) {
         input_file >> buffer;
         OutputPoints[exemplar][output] = 0.0;
         OutputPoints[exemplar][output] = buffer;
         //cout << OutputPoints[exemplar][output] << endl;
      }
   }
   input_file.close();
   //**************************** From Input File **************

   //test_it.select_approximation_exemplars();
   test_it.normalize_data(InputPoints,OutputPoints,norm_bound);
   test_it.set_aside_test_exemplars(percent);
   test_it.build_approximation(svdfactor,test_it.numExemplars-1);

   //**********************
   ofstream output_file;
   output_file.open("quad_output.m");
   output_file << "ann=[" << endl;
   AnnInput.resize(num_in);
   AnnOutput.resize(num_out);
   for(j=0;j<num_exemplars;j++) {
      for(k=0;k<num_in;k++) {
         AnnInput[k] = InputPoints[j][k];
	 output_file << AnnInput[k] << "   ";
      }
      test_it.map(AnnInput,&AnnOutput);
      for(k=0;k<num_out;k++) {
         output_file << AnnOutput[k] << "   ";
         output_file << OutputPoints[j][k] << "  ";
      }
      output_file << endl;
   }
   output_file << "];" << endl;
   output_file.close();
   //**********************
   
   return 0;
}

