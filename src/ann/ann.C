// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

// ***************************
// ann.C - implementation file
// ***************************

// DAKOTA includes
//#include "data_types.h"
#include "surfpack.h"
#include "system_defs.h"
#include <vector>
#include <sstream>
#include "vector_enhancements.h"

// ANN includes
#include "utilities.h"
#include "ann.h"
#include "random.h"
#include "convert.h"
//extern "C" {
#include "allocate.h"
#include "dsvd2.h"
//}

using namespace std;
using namespace surfpack;


//#define TEST_PERCENT  0.05  //defined in call to set_aside_test_exemplars

// ********************************************************
// Last Modified:  6/28/99
// This function used to get the examplars to be used for training
// from a file, but is not currently used with the DAKOTA
// implementation.  Function may be used in the future.
// ********************************************************
int  ANNApprox::select_approximation_exemplars(void) 
{ return(0); }

// ********************************************************
// Last Modified:  6/28/99
// This function will normalize the data to be used for training.
// The function requires a "low" and "high" argument that
// specifies the lower and upper range to which the data will
// be normalized. e.g. (-0.8.0.8)
// Currently, the function calls for a "norm_bound" factor, like 0.8,
// to set the bounds for the normalization.  The code will used the
// negative value of "norm_bound" to set the limits.
// ********************************************************
int  ANNApprox::normalize_data(vector< vector<double> > Inputs,
			       vector< vector<double> > Outputs, double norm_bound)
{
   int j,k,in_exemplars,out_exemplars;
   double x,delta,low;
   
   normalize_factor = norm_bound;
   low = -normalize_factor;

   numInputs = num_columns(Inputs);
   numOutputs = num_columns(Outputs);
   in_exemplars = num_rows(Inputs);
   out_exemplars = num_rows(Outputs);
   if(in_exemplars == out_exemplars)
     numExemplars = in_exemplars;
   else {
     cerr << "ERROR: Row size of INPUT matrix must equal Row size of OUTPUT "
	  << "matrix" << endl;
     exit(-1);
   }

   minInputs.resize(numInputs);
   maxInputs.resize(numInputs);
   deltaInputs.resize(numInputs);
   minOutputs.resize(numOutputs);
   maxOutputs.resize(numOutputs);
   deltaOutputs.resize(numOutputs);

   maxInputs  = FindMax(Inputs);
   minInputs  = FindMin(Inputs);
   maxOutputs = FindMax(Outputs);
   minOutputs = FindMin(Outputs);

   for(k=0;k<numInputs;k++) {
      deltaInputs[k] = maxInputs[k] - minInputs[k];
   }
   for(k=0;k<numOutputs;k++) {
      deltaOutputs[k] = maxOutputs[k] - minOutputs[k];
   }
   reshape_2d(trainInputs,numExemplars,numInputs);
   reshape_2d(trainOutputs,numExemplars,numOutputs);

   //Normalizing Input and Output Matrices
   delta = normalize_factor-low;
   for(j=0;j<numExemplars;j++) {
      for(k=0;k<numInputs;k++) {
         x = Inputs[j][k];
         trainInputs[j][k] = delta*((x-minInputs[k])/deltaInputs[k])+low;
      }
      for(k=0;k<numOutputs;k++) {
         x = Outputs[j][k];
         trainOutputs[j][k] = delta*((x-minOutputs[k])/deltaOutputs[k])+low;
      }
   }
	//printf("ANNApprox::normalize_data           ...complete\n");
   return(0);
}

// ********************************************************
// Last Modified:  4/7/98
// This function will set asinde TEST_PERCENT of the specified training
// data for testing.  The function will randomly pick TEST_PERCENT of the
// training exemplars in an attempt to get a reprsentative set
// for testing.
// ********************************************************
int  ANNApprox::set_aside_test_exemplars(double TEST_PERCENT)
{
   int j,k;
   vector< vector<double> > temp_in,temp_out;
   vector<int> index;

   reshape_2d(temp_in,numExemplars,numInputs);
   reshape_2d(temp_out,numExemplars,numOutputs);
   index.resize(numExemplars);

   for(j=0;j<numExemplars;j++) {
      index[j] = j;
      for(k=0;k<numInputs;k++) {
         temp_in[j][k] = trainInputs[j][k];
      }
      for(k=0;k<numOutputs;k++) {
         temp_out[j][k] = trainOutputs[j][k];
      }
   }
   tumble_idx(numExemplars,index);

   numTestcases = (int)floor(TEST_PERCENT*numExemplars);
   numExemplars = numExemplars - numTestcases;

   reshape_2d(trainInputs,numExemplars,numInputs);
   reshape_2d(trainOutputs,numExemplars,numOutputs);
   reshape_2d(testInputs,numTestcases,numInputs);
   reshape_2d(testOutputs,numTestcases,numOutputs);

   //printf("Test Exemplars...\n");
   for(j=0;j<numTestcases;j++) {
      for(k=0;k<numInputs;k++) {
         testInputs[j][k] = temp_in[index[j]][k];
         //printf("%2d in[%2d][%d] = %9.5f  ",index[j],j,k,testInputs[j][k]);
      }
      for(k=0;k<numOutputs;k++) {
         testOutputs[j][k] = temp_out[index[j]][k];
         //printf("  out[%2d][%d] = %9.5f  \n",j,k,testOutputs[j][k]);
      }
   }
   //printf("Train Exemplars...\n");
   for(j=0;j<numExemplars;j++) {
      for(k=0;k<numInputs;k++) {
         trainInputs[j][k] = temp_in[index[j+numTestcases]][k];
         //printf("%2d in[%2d][%d] = %9.5f  ",index[j+numTestcases],j,k,trainInputs[j][k]);
      }
      for(k=0;k<numOutputs;k++) {
         trainOutputs[j][k] = temp_out[index[j+numTestcases]][k];
         //printf("  out[%2d][%d] = %9.5f  \n",j,k,trainOutputs[j][k]);
      }
   }
	//printf("ANNApprox::set_aside_test_exemplars ...complete\n");
   return(0);
}

// ********************************************************
// Last Modified:  4/28/98
// This function will train the neural network.  All variables
// necessary for implementation of the network are native to
// DAKOTA, i.e. hiddenWeights is a DakotaRealMatrix.  Those
// variables used in the SVD algorithm with temporarily be
// copied to double** and double* variables.  Once the SVD
// algorithm is complete, the data will be copied back to
// vector< vector<double> > or DakotaValArrays.
// ********************************************************
int  ANNApprox::build_approximation(double svdfactor, int Neurons)
{
   int j,k,l,index;
   double low,high;
   double **SVD_MATRIX, **U, **V, *W, sum, sumW;
   vector< vector<double> > temp1,temp2,r,winv,q,u,v,svd_matrix;
   vector< vector<double> > temp3,temp4,y,error;
   vector<double> w,ones,ave_error,std_error,rms_error;

   //Check to make sure the number of neurons specified in the 
   //approximation is NOT greater than the number of exemplars-1.
   if(Neurons > numExemplars - 1) {
     Neurons = numExemplars -1;
     printf("ANNApprox::build_approximation      ***WARNING***\n");
     printf("  Number of neurons must be less than the number of training exemplars!\n");
   }
   //Check to make sure the number of neurons are not too big.
   //If the number of neurons are too large, the convergence on the SVD
   //calculation my be limited and the program will exit.
   if(Neurons > 100) {
     Neurons = 100;
     printf("  Number of neurons being set to 100.\n");
   }
   //Reshape to proper dimensions
   reshape_2d(hiddenWeights,numInputs,Neurons);
   hiddenOffsets.resize(Neurons);
   reshape_2d(outputWeights,Neurons,numOutputs);
   outputOffsets.resize(numOutputs);
   reshape_2d(temp1,Neurons+1,Neurons+1);
   reshape_2d(temp2,Neurons+1,numExemplars);
   reshape_2d(r,numExemplars,Neurons);
   reshape_2d(winv,Neurons+1,Neurons+1);
   reshape_2d(q,Neurons+1,numOutputs);
   reshape_2d(u,numExemplars,Neurons+1);
   w.resize(Neurons+1);
   reshape_2d(v,Neurons+1,Neurons+1);
   reshape_2d(svd_matrix,numExemplars,Neurons+1);

   low  = -1.0;
   high =  1.0;

   numNeurons = Neurons;
   svd_factor = svdfactor;
   if(svd_factor <= 0.0 || svd_factor >= 1.0) {
      cout << "svd_factor is not within (0,1)...  being set to 0.999\n";
      svd_factor = 0.999;
   }
   //Currently using srand and rand in RandomReal functions.
   //Use
   RandomizeSeed(91);
   //if you use drand48 in random.C in the future.
   
   srand(91);
   for(j=0;j<numInputs;j++) {
      for(k=0;k<numNeurons;k++) {
         hiddenWeights[j][k] = RandomReal(low,high);
      }
   }
   for(j=0;j<numNeurons;j++) {
      hiddenOffsets[j] = RandomReal(low,high);
   }

   for(j=0;j<numExemplars;j++) {
      for(k=0;k<numNeurons;k++) {
         sum = 0.0;
         for(l=0;l<numInputs;l++) {
            sum += trainInputs[j][l]*hiddenWeights[l][k];
         }
         r[j][k] = sum;
      }
   }
   for(j=0;j<numExemplars;j++) {
      for(k=0;k<numNeurons;k++) {
         r[j][k] += hiddenOffsets[k];
         r[j][k] = tanh(r[j][k]);
      }
   }
   for(j=0;j<numExemplars;j++) {
      for(k=0;k<numNeurons;k++) {
         svd_matrix[j][k] = r[j][k];
      }
      svd_matrix[j][numNeurons] = 1.0;
   }

   //**********************
   DakotaToPtrMatrix(svd_matrix,&SVD_MATRIX);
   DakotaToPtrMatrix(u,&U);
   DakotaToPtrMatrix(v,&V);
   DakotaToPtrArray (w,&W);
      dsvd2(SVD_MATRIX,numExemplars,numNeurons+1,U,W,V);
   PtrMatrixToDakota(&SVD_MATRIX,svd_matrix);
   PtrMatrixToDakota(&U,u);
   PtrMatrixToDakota(&V,v);
   PtrArrayToDakota (&W,w);
   //**********************

   sumW = 0.0;
   for(j=0;j<numNeurons+1;j++) {
      sumW += w[j];
   }
   sum  = 0.0;
   index = numNeurons+1;
   //cout << "Total number of terms possible is " << numNeurons << endl;
   for(j=0;j<numNeurons+1;j++) {
      sum += w[j];
      index = j;
      //printf("w[%2d] = %3.9e \tsum/sumW = %3.9e \n",j,w[j],sum/sumW);
      if(sum/sumW > svd_factor) break;
   }
   index = index+1;
   //cout << index << " Terms Retained (" << sum/sumW << ")" << endl;

   //Zeroing Terms not retained...
   for(j=0;j<numExemplars;j++) {
      for(k=index;k<numNeurons+1;k++) {
         u[j][k] = 0.0;
      }
   }
   for(j=index;j<numNeurons+1;j++) {
      w[j] = 0.0;
   }
   for(j=0;j<numNeurons+1;j++) {
      for(k=index;k<numNeurons+1;k++) {
         v[j][k] = 0.0;
      }
   }

   for(j=0;j<numNeurons+1;j++) {
      for(k=0;k<numNeurons+1;k++) {
         winv[j][k] = 0.0;
      }
   }

   for(j=0;j<index;j++) {
      winv[j][j] = 1.0/w[j];
   }

   for(j=0;j<numNeurons+1;j++) {
      for(k=0;k<numNeurons+1;k++) {
         sum = 0.0;
         for(l=0;l<numNeurons+1;l++) {
            sum += v[j][l]*winv[l][k];
         }
         temp1[j][k] = sum;
      }
   }
   for(j=0;j<numNeurons+1;j++) {
      for(k=0;k<numExemplars;k++) {
         sum = 0.0;
         for(l=0;l<numNeurons+1;l++) {
            sum += temp1[j][l]*u[k][l];
         }
         temp2[j][k] = sum;
      }
   }
   for(j=0;j<numNeurons+1;j++) {
      for(k=0;k<numOutputs;k++) {
         q[j][k] = 0.0;
         sum = 0.0;
         for(l=0;l<numExemplars;l++) {
            sum += temp2[j][l]*atanh(trainOutputs[l][k]);
         }
         q[j][k] = sum;
      }
   }

   for(j=0;j<numNeurons;j++) {
      for(k=0;k<numOutputs;k++) {
         outputWeights[j][k] = q[j][k];
      }
   }
   for(k=0;k<numOutputs;k++) {
      outputOffsets[k] = q[numNeurons][k];
   }
   //Begin mapping for training error
   ones.resize(numExemplars);
   ones.assign(ones.size(), 1.0);
   reshape_2d(error,numExemplars,numOutputs);
   ave_error.resize(numOutputs);
   std_error.resize(numOutputs);
   rms_error.resize(numOutputs);
   reshape_2d(temp1,numExemplars,numNeurons);
   reshape_2d(temp2,numExemplars,numNeurons);
   reshape_2d(temp3,numExemplars,numOutputs);
   reshape_2d(temp4,numExemplars,numOutputs);
       reshape_2d(r,numExemplars,numNeurons);
       reshape_2d(y,numExemplars,numOutputs);

   temp1 = multiplyMM(trainInputs,hiddenWeights);
   for(j=0;j<numExemplars;j++) {
      for(k=0;k<numNeurons;k++) {
         temp2[j][k] = ones[j]*hiddenOffsets[k];
         r[j][k] = tanh(temp1[j][k] + temp2[j][k]);
      }
   }
   temp3 = multiplyMM(r,outputWeights);
   for(j=0;j<numExemplars;j++) {
      for(k=0;k<numOutputs;k++) {
         temp4[j][k] = ones[j]*outputOffsets[k];
         y[j][k] = tanh(temp3[j][k] + temp4[j][k]);
         error[j][k] = trainOutputs[j][k] - y[j][k];
      }
   }
   for(k=0;k<numOutputs;k++) {
         ave_error[k] = 0.0;
         std_error[k] = 0.0;
         rms_error[k] = 0.0;
   }
   for(j=0;j<numExemplars;j++) {
      for(k=0;k<numOutputs;k++) {
         ave_error[k] += error[j][k];
         rms_error[k] += error[j][k]*error[j][k];
      }
   }
   for(k=0;k<numOutputs;k++) {
         ave_error[k] = ave_error[k]/double(numExemplars);
         std_error[k] = sqrt((rms_error[k] - double(numExemplars)*
         		ave_error[k]*ave_error[k])/(double(numExemplars)-1.0));
         rms_error[k] = sqrt(rms_error[k]/double(numExemplars));
   }
   //printf("\nNormalized errors associated with %d training exemplars\n",numExemplars);
   for(k=0;k<numOutputs;k++) {
     //printf("Output =%2d   ave_error = %7.4f  std_error = %7.4f  rms_error = %7.4f\n",k,ave_error[k],std_error[k],rms_error[k]);
   }

   //Begin mapping for testing error
   if(numTestcases > 0) {
     ones.resize(numTestcases);
     ones.assign(ones.size(), 1.0);
     reshape_2d(error,numTestcases,numOutputs);
     reshape_2d(temp1,numTestcases,numNeurons);
     reshape_2d(temp2,numTestcases,numNeurons);
     reshape_2d(temp3,numTestcases,numOutputs);
     reshape_2d(temp4,numTestcases,numOutputs);
     reshape_2d(r,numTestcases,numNeurons);
     reshape_2d(y,numTestcases,numOutputs);

     temp1 = multiplyMM(testInputs,hiddenWeights);
     for(j=0;j<numTestcases;j++) {
       for(k=0;k<numNeurons;k++) {
	 temp2[j][k] = ones[j]*hiddenOffsets[k];
	 r[j][k] = tanh(temp1[j][k] + temp2[j][k]);
       }
     }
     temp3 = multiplyMM(r,outputWeights);
     for(j=0;j<numTestcases;j++) {
       for(k=0;k<numOutputs;k++) {
	 temp4[j][k] = ones[j]*outputOffsets[k];
	 y[j][k] = tanh(temp3[j][k] + temp4[j][k]);
	 error[j][k] = testOutputs[j][k] - y[j][k];
       }
     }
     for(k=0;k<numOutputs;k++) {
       ave_error[k] = 0.0;
       std_error[k] = 0.0;
       rms_error[k] = 0.0;
     }
     for(j=0;j<numTestcases;j++) {
       for(k=0;k<numOutputs;k++) {
	 ave_error[k] += error[j][k];
	 rms_error[k] += error[j][k]*error[j][k];
       }
     }
     for(k=0;k<numOutputs;k++) {
       ave_error[k] = ave_error[k]/double(numTestcases);
       std_error[k] = sqrt((rms_error[k] - double(numTestcases)*
		      ave_error[k]*ave_error[k])/(double(numTestcases)-1.0));
       rms_error[k] = sqrt(rms_error[k]/double(numTestcases));
     }
     //printf("\nNomralized errors associated with %d testing exemplars\n",numTestcases);
     for(k=0;k<numOutputs;k++) {
       //printf("Output =%2d   ave_error = %7.4f  std_error = %7.4f  rms_error = %7.4f\n",k,ave_error[k],std_error[k],rms_error[k]);
     }
   }
   else {
     //printf("\nNormalized errors associated with %d testing exemplars\nN/A",0);
   }
   //printf("ANNApprox::build_approximation      ...complete\n");
   return(0);
}

// ********************************************************
// Last Modified:  6/28/99
// This function will use trained parameters for the neural
// network and perform an input/output mapping for a given
// set of input values sent to the function in the ANNInput
// array.  The "mapped" output from the neural network will
// be returned in the array ANNOutput.
// ********************************************************
int  ANNApprox::map(vector<double> ANNInput, vector<double> *ANNOutput)
{
   int j,k,exemplars;
   double delta,x,low;
   vector<double> ones;
   vector< vector<double> > inputs,outputs,temp1,temp2,temp3,temp4,r,y;

   low = -normalize_factor;
   //this function could be implemented by mapping a number of data points
   //or exemplars at one time,but is implemented here assuming one data point
   //is being mapped at a time.
   exemplars = 1; //number of exemplars to be mapped  
  
   // set array and matrices to appropriate sizes
   ones.resize(exemplars);
   ones.assign(ones.size(), 1.0);
   //temp1.resize(exemplars);
   //temp2.resize(exemplars);
   //temp3.resize(exemplars);
   //temp4.resize(exemplars);
   //r.resize(exemplars);
   //y.resize(exemplars);
   //inputs.resize(exemplars);
   //outputs.resize(exemplars);

   //for (int i = 0; i < exemplars; i++) {
   //  temp1[i].resize(numNeurons);
   //  temp2[i].resize(numNeurons);
   //  temp3[i].resize(numOutputs);
   //  temp4[i].resize(numOutputs);
   //  r[i].resize(numNeurons);
   //  y[i].resize(numOutputs);
   //  inputs[i].resize(numInputs);
   //  outputs[i].resize(numOutputs);
   //}

   // new reshaping code using vector_enhancements methods instead of DAKOTA
   // data structures 
   ones.resize(exemplars);
   ones.assign(ones.size(), 1.0);
   reshape_2d(temp1,exemplars,numNeurons);
   reshape_2d(temp2,exemplars,numNeurons);
   reshape_2d(temp3,exemplars,numOutputs);
   reshape_2d(temp4,exemplars,numOutputs);
       reshape_2d(r,exemplars,numNeurons);
       reshape_2d(y,exemplars,numOutputs);
   reshape_2d(inputs,exemplars,numInputs);
   reshape_2d(outputs,exemplars,numOutputs);
   

   // old reshaping code from DAKOTA
   //ones.reshape(exemplars);
   //ones = 1.0;
   //temp1.reshape_2d(exemplars,numNeurons);
   //temp2.reshape_2d(exemplars,numNeurons);
   //temp3.reshape_2d(exemplars,numOutputs);
   //temp4.reshape_2d(exemplars,numOutputs);
   //    r.reshape_2d(exemplars,numNeurons);
   //    y.reshape_2d(exemplars,numOutputs);
   //inputs.reshape_2d(exemplars,numInputs);
   //outputs.reshape_2d(exemplars,numOutputs);





   delta = normalize_factor-low;

   for(j=0;j<exemplars;j++) {
      for(k=0;k<numInputs;k++) {
         x = ANNInput[k];
         x = delta*((x-minInputs[k])/deltaInputs[k])+low;
         inputs[j][k] = x;
      }
      //for(k=0;k<numOutputs;k++) {
         //x = ANNOutput[k];
         //x = delta*((x-minOutputs[k])/deltaOutputs[k])+low;
         //outputs[j][k] = x;
      //}
   }
   temp1 = multiplyMM(inputs,hiddenWeights);
   for(j=0;j<exemplars;j++) {
      for(k=0;k<numNeurons;k++) {
         temp2[j][k] = ones[j]*hiddenOffsets[k];
         r[j][k] = tanh(temp1[j][k] + temp2[j][k]);
      }
   }
   temp3 = multiplyMM(r,outputWeights); 
   for(j=0;j<exemplars;j++) {
      for(k=0;k<numOutputs;k++) {
         temp4[j][k] = ones[j]*outputOffsets[k];
         y[j][k] = tanh(temp3[j][k] + temp4[j][k]);
         (*ANNOutput)[k] = minOutputs[k]+((y[j][k]-low)/delta)*deltaOutputs[k];
      }
   }
//	printf("ANNApprox::map                      ...complete\n");
   return(0);
}

//********************************************************
// Last Modified:  4/7/98
int  ANNApprox::update_approximation(void)
{ return(0); }
//********************************************************

void ANNApprox::writeBinary(std::ostream& os)
{
  os.write(reinterpret_cast<char*>(&normalize_factor),sizeof(normalize_factor));
  os.write(reinterpret_cast<char*>(&numInputs),sizeof(numInputs));
  os.write(reinterpret_cast<char*>(&numNeurons),sizeof(numNeurons));
  os.write(reinterpret_cast<char*>(&numOutputs),sizeof(numOutputs));
  int i,j;
  for(i = 0; i < numInputs; i++) { 
    os.write(reinterpret_cast<char*>(&minInputs[i]),sizeof(minInputs[i]));
  }
  for(i = 0; i < numInputs; i++) { 
    os.write(reinterpret_cast<char*>(&deltaInputs[i]),sizeof(deltaInputs[i]));
  }
  for(i = 0; i < numNeurons; i++) { 
    os.write(reinterpret_cast<char*>(&hiddenOffsets[i]),sizeof(hiddenOffsets[i]));
  }
  for(i = 0; i < numInputs; i++) { 
    for(j = 0; j < numNeurons; j++) {
      os.write(reinterpret_cast<char*>(&hiddenWeights[i][j]),sizeof(hiddenWeights[i][j]));
    }
  }
  for(i = 0; i < numOutputs; i++) { 
    os.write(reinterpret_cast<char*>(&minOutputs[i]),sizeof(minOutputs[i]));
  }
  for(i = 0; i < numOutputs; i++) { 
    os.write(reinterpret_cast<char*>(&deltaOutputs[i]),sizeof(deltaOutputs[i]));
  }
  for(i = 0; i < numOutputs; i++) { 
    os.write(reinterpret_cast<char*>(&outputOffsets[i]),sizeof(outputOffsets[i]));
  }
  for(i = 0; i < numNeurons; i++) { 
    for(j = 0; j < numOutputs; j++) {
      os.write(reinterpret_cast<char*>(&outputWeights[i][j]),sizeof(outputWeights[i][j]));
    }
  }
}

void ANNApprox::writeText(std::ostream& os)
{
  // ios_base::flags return object of type ios::fmtflags but OSF compiler doesn't like it
  long old_flags = os.flags();
  unsigned old_precision = os.precision(surfpack::output_precision);
  os.setf(ios::scientific);
  
  os << normalize_factor << " normalizing factor" << endl;
  os << numInputs << " number of inputs" << endl;
  os << numNeurons << " number of neruons" << endl;
  os << numOutputs << " number of outputs" << endl;
  int i,j;
  for(i = 0; i < numInputs; i++) { 
    os << minInputs[i] << " minInputs[" << i << "]" << endl;
  }
  for(i = 0; i < numInputs; i++) { 
    os << deltaInputs[i] << " deltaInputs[" << i << "]" << endl;
  }
  for(i = 0; i < numNeurons; i++) { 
    os << hiddenOffsets[i] << " hiddenOffsets[" << i << "]" << endl;
  }
  for(i = 0; i < numInputs; i++) { 
    for(j = 0; j < numNeurons; j++) {
      os << hiddenWeights[i][j] 
         << " hiddenWeights[" << i << "][" << j << "]" << endl;
    }
  }
  for(i = 0; i < numOutputs; i++) { 
    os << minOutputs[i] << " minOutputs[" << i << "]" << endl;
  }
  for(i = 0; i < numOutputs; i++) { 
    os << deltaOutputs[i] << " deltaOutputs[" << i << "]" << endl;
  }
  for(i = 0; i < numOutputs; i++) { 
    os << outputOffsets[i] << " outputOffsets[" << i << "]" << endl;
  }
  for(i = 0; i < numNeurons; i++) { 
    for(j = 0; j < numOutputs; j++) {
      os << outputWeights[i][j] 
         << " outputWeights[" << i << "][" << j << "]" << endl;
    }
  }
  os.precision(old_precision);
  os.flags(old_flags);
}

void ANNApprox::readBinary(std::istream& is)
{
  is.read(reinterpret_cast<char*>(&normalize_factor),sizeof(normalize_factor));
  is.read(reinterpret_cast<char*>(&numInputs),sizeof(numInputs));
  is.read(reinterpret_cast<char*>(&numNeurons),sizeof(numNeurons));
  is.read(reinterpret_cast<char*>(&numOutputs),sizeof(numOutputs));
  int i,j;
  minInputs.resize(numInputs);
  for(i = 0; i < numInputs; i++) { 
    is.read(reinterpret_cast<char*>(&minInputs[i]),sizeof(minInputs[i]));
  }
  deltaInputs.resize(numInputs);
  for(i = 0; i < numInputs; i++) { 
    is.read(reinterpret_cast<char*>(&deltaInputs[i]),sizeof(deltaInputs[i]));
  }
  hiddenOffsets.resize(numNeurons);
  for(i = 0; i < numNeurons; i++) { 
    is.read(reinterpret_cast<char*>(&hiddenOffsets[i]),sizeof(hiddenOffsets[i]));
  }
  reshape_2d(hiddenWeights,numInputs,numNeurons);
  for(i = 0; i < numInputs; i++) { 
    for(j = 0; j < numNeurons; j++) {
      is.read(reinterpret_cast<char*>(&hiddenWeights[i][j]),sizeof(hiddenWeights[i][j]));
    }
  }
  minOutputs.resize(numOutputs);
  for(i = 0; i < numOutputs; i++) { 
    is.read(reinterpret_cast<char*>(&minOutputs[i]),sizeof(minOutputs[i]));
  }
  deltaOutputs.resize(numOutputs);
  for(i = 0; i < numOutputs; i++) { 
    is.read(reinterpret_cast<char*>(&deltaOutputs[i]),sizeof(deltaOutputs[i]));
  }
  outputOffsets.resize(numOutputs);
  for(i = 0; i < numOutputs; i++) { 
    is.read(reinterpret_cast<char*>(&outputOffsets[i]),sizeof(outputOffsets[i]));
  }
  reshape_2d(outputWeights,numNeurons,numOutputs);
  for(i = 0; i < numNeurons; i++) { 
    for(j = 0; j < numOutputs; j++) {
      is.read(reinterpret_cast<char*>(&outputWeights[i][j]),sizeof(outputWeights[i][j]));
    }
  }
  
}

void ANNApprox::readText(std::istream& is)
{
  string sline;
  istringstream streamline;
  getline(is,sline); streamline.str(sline);
  streamline >> normalize_factor;
  getline(is,sline); streamline.str(sline);
  streamline >> numInputs;
  getline(is,sline); streamline.str(sline);
  streamline >> numNeurons;
  getline(is,sline); streamline.str(sline);
  streamline >> numOutputs;
  int i,j;
  minInputs.resize(numInputs);
  for(i = 0; i < numInputs; i++) { 
    getline(is,sline); streamline.str(sline);
    streamline >> minInputs[i];
  }
  deltaInputs.resize(numInputs);
  for(i = 0; i < numInputs; i++) { 
    getline(is,sline); streamline.str(sline);
    streamline >> deltaInputs[i];
  }
  hiddenOffsets.resize(numNeurons);
  for(i = 0; i < numNeurons; i++) { 
    getline(is,sline); streamline.str(sline);
    streamline >> hiddenOffsets[i];
  }
  reshape_2d(hiddenWeights,numInputs,numNeurons);
  for(i = 0; i < numInputs; i++) { 
    for(j = 0; j < numNeurons; j++) {
      getline(is,sline); streamline.str(sline);
      streamline >> hiddenWeights[i][j];
    }
  }
  minOutputs.resize(numOutputs);
  for(i = 0; i < numOutputs; i++) { 
    getline(is,sline); streamline.str(sline);
    streamline >> minOutputs[i];
  }
  deltaOutputs.resize(numOutputs);
  for(i = 0; i < numOutputs; i++) { 
    getline(is,sline); streamline.str(sline);
    streamline >> deltaOutputs[i];
  }
  outputOffsets.resize(numOutputs);
  for(i = 0; i < numOutputs; i++) { 
    getline(is,sline); streamline.str(sline);
    streamline >> outputOffsets[i];
  }
  reshape_2d(outputWeights,numNeurons,numOutputs);
  for(i = 0; i < numNeurons; i++) { 
    for(j = 0; j < numOutputs; j++) {
      getline(is,sline); streamline.str(sline);
      streamline >> outputWeights[i][j];
    }
  }
}
