/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#include "surfpack.h"
#include "RBFNetSurface.h"
#include "SurfData.h"

using std::cerr;
using std::cout;
using std::endl;
using std::find;
using std::ios;
using std::istream;
using std::istringstream;
using std::ostream;
using std::random_shuffle;
using std::stack;
using std::setw;
using std::string;
using std::vector;

//_____________________________________________________________________________
// Data members 
//_____________________________________________________________________________

const string RBFNetSurface::name = "RBFNet";

//_____________________________________________________________________________
// Basis Function Helper Classes 
//_____________________________________________________________________________

BasisFunction::BasisFunction()
  : center(), weight(1.0), radii(1)
{

}

BasisFunction::BasisFunction(unsigned dims)
  : center(), weight(1.0), radii(dims)
{
  center.resize(dims);
}

void BasisFunction::resize(unsigned dims)
{
  center.resize(dims);
  radii.resize(dims);
}

double BasisFunction::evaluate(const vector<double>& x)
{
  assert(x.size() == center.xSize());
  double accumulator = 0.0;
  for (unsigned i = 0; i < center.xSize(); i++) {
    double temp = center[i] - x[i];
    //temp *= temp / radii[i];
    temp *= temp ;
    accumulator += temp;
  }
  //return exp(-accumulator);
  return sqrt(radii[0] + accumulator);
}

double BasisFunction::weightedEvaluate(const vector<double>& x)
{
  return weight*evaluate(x);
}

void BasisFunction::setCenter(vector<double>& center_)
{
  center = SurfPoint(center_);
  radii.resize(center.xSize());
}

void BasisFunction::setRadii(vector<double>& radii_)
{
  assert(radii_.size() == center.xSize());
  radii = radii_;
}

void BasisFunction::print(ostream& os)
{
  for (unsigned i = 0; i < center.xSize(); i++) {
    os << center[i] << " ";
  }
  os << " | ";
  for (unsigned i = 0; i < radii.size(); i++) {
    os << radii[i] << " ";
  }
  os << " weight: " << weight;
  os << endl;
}
//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

RBFNetSurface::RBFNetSurface(SurfData* sd) : Surface(sd)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  init();
}

RBFNetSurface::RBFNetSurface(const string filename) : Surface(0)
{
#ifdef __TESTING_MODE__
  constructCount++;
#endif
  read(filename);
}
RBFNetSurface::~RBFNetSurface()
{
#ifdef __TESTING_MODE__
  destructCount++;
#endif
  for (unsigned i = 0; i < bfs.size(); i++) {
    delete bfs[i]; bfs[i] = 0;
  }
}

void RBFNetSurface::init()
{
  // .1 is a magic number.  More intelligent algorithm for choosing radii 
  // is needed 
  radius = 0.1;
  free_param = 0.0;
  minPartitionSize = 2;
}

//_____________________________________________________________________________
// Overloaded Operators 
//_____________________________________________________________________________

//_____________________________________________________________________________
// Queries
//_____________________________________________________________________________

const string RBFNetSurface::surfaceName() const
{
  return name;
}

unsigned RBFNetSurface::minPointsRequired() const
{
  if (sd) {
    return sd->xSize();
  } else {
    cerr << "Cannot determine min points required without data" << endl;
    return INT_MAX;
  }
}

//double RBFNetSurface::evaluate(const vector<double>& x)
//{
//  //cout << "Evaluate----------------------------------" << endl;
//  double result = 0.0;
//  double temp;
//  for (unsigned i = 0; i < centers.size(); i++) {
//    temp = surfpack::euclideanDistance(x,centers[i].X());
//    //cout << "dist: " << temp;
//    // distance squared
//    temp *= temp;
//    // scale by size
//    temp /= sizes[i];
//    //cout << " center: " << i
//    //     << " weight: " << weights[i]
//    //     << " rbf: " << exp(-temp)
//    //     << endl;
//    result += weights[i]*exp(-temp);
//  }
//  //cout << "End Evaluate: " << result << "-------" << endl;
//  return result;
//}

double RBFNetSurface::evaluate(const vector<double>& x)
{
  //cout << "Evaluate----------------------------------" << endl;
  double result = free_param;
  for (unsigned i = 0; i < basis_functions.size(); i++) {
    result += basis_functions[i]->weightedEvaluate(x);
  }
  //cout << "End Evaluate: " << result << "-------" << endl;
  return result;
}

//_____________________________________________________________________________
// Commands 
//_____________________________________________________________________________

void RBFNetSurface::build(SurfData& surfData)
{
  partition(surfData, minPartitionSize);
  //generateManyOptions(surfData);
  selectModelBasisFunctions(surfData);
  return;
  centers.clear();
  for (unsigned i = 0; i < surfData.size(); i++) {
    centers.push_back(SurfPoint(surfData[i].X()));
  }
  sizes.resize(surfData.size());
  for (unsigned j = 0; j < sizes.size(); j++) {
    sizes[j] = radius; 
      
  }
  int numpts = static_cast<int>(surfData.size());
  int numcenters = static_cast<int>(centers.size());
  //cout << "numpts: " << numpts
  //     << " numcenters: " << numcenters
  //     << endl;
  
  // allocate space for the matrices
  double* hMatrix = new double[numpts*numcenters];
  double* responseVector = new double[numpts];
  double* hTransHMatrix = new double[numcenters*numcenters];
  double* resultVector = new double[numcenters];

  // populate the H and y matrices 
  for (unsigned point = 0; point < surfData.size(); point++) {
    for(int centerIndex = 0; centerIndex < numcenters; centerIndex++) 
    {
      double temp = surfpack::euclideanDistance(surfData[point].X(),
        centers[centerIndex].X());
      temp *= temp;
      temp /= sizes[centerIndex];
      temp = exp(-temp);
      hMatrix[point+centerIndex*numpts] = temp; 
    }
    responseVector[point] = surfData.getResponse(point);
  }

  //writeMatrix("H Matrix",hMatrix,numpts,numcenters,cout);
  //writeMatrix("response Vector",responseVector,numpts,1,cout);
  // The equation that we want to solve is H'Hw = H'y (solving for w)
  // We just filled H and y
  // We first need to compute H'y
  
  // increment (in memory) for the vectors x & y in alpha*A*x + beta*y
  int inc = 1; 

  // specify transpose ('T') or no transpose ('N')
  // Since we want H'y, we need to do transpose
  char trans = 'T';

  // dgemv performs operation y <- alpha*A*x + beta*y (A can be transposed)
  // We have alpha=1 beta=0 (no scaling), so y <- A'x 
  // We want H'y, so for us it's y <- H'y; we need to put our y vector in both the x & y positions
  
  // no scaling
  double alpha = 1.0;

  // don't perform the optional addition
  double beta = 0.0; 

  DGEMV_F77(trans,numpts,numcenters,alpha,hMatrix,numpts,responseVector,inc,
	    beta,resultVector,inc);
  //writeMatrix("result vector after dgemv",resultVector,numcenters,1,cout);
  
  // now responseVector holds H'y

  // now we need to compute H'H
  // dgemm computes (ignoring the scalars) C <- A*B+C, where A is m x k, B is k x n, C is m x n
  // in our case A (hMatrix) is numpts x numcenters and hTransHMatrix is numcenters*numcenters
  int m = static_cast<int>(numcenters);
  int n = static_cast<int>(numcenters);
  int k = static_cast<int>(numpts);
  char transa = 'T';
  char transb = 'N';
  DGEMM_F77(transa,transb,m,n,k,alpha,hMatrix,k,hMatrix,k,beta,hTransHMatrix,m);

  //writeMatrix("hTransHMatrix",hTransHMatrix,numcenters,numcenters,cout);
  // values must be passed by reference to Fortran, so variables must be declared for info, nrhs, trans
  int info;
  int nrhs=1;
  trans = 'N';
  //cout << "A Matrix: " << endl;
  ////writeMatrix(a,pts,coefficients.size(),cout);
  //cout << "B Vector: " << endl;
  ////writeMatrix(b,pts,1,cout);
  int lwork = numcenters * numcenters * 2;
  double* work = new double[lwork];
  DGELS_F77(trans,numcenters,numcenters,nrhs,hTransHMatrix,numcenters,
	    resultVector,numcenters,work,lwork,info);
  //writeMatrix("weights after dgels",resultVector,numcenters,1,cout);
  
  weights.resize(numcenters);
  for (int i = 0; i < numcenters; i++) {
    weights[i] = resultVector[i];
  }
  //cout << "A Matrix after: " << endl;
  ////writeMatrix(a,pts,numCoeff,cout);
  if (info < 0) {
          cerr << "dgels_ returned with an error" << endl;
  }
  delete [] work;
  delete [] resultVector;
  delete [] hMatrix;
  delete [] hTransHMatrix;
  delete [] responseVector;
}

void RBFNetSurface::buildCandidate(SurfData& surfData, 
  vector<BasisFunction*>& cand_bfs)
{
  int numpts = static_cast<int>(surfData.size());
  int numcenters = static_cast<int>(cand_bfs.size())+1;
  // debug code
  for (unsigned dd = 0; dd < cand_bfs.size(); dd++) {
    assert(cand_bfs[dd]);
  }
  
  // allocate space for the matrices
  double* hMatrix = new double[numpts*numcenters];
  double* responseVector = new double[numpts];
  double* hTransHMatrix = new double[numcenters*numcenters];
  double* resultVector = new double[numcenters];

  // populate the H and y matrices 
  for (unsigned point = 0; point < surfData.size(); point++) {
    for(int centerIndex = 0; centerIndex < numcenters-1; centerIndex++) 
    {
      hMatrix[point+centerIndex*numpts] 
        = cand_bfs[centerIndex]->evaluate(surfData[point].X()); 
    }
    hMatrix[point+(numcenters-1)*numpts] = 1.0;
    responseVector[point] = surfData.getResponse(point);
  }

  //writeMatrix("H Matrix",hMatrix,numpts,numcenters,cout);
  //writeMatrix("response Vector",responseVector,numpts,1,cout);
  // The equation that we want to solve is H'Hw = H'y (solving for w)
  // We just filled H and y
  // We first need to compute H'y
  
  // increment (in memory) for the vectors x & y in alpha*A*x + beta*y
  int inc = 1; 

  // specify transpose ('T') or no transpose ('N')
  // Since we want H'y, we need to do transpose
  char trans = 'T';

  // dgemv performs operation y <- alpha*A*x + beta*y (A can be transposed)
  // We have alpha=1 beta=0 (no scaling), so y <- A'x 
  // We want H'y, so for us it's y <- H'y; we need to put our y vector in both the x & y positions
  
  // no scaling
  double alpha = 1.0;

  // don't perform the optional addition
  double beta = 0.0; 

  DGEMV_F77(trans,numpts,numcenters,alpha,hMatrix,numpts,responseVector,inc,
	    beta,resultVector,inc);
  //writeMatrix("result vector after dgemv",resultVector,numcenters,1,cout);
  
  // now responseVector holds H'y

  // now we need to compute H'H
  // dgemm computes (ignoring the scalars) C <- A*B+C, where A is m x k, B is k x n, C is m x n
  // in our case A (hMatrix) is numpts x numcenters and hTransHMatrix is numcenters*numcenters
  int m = static_cast<int>(numcenters);
  int n = static_cast<int>(numcenters);
  int k = static_cast<int>(numpts);
  char transa = 'T';
  char transb = 'N';
  //cout << "m : " << m << endl;
  DGEMM_F77(transa,transb,m,n,k,alpha,hMatrix,k,hMatrix,k,beta,hTransHMatrix,m);

  //writeMatrix("hTransHMatrix",hTransHMatrix,numcenters,numcenters,cout);
  // values must be passed by reference to Fortran, so variables must be declared for info, nrhs, trans
  int info;
  int nrhs=1;
  trans = 'N';
  //cout << "A Matrix: " << endl;
  ////writeMatrix(a,pts,coefficients.size(),cout);
  //cout << "B Vector: " << endl;
  ////writeMatrix(b,pts,1,cout);
  int lwork = numcenters * numcenters * 2;
  double* work = new double[lwork];
  DGELS_F77(trans,numcenters,numcenters,nrhs,hTransHMatrix,numcenters,
    resultVector,numcenters,work,lwork,info);
  //writeMatrix("weights after dgels",resultVector,numcenters,1,cout);
  
  for (int i = 0; i < numcenters-1; i++) {
    //weights[i] = resultVector[i];
    cand_bfs[i]->weight = (resultVector[i]);
  }
  this->free_param = resultVector[numcenters-1];
  //cout << "free_param: " << free_param << endl;
  this->basis_functions = cand_bfs;
  //cout << "A Matrix after: " << endl;
  ////writeMatrix(a,pts,numCoeff,cout);
  if (info < 0) {
          cerr << "dgels_ returned with an error" << endl;
  }
  delete [] work;
  delete [] resultVector;
  delete [] hMatrix;
  delete [] hTransHMatrix;
  delete [] responseVector;
}

void RBFNetSurface::printSet(string header, vector<BasisFunction*>& set)
{
    cout << header << endl;
    for (unsigned k = 0; k < set.size(); k++) {
      set[k]->print(cout);
    }
}
void RBFNetSurface::selectModelBasisFunctions(SurfData& surfData)
{
  assert(surfData.size() > 0);
  int datadim = surfData[0].xSize();
  assert(partition_nodes.size() != 0);
  // Debug code begin
  //SurfData testData("linetest.txt");
  // End debug code

  vector<BasisFunction*> bestSet;
  vector<BasisFunction*> baseSet;
  vector<BasisFunction*> currentSet;
  vector<PartitionNode*> activeNodes;
  // List of active nodes initialized with root rbf
  activeNodes.push_back(partition_nodes[0]);
  bestSet.push_back(partition_nodes[0]->basis_function);
  // baseline model consists of only the root rbf
  buildCandidate(surfData,bestSet);
  builtOK = true; dataModified = false;
  double tmse = this->goodnessOfFit(string("mean_squared"),&surfData);
  double bestMetric = tmse;
  //double bestMetric = surfData.size()*log(tmse) + 
   // bestSet.size()*log(2.0*(double)surfData.size());
  //if (datadim == 1) this->getValue(testData);
  cout << "Initial metric: " << bestMetric << endl;
  // while there are active nodes, iterate through the list, and for each node
  // in the list, try adding either or both of its children.  Then, if the
  // parent of the current node is present in the model, try removing it.  
  // Try the model without the parent, but with either or both of the children.
  while (!activeNodes.empty()) {
    //cout << "Active nodes: " << activeNodes.size() << endl;
    // -1 signifies that no modification tried thus far involving the children 
    // of nodes on the active list improves the metric.
    int bestIndex = -1; 
    baseSet = bestSet;
    //printSet(string("Base Set"),baseSet);
    for (unsigned i = 0; i < activeNodes.size(); i++) {
      //cout << "Considering active node " << i << ":";
      //activeNodes[i]->basis_function->print(cout);
      PartitionNode& currentNode = *activeNodes[i];
      currentSet = baseSet;
      // add left child
      //cout << "Trying +L" << endl;
      currentSet.push_back(currentNode.left_child->basis_function);
      tryModel(surfData,i,bestIndex,bestMetric,currentSet,bestSet);
  //if (datadim == 1) this->getValue(testData);
      // add right child
      //cout << "Trying +LR" << endl;
      currentSet.push_back(currentNode.right_child->basis_function);
      tryModel(surfData,i,bestIndex,bestMetric,currentSet,bestSet);
  //if (datadim == 1) this->getValue(testData);
      // remove left child
      //cout << "Trying +R" << endl;
      vector<BasisFunction*>::iterator ritr = find(currentSet.begin(),
        currentSet.end(), currentNode.left_child->basis_function);
      assert(ritr != currentSet.end());
      currentSet.erase(ritr);
      tryModel(surfData,i,bestIndex,bestMetric,currentSet,bestSet);
  //if (datadim == 1) this->getValue(testData);
      // try four additional options if parent of current node is present
      ritr = find(currentSet.begin(),currentSet.end(),
        currentNode.basis_function);
      if (ritr != currentSet.end()) {
        // remove parent
        //cout << "Trying -P +R" << endl;
        currentSet.erase(ritr);
      tryModel(surfData,i,bestIndex,bestMetric,currentSet,bestSet);
  //if (datadim == 1) this->getValue(testData);
        // remove right child
        //cout << "Trying -P " << endl;
        ritr = find(currentSet.begin(),currentSet.end(),
          currentNode.right_child->basis_function);
 	assert(ritr != currentSet.end());
        currentSet.erase(ritr);
      tryModel(surfData,i,bestIndex,bestMetric,currentSet,bestSet);
  //if(datadim == 1) this->getValue(testData);
        // add left child
        //cout << "Trying -P +L" << endl;
        currentSet.push_back(currentNode.left_child->basis_function);
      tryModel(surfData,i,bestIndex,bestMetric,currentSet,bestSet);
  //if(datadim == 1) this->getValue(testData);
        // add right child
        //cout << "Trying -P +LR" << endl;
        currentSet.push_back(currentNode.right_child->basis_function);
      tryModel(surfData,i,bestIndex,bestMetric,currentSet,bestSet);
  //if(datadim == 1) this->getValue(testData);
      }
    }
    if (bestIndex == -1) {
      // None of the active nodes yielded an improvement, so we will remove
      // the first node in the list and replace it with its children.
      bestIndex = 0;
    } 
    assert(activeNodes[bestIndex]->right_child);
    assert(activeNodes[bestIndex]->left_child);
    //cout << "Replacing node " << bestIndex << " with its children." << endl;
    // Only put nodes on the active list if they have children
    if (activeNodes[bestIndex]->left_child->left_child) {
      assert(activeNodes[bestIndex]->left_child->right_child);
      activeNodes.push_back(activeNodes[bestIndex]->left_child);
    }
    if (activeNodes[bestIndex]->right_child->left_child) {
      assert(activeNodes[bestIndex]->right_child->right_child);
      activeNodes.push_back(activeNodes[bestIndex]->right_child);
    }
    activeNodes.erase(activeNodes.begin());
  }
  //cout << "Best: " << bestMetric << endl;
  buildCandidate(surfData,bestSet);
  //if (datadim == 1) this->getValue(testData);
  //cout << "sse: " << this->goodnessOfFit("sum_squared",&surfData) << endl;
  //cout << "mse: " << this->goodnessOfFit("mean_squared",&surfData) << endl;
  //cout << "rmae: " << this->goodnessOfFit("max_relative",&surfData) << endl;
  //cout << "rsquared: " << this->goodnessOfFit("rsquared",&surfData) << endl;

    // Begin debug code
  //SurfData dataCopy(testData);
  //dataCopy.write("considered.txt");
    // End debug code
  
}

bool RBFNetSurface::tryModel(SurfData& surfData, int currentIndex, int& bestIndex, 
    double& bestMetric, vector< BasisFunction* >& currentSet, 
    vector< BasisFunction* >& newBestSet)
{
  bool updated = false;
  //printSet("Trying set", currentSet);
  buildCandidate(surfData, currentSet);
  double tmse = this->goodnessOfFit(string("mean_squared"),&surfData);
  //cout << "Mean squared error" << tmse << endl;
  //double metric = surfData.size()*log(tmse) + 
  //  currentSet.size()*log(2.0*(double)surfData.size());
  double metric = tmse;
  //cout << "Metric: " << metric << endl;
  if (metric < bestMetric) {
    //cout << "New best metric: " << bestMetric << " RBFs: " << currentSet.size() <<endl;
    bestIndex = currentIndex;
    bestMetric = metric;
    newBestSet = currentSet;
    //cout << "New best set: " << endl;
    //for (unsigned k = 0; k < newBestSet.size(); k++) {
    //  newBestSet[k]->print(cout);
    //}
    updated = true;
  }
  //printSet("best set seen so far",newBestSet);
  return updated;
}

void RBFNetSurface::generateManyOptions(SurfData& surfData)
{
  // Debug code begin
  //SurfData testData("linetest.txt");
  // End debug code

  vector<BasisFunction*> bestBases;
  double best = 1e6;
  builtOK = true; dataModified = false;
  for (unsigned i = 0; i < 100; i++) {
    vector<BasisFunction*> basesToUse;
    //for (unsigned j = 0; j < bfs.size(); j++) {
    //  if (rand() % 3 == 0) {
    //    basesToUse.push_back(bfs[j]);
    //  }
    //}
    int numToUse = 1 + rand() % bfs.size();
    //int numToUse = 3;
    random_shuffle(bfs.begin(),bfs.end());
    for (unsigned j = 0; j < numToUse;j++) basesToUse.push_back(bfs[j]);
    buildCandidate(surfData,basesToUse);
    double tmse = this->goodnessOfFit(string("mean_squared"),&surfData);
    // Begin debug code
    //if (xsize == 1) this->getValue(testData);
    // End debug code
    double bic = surfData.size()*log(tmse) + numToUse*log(2.0*(double)surfData.size());
    if (bic < best) {
      best = bic;
      bestBases = basesToUse;
      cout << "Best: " << best << " mse: " << tmse << endl;
    }
    //cout << tmse << endl;
  }
  buildCandidate(surfData,bestBases);
  //if (xsize == 1) this->getValue(testData);
  //cout << "sse: " << this->goodnessOfFit("sse",&surfData) << endl;
  //cout << "mse: " << this->goodnessOfFit("mse",&surfData) << endl;
  //cout << "mrae: " << this->goodnessOfFit("mrae",&surfData) << endl;
  //cout << "rsquared: " << this->goodnessOfFit("rsquared",&surfData) << endl;

  // Begin debug code
  //SurfData dataCopy(testData);
  //dataCopy.write("candidateFunctions.txt");
  // End debug code
 
}

void RBFNetSurface::config(const Arg& arg)
{
  string argname = arg.name;
  if (argname == "radius") {
    radius = arg.getRVal()->getReal();
  } else if (argname == "min_partition") {
    minPartitionSize = arg.getRVal()->getInteger();    
  } else {
    Surface::config(arg);
  }
}
/// Create a surface of the same type as 'this.'  This objects data should
/// be replaced with the dataItr passed in, but all other attributes should
/// be the same (e.g., a second-order polynomial should return another 
/// second-order polynomial.  Surfaces returned by this method can be used
/// to compute the PRESS statistic.
RBFNetSurface* RBFNetSurface::makeSimilarWithNewData(SurfData* surf_data)
{
  return new RBFNetSurface(surf_data);
}

//_____________________________________________________________________________
// Helper methods 
//_____________________________________________________________________________

double RBFNetSurface::computeMetric(vector<double>& left, 
  vector<double>& right)
{
  double leftsse = surfpack::sum_squared_deviations(left);
  double rightsse = surfpack::sum_squared_deviations(right);
  double sum = leftsse + rightsse;
  //cout << " Left: " << setw(15) << leftsse
  //     << " Right: " << setw(15) << rightsse
  //     << " Sum: " << setw(15) << sum;
  return sum;
//  return surfpack::sum_squared_deviations(left) + 
//	 surfpack::sum_squared_deviations(right);
}

void RBFNetSurface::partition(SurfData& surf_data, unsigned min_partition_size)
{
  cout << "Partition: " << min_partition_size << endl;
  // We will only look at the values for one response variable in the data
  unsigned response = surf_data.getDefaultIndex();
  // The first set to look at is the set of all points in the given data set
  vector<const SurfPoint*> allpts(surf_data.size());
  for (unsigned i = 0; i < surf_data.size(); i++) {
    allpts[i] = &surf_data[i];
  }
  // stash: each set of points that will lead to a candidate basis function
  // will be stored on the stash after processing
  vector< PartitionNode* > stash;
  // Collections of points that may still be subdivided are included in the var
  // sets 
  stack< PartitionNode* > sets;
  sets.push(new PartitionNode(allpts,0));
  while (!sets.empty()) {
    vector<const SurfPoint* >& currentSet = sets.top()->set;
    int currentSetSize = currentSet.size();
    //cout << "Points in current set: " << currentSet.size() << endl;
    // bestMetricValue: the best sum-of-squared-deviations value for
    // any of the split points for this set that have been considered so far
    double bestMetricValue = 1e200;
    // The dimension for splitting that leads to the bestMetricValue
    unsigned bestSplitDim = 0;
    // The value along the bestSplitDim that leads to the bestMetricValue 
    double bestSplitValue = 0;
    // Iterate over all dimensions, looking for the best choice to split on
    for (unsigned dim = 0; dim < currentSet[0]->xSize(); dim++) {
      // Iterate over all points in the set-- paying attention only to the
      // current dim-- looking for the best value along that dimension to 
      // split on
      for (unsigned thresholdPoint = 0; thresholdPoint < currentSet.size();
	   thresholdPoint++) {
        double currentSplitValue = currentSet[thresholdPoint]->X()[dim];
	// left: collection of response values for points where the value along 
	// current dimension is less than or equal tothe currentSplitValue
        vector<double> left;
	// right: collection of response values for points where the value along
	// current dimension (dim) is greater than the currentSplitValue
        vector<double> right;
	// Iterate over all the points in the set; put the response values
	// of the points in the 'left' or 'right' collections
        for (unsigned j = 0; j < currentSet.size(); j++) {
          const SurfPoint& sp = *currentSet[j];
          if (sp[dim] <= currentSplitValue) {
            left.push_back(sp.F(response));
          } else {   
            right.push_back(sp.F(response));
          }
        }
	// Compute the mean and standard deviation for values in the 'left'
	// collection; then do the same for the one on the right.  The metric
	// value is the sum of the standard deviations for the left and right. 
        //cout << "ThreshInd: " << setw(5) << thresholdPoint << " Thresh: " << currentSplitValue;
        double metric = computeMetric(left,right);
        //cout << endl;
        //cout << "left: " << left.size()
	//     << " right: " << right.size()
	//     << " Metric for dim " << dim 
        //     << " val " << currentSplitValue
	//     << ": " << metric << endl;
	// If the result is the best seen so far for this set of data,
	// update the values for best dim, split value, and metric value
        if (metric < bestMetricValue) {
          //cout << "New best metric: " << metric << endl;
          bestSplitDim = dim;
          bestSplitValue = currentSplitValue;
	  bestMetricValue = metric;
        }
      }
    }
    // Now that the best split dimension and split value have been determined
    // we need to actually create the new set using these parameters
    vector<const SurfPoint*> newLeft;
    vector<const SurfPoint*> newRight;
//    cout << "bestSplitDim: " << bestSplitDim
//	 << " bestSplitValue: " << bestSplitValue;
    for (unsigned k = 0; k < currentSet.size(); k++) {
      if (currentSet[k]->X()[bestSplitDim] <= bestSplitValue) {
        newLeft.push_back(currentSet[k]);
      } else {   
        newRight.push_back(currentSet[k]);
      }
    }
//        cout << " left size: " << newLeft.size();
 //       cout << " right size: " << newRight.size() << endl;

    // Move the set under consideration from the stack to the stash
    stash.push_back(sets.top());
    PartitionNode* currentNode = sets.top();
    sets.pop();
    // If the new 'left' and 'right' sets are big enough that they can be split
    // again, push them onto the stack; otherwise, put them on the stash
    PartitionNode* newRightNode = new PartitionNode(newRight,currentNode);
    PartitionNode* newLeftNode = new PartitionNode(newLeft,currentNode);
    currentNode->right_child = newRightNode;
    currentNode->left_child = newLeftNode;
    if (newRight.size() > min_partition_size) {
      //cout << "New set (right): " << newRight.size() << endl;
      sets.push(newRightNode);
    } else {
      stash.push_back(newRightNode);
    }
    if (newLeft.size() > min_partition_size) {
      //cout << "New set (left): " << newLeft.size() << endl;
      sets.push(newLeftNode);
    } else {
      stash.push_back(newLeftNode);
    }
  } // while !empty
  //for (unsigned i = 0; i < stash.size(); i++) {
    //vector< const SurfPoint* >& current = stash[i];
  //  cout << "Stash " << i << endl;
    //for (unsigned j = 0; j < current.size(); j++) {
    //  cout << *current[j] << endl;
    //}
  //}
  //exit(0);
  //cout << "stash size: " << stash.size() << endl;
  //cout << "data size: " << surf_data.size() << endl;
  computeRBFCenters(stash);
}

void RBFNetSurface::computeRBFCenters(
  vector< PartitionNode* >& partitions)
{
  // Debug code begin
  //SurfData testData("linetest.txt");
  vector<BasisFunction*> unibasis(1);
  // End debug code

  SurfData centersRBF;
  SurfData sizesRBF;
  // Iterate over all the sets, creating a new basis function
  // for each one
  for (unsigned i = 0; i < partitions.size(); i++) {
    //cout << "partition: " << i 
    //     << " size: " << partitions[i]->set.size()
    //     << " left child: " << (partitions[i]->left_child ? partitions[i]->left_child->set.size() : 0)
    //     << " right child: " << (partitions[i]->right_child ? partitions[i]->right_child->set.size() : 0)
    //     << endl;
    // Identify set that will be examined
    vector< const SurfPoint*>& currentSet = partitions[i]->set;

    // Create a radial basis function with appropriate center and radius
    const SurfPoint& firstPt = *currentSet[0];
    vector<double> newcenter(firstPt.xSize());
    vector<double> newradius(firstPt.xSize());
    // Don't worry about partitions of one point
    if (currentSet.size() == 1) {
        newcenter = currentSet[0]->X();
        // If this is the node, there should be more than one point;
        // otherwise, this node should have a parent
        assert(partitions[i]->parent);
        for (unsigned d = 0; d < newradius.size(); d++) {
          newradius[d] = partitions[i]->parent->basis_function->radii[d]/2;
        }
    } else { 
        vector<double> maxes(firstPt.xSize());
        vector<double> mins(firstPt.xSize());
        // Iterate over each dimension, initializing min, max
        for (unsigned d = 0; d < firstPt.xSize(); d++) {
          maxes[d] = mins[d] = firstPt[d];
        }
        // Iterate over all points in the set, updating maxes, mins
        for (unsigned p = 1; p < currentSet.size(); p++) {
          // Iterate over each dimension, updating max,min for that dim
          const SurfPoint currentPt = *currentSet[p];
          for (unsigned d = 0; d < firstPt.xSize(); d++) {
            if (currentPt[d] > maxes[d]) maxes[d] = currentPt[d];
            else if (currentPt[d] < mins[d]) mins[d] = currentPt[d];
          }
        }

        double radiusavg = 0.0;
        for (unsigned d = 0; d < newcenter.size(); d++) {
          newcenter[d] = (maxes[d] + mins[d]) / 2;
          newradius[d] = 100.0*fabs(maxes[d] - mins[d]);
          radiusavg += newradius[d];
        }
        radiusavg /= newcenter.size();
        for (unsigned d = 0; d < newcenter.size(); d++) {
          if (fabs(newradius[d]) < .01 * radiusavg) {
            newradius[d] = radiusavg;
          } 
        }
    }
    BasisFunction* newbf = new BasisFunction();
    newbf->setCenter(newcenter);
    newbf->setRadii(newradius);
    partitions[i]->basis_function = newbf;
    bfs.push_back(newbf);
    //newbf->print(cout);

    // Begin debug code
    builtOK = true; dataModified = false;
    unibasis[0] = newbf;
    this->basis_functions = unibasis;
    //if (firstPt.xSize() == 1) this->getValue(testData);
    // End debug code
  }
    // Begin debug code
  //SurfData dataCopy(testData);
  //dataCopy.write("testResults.txt");
    // End debug code
  partition_nodes = partitions;
    // Begin debug code
  for (unsigned x = 0; x < partition_nodes.size(); x++) {
    assert(partition_nodes[x]->basis_function);
  }
    // End debug code
  
}
  
//_____________________________________________________________________________
// I/O 
//_____________________________________________________________________________

void RBFNetSurface::writeBinary(ostream& os)
{
  unsigned numcenters = centers.size();
  unsigned n_dimensions = centers[0].xSize();
  //unsigned fsize = centers[0].fSize();
  os.write(reinterpret_cast<char*>(&numcenters),sizeof(numcenters));
  os.write(reinterpret_cast<char*>(&n_dimensions),sizeof(n_dimensions));
  //os.write(reinterpret_cast<char*>(&fsize),sizeof(fsize));
  for (unsigned j = 0; j < numcenters; j++) {
    centers[j].writeBinary(os);
    os.write(reinterpret_cast<char*>(&sizes[j]),sizeof(sizes[j]));
    os.write(reinterpret_cast<char*>(&weights[j]),sizeof(weights[j]));
  }
}

void RBFNetSurface::writeText(ostream& os)
{
  ios::fmtflags old_flags = os.flags();
  unsigned old_precision = os.precision(surfpack::output_precision);
  os.setf(ios::scientific);
  os << centers.size() << " Number of centers" << endl
     << centers[0].xSize() << " Number of dimensions" << endl;
     //<< centers[0].fSize() << " Number of responses" << endl;
  for (unsigned j = 0; j < centers.size() ; j++) {
    centers[j].writeText(os);
    os << sizes[j] << " Radius of rbf #" << j << endl
       << weights[j] << " Weight of rbf #" << j << endl;
  }
  os.flags(old_flags);
  os.precision(old_precision);
}

void RBFNetSurface::readBinary(istream& is)
{
  unsigned numcenters;  
  unsigned n_dimensions;
  //unsigned fsize;
  is.read(reinterpret_cast<char*>(&numcenters),sizeof(numcenters));
  is.read(reinterpret_cast<char*>(&n_dimensions),sizeof(n_dimensions));
  //is.read(reinterpret_cast<char*>(&fsize),sizeof(fsize));
  double size;
  double weight;
  centers.clear();
  for (unsigned j = 0; j < numcenters; j++) {
    centers.push_back(SurfPoint(n_dimensions,0,is));
    is.read(reinterpret_cast<char*>(&size),sizeof(size));
    sizes.push_back(size);
    is.read(reinterpret_cast<char*>(&weight),sizeof(weight));
    weights.push_back(weight);
  }
}

void RBFNetSurface::readText(istream& is)
{
  unsigned numcenters;  
  unsigned n_dimensions;
  // read numcenters
  string sline;
  getline(is,sline);
  istringstream streamline;
  streamline.str(sline);
  streamline >> numcenters;
  //cout << "numcenters: " << numcenters << endl;
  // read n_dimensions (dimensionality of the rbf centers)
  getline(is,sline);
  streamline.str(sline);
  string buffer = streamline.str();
  streamline >> n_dimensions;
  //cout << "n_dimensions: " << n_dimensions << endl;
  //cout << "buffer: " << buffer << endl;

  // read sizes and weights
  double size;
  double weight;
  for (unsigned j = 0; j < numcenters; j++) {
    // read the location of a center in as a surfpoint
    getline(is,sline);
    centers.push_back(SurfPoint(n_dimensions,0,sline));
    //centers[j].writeText(cout);
    //cout << endl;
    // read size of rbf center;
    getline(is,sline);
    streamline.str(sline);
    buffer = streamline.str();
    streamline >> size;
    //cout << "size: " << size << endl;
    //cout << "buffer: " << buffer << endl;
    
    sizes.push_back(size);
    // read weight of rbf center;
    getline(is,sline);
    streamline.str(sline);
    buffer = streamline.str();
    streamline >> weight;
    //cout << "weight: " << weight << endl;
    //cout << "buffer: " << buffer << endl;
    weights.push_back(weight);
  }
}

//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

#ifdef __TESTING_MODE__
  int RBFNetSurface::constructCount = 0;
  int RBFNetSurface::destructCount = 0;
#endif
