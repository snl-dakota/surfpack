#ifndef __RBF_NET_SURFACE_H__ 
#define __RBF_NET_SURFACE_H__ 
#include "surfpack_config.h"
#include "surfpack_system_headers.h"
#include "SurfPoint.h"
#include "Surface.h"

typedef float real;

//_____________________________________________________________________________
// Basis Function Helper Classes 
//_____________________________________________________________________________

class BasisFunction
{
  public:
    BasisFunction();
    BasisFunction(unsigned dims);
    void resize(unsigned dims);
    double evaluate(const std::vector<double>& x);
    double weightedEvaluate(const std::vector<double>& x);
    SurfPoint center;
    double weight;
    std::vector<double> radii;
    void setCenter(std::vector<double>& radii_);
    void setRadii(std::vector<double>& radii_);
    void print(std::ostream& os);
};
     

//_____________________________________________________________________________
// Partition Node Helper Class 
//_____________________________________________________________________________

class PartitionNode
{
  public:
    std::vector< const SurfPoint* > set;
    BasisFunction* basis_function;
    PartitionNode* parent;
    PartitionNode* left_child;
    PartitionNode* right_child;
    PartitionNode(std::vector< const SurfPoint* >& set_, 
      PartitionNode* parent_) : set(set_),
      parent(parent_),left_child(0),right_child(0), basis_function(0) {}
    ~PartitionNode() {}
};




class RBFNetSurface : public Surface
{
//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

public:
  RBFNetSurface(SurfData* sd);
  RBFNetSurface(const std::string filename);
  ~RBFNetSurface();
  void init();

//_____________________________________________________________________________
// Overloaded Operators 
//_____________________________________________________________________________

//_____________________________________________________________________________
// Queries
//_____________________________________________________________________________

  virtual const std::string surfaceName() const;
  
  virtual unsigned minPointsRequired() const;
  
  virtual double evaluate(const std::vector<double>& x);
//_____________________________________________________________________________
// Commands 
//_____________________________________________________________________________

  virtual void build(SurfData& data);
  virtual void buildCandidate(SurfData& data, std::vector<BasisFunction*>& cbfs);
  virtual void generateManyOptions(SurfData& surfData);
  virtual void selectModelBasisFunctions(SurfData& surfData);
  
  virtual void config(const Arg& arg);
  /// Create a surface of the same type as 'this.'  This objects data should
  /// be replaced with the dataItr passed in, but all other attributes should
  /// be the same (e.g., a second-order polynomial should return another 
  /// second-order polynomial.  Surfaces returned by this method can be used
  /// to compute the PRESS statistic.
  virtual RBFNetSurface* makeSimilarWithNewData(SurfData* surf_data);

//_____________________________________________________________________________
// Helper methods 
//_____________________________________________________________________________

  virtual double computeMetric(std::vector<double>& left, 
    std::vector<double>& right);

  virtual void partition(SurfData& sd);

  virtual void computeRBFCenters(
    std::vector< PartitionNode* >& partitions);

  virtual bool tryModel(SurfData& surfData, int currentIndex, int& bestIndex, 
    double& bestMetric, std::vector< BasisFunction* >& currentSet, 
    std::vector< BasisFunction* >& newBestSet);
  virtual void printSet(std::string header, std::vector<BasisFunction*>& set);

//_____________________________________________________________________________
// I/O 
//_____________________________________________________________________________

  virtual void writeBinary(std::ostream& os);
  virtual void writeText(std::ostream& os);
  virtual void readBinary(std::istream& is);
  virtual void readText(std::istream& is);
//_____________________________________________________________________________
// Data members 
//_____________________________________________________________________________
protected:
  static const std::string name;
  std::vector<double> weights;
  std::vector<double> sizes;
  std::vector<SurfPoint> centers;
  double radius;
  std::vector<BasisFunction*> bfs;
  double free_param;
  std::vector<BasisFunction*> basis_functions;
  std::vector<PartitionNode*> partition_nodes;

//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

#ifdef __TESTING_MODE__
public:
  static int constructCount;
  static int destructCount;
#endif

};

#endif
