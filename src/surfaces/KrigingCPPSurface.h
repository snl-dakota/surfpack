#include "surfpack_config.h"

#ifndef __KRIGING_CPP_SURFACE_H__
#define __KRIGING_CPP_SURFACE_H__

#include "SurfpackParser.h"
#include "SurfpackMatrix.h"

/// Validator class for KrigingSurface
class KrigingCPPSurface : public Surface
{

//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

public:
  KrigingCPPSurface(SurfData* sd);
  KrigingCPPSurface(const std::string filename);
  ~KrigingCPPSurface(); 
private:
  /// Explicitly disallow default constructor  
  KrigingCPPSurface();
//_____________________________________________________________________________
// Overloaded operators 
//_____________________________________________________________________________


//_____________________________________________________________________________
// Queries 
//_____________________________________________________________________________
public:
  
  virtual const std::string surfaceName() const;
  static unsigned minPointsRequired(unsigned hypothetical_xsize);
  virtual unsigned minPointsRequired() const;
  virtual double evaluate(const std::vector<double>& x);

//_____________________________________________________________________________
// Commands 
//_____________________________________________________________________________

  void setConminThetaVars(const std::vector<double>& vals);
  void useUniformCorrelationValue(double correlation);
  void usePreComputedCorrelationVector(const std::vector<double>& vals);
  void build(SurfData& data);
  virtual void config(const SurfpackParser::Arg& arg);
  
  /// Create a surface of the same type as 'this.'  This objects data should
  /// be replaced with the dataItr passed in, but all other attributes should
  /// be the same (e.g., a second-order polynomial should return another 
  /// second-order polynomial.  Surfaces returned by this method can be used
  /// to compute the PRESS statistic.
  virtual KrigingCPPSurface* makeSimilarWithNewData(SurfData* surfData);

//_____________________________________________________________________________
// Helper methods 
//_____________________________________________________________________________

  void buildModel(SurfData& data);
  double correlation_function(const std::vector<double>& correlations,
    const std::vector<double>& pt1, const std::vector<double>& pt2);
  std::vector<double> useConminToFindCorrelationParams();
  double likelihoodEstimation();

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
  bool runConminFlag;
  std::vector<double> rhs;
  std::vector<double> correlationVector;
  std::vector<double> conminCorrelations;
  double likelihood;
  double betaHat;
  double determinantCorrMatrix;
//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

#ifdef __TESTING_MODE__ 
  friend class KrigingCPPUnitTest;
  friend class SurfaceUnitTest;
#endif
};

#endif
