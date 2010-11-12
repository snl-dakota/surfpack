#ifndef __SURFDATA_HPP__
#define __SURFDATA_HPP__
#include <stdlib.h>
#include <stdio.h>
#include "NKM_SurfMat.hpp"
#include <math.h>
//#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

//#define __SURFDATA_ERR_CHECK__

namespace nkm {

using std::istream;
using std::ios;
using std::ostream;
using std::istringstream;
using std::string;

class SurfPackModel;
class SurfData;
class SurfDataScaler;


class SurfData
{
private:
  ///the number of data points
  int npts;

  ///the number of real input variables
  int nvarsr;

  ///the number of integer input variables
  int nvarsi;

  ///the number of outputs
  int nout;
  
  ///indicates the output that you are building the emulator for
  int jout;

  bool ifHaveMinMaxXr;

  MtxDbl minMaxXr;

  ///indicates groups of locked aspect ratio (group scaling) real input variables.  It is an empty matrix if you don't want to lock aspect ratios during scaling.  When it is not an empty matrix it must be a matrix with 2 rows and nvarsr columns, the top row holds the group number the bottom row holds index to the real input variable dimension (0,1,2,...,nvarsr-1), this matrix should be column sorted in the constructor before calling the scaleToDefault function, so that the dimensions with the same group will be listed sequentially.  lockx should only need to be accessed by constructors and the scaleToDefault() and groupScale() functions and the assignment operator.
  MtxInt lockxr;

  ///an npts by nout matrix, haved1y(ipt,j)=0 OR 1 => we don't OR do have the d1y for y(ipt,j) 
  MtxInt haved1y; 

  ///an npts by nout matrix, haved2y(ipt,j)=0 OR 1 => we don't OR do have the d2y for y(ipt,j) 
  MtxInt haved2y;

  ///Nd1y=nchoosek(1-1+nvarsr,1)=nvarsr =# of mixed partial first derivatives of 1 output with respect to all real inputs, I use this instead of nvarsr for clarity and consistency with the second derivative case
  int Nd1y; 

  //Nd2y=nchoosek(2-1+nvarsr,2)=# of mixed partial second derivatives of 1 output with respect to all real inputs (i.e. the number of unique terms in a hessiaon)

  ///an Nd1y by nvarsr matrix that contains the individual dimension's orders of mixed partial derivatives of the first derivative (gradient) produced by calling multi_dim_poly_power(d1ydeg,nvarsr,-1);  
  MtxInt d1ydeg;

  ///an nchoosek(2-1+nvarsr,2) by nvarsr matrix that contains the individual dimension's orders of mixed partial derivatives of total order exactly equal to 2 (i.e. the components of thee "hessian") produced by calling multi_dim_poly_power(d2ydeg,nvarsr,-2);  
  MtxInt d2ydeg;

  ///an npts by (Nd1y*nout) matrix that contains the partial 1st derivatives ("gradients") of outputs for individual points d1y(ipt,jder+jy*Nd1y)=d(y(x(ipt,:)))/d(x(*,jder)) i.e. first derivative of output variable jy with respect to input variable jder evaluated at point "ipt"
  MtxDbl d1y;

  ///an npts by (Nd2y*nout) matrix that contains the mixed partial 2nd derivatives (components of the hessian) of outputs for individual points d2y(ipt,jdeg+jy*Nd2y) = d^2y(x(ipt,:))/product(dx(*,ivarr)^d2ydeg(jdeg,ivarr);ivarr=0->nvarsr-1))
  MtxDbl d2y;


  ///Says how to unscale the real input variables: that us the unscaledxr(i,j)=xr(i,j)*fabs(unscalexr(0,j))+unscalexr(1,j).  unscalexr is a protected rather than private variable to make it easy to access the length scaling (without also shifting).  An unscalexr(0,j)<0 indicates that real input dimension j is "singular", i.e that there is only one unique value in dimension j, which can be used as a flag to ignore this dimension when constructing an emulator.  For singular dimension j, scaleToDefault() sets unscalexr(0,j) to be -1.0 and unscalexr(1,j) to be that dimensions single value, i.e. the scaled xr(:,j) will be all 0.0's since the default domain is either a hypercube or hyper-rectangle (if the aspect ratio of a group of dimensions is locked) that has volume 1 AND IS CENTERED AT 0.0
  MtxDbl unscalexr;

  ///how to unscale the output variables, these are protected instead of private to make it easy to access the length scaling (without also shifting) and unscaley(0,j)<0 indicates that output dimension j is singular (only one unique value in dimension j) which can be used as a flag to NOT construct the emulator and just return the single value for this singular output, for the "default" scaling unscaley(0,j) will be -1.0 and unscaley(1,j) will be the single value for the singular dimension i.e. the scaled y(:,j) will be all 0.0's
  MtxDbl unscaley;
  
  ///labels for the real input variables
  std::vector< std::string > xrLabels;

  ///labels for the integer input variables
  std::vector< std::string > xiLabels;

  ///labels for the output variables
  std::vector< std::string > yLabels;

  /// this just allocates space for unscalexr and unscaley and initializes them to the identity to "unscale" (it xr==unscale(xr) is true, and y==unscale(y) is true), it is a necessary feature and is named "dontScale" forseeing the likely possibility that a constructor that doesn't automatically scale will be added
  inline void dontScale(){
    //void dontScale(){
    assert((unscalexr.getNElems()==0)&&(unscaley.getNElems()==0));
    unscalexr.newSize(2,nvarsr); 
    unscaley.newSize(2,nout);
    int j;
    for(j=0;j<nvarsr;j++){
      unscalexr(0,j)=1.0;
      unscalexr(1,j)=0.0;
    }
    for(j=0;j<nout;j++){
      unscaley(0,j)=1.0;
      unscaley(1,j)=0.0;
    }
    return;
  };

  ///individually scale each dimension to length 1 centered at zero, (a "singular" dimension, i.e. a dimension in which all points share the same value, receives unscalea(0)=-1.0 and unscalea(1) is set to the single constant value, a negative unscalea(0) is a flag that says to exclude this dimension from the analysis) at the time of this writing this function is intended to only be called by the scaleToDefault() function
  void indivScale(MtxDbl& a, MtxDbl& unscalea, const MtxDbl& minmaxa, bool have_minmaxa);

  ///scale _this_ group of dimensions to a hyperrectangle of volume 1 centered at zero while preserving/locking the aspect ratio of the group of dimensions (note that this should only be called for __a__ __group__ (as in 1 group at a time) of real input variables, (a "singular" dimension j, i.e. a dimension in which all points share the same value, is not counted as part of the group, instead it is mapped to 0.0, has unscalea(0,j)=-1.0 and unscalea(1,j) is set to the single constant value) at the time of this writing this function is intended to only be called by the scaleToDefault() function
  void groupScale(MtxDbl& a, MtxDbl& unscalea, const MtxDbl& minmaxa, bool have_minmaxa);

  ///returns 1 if the jth output is singular (and sets singular_y to the single y value) and returns 0 otherwise
  int isYSingular(int j,double& singular_y) {
    singular_y=unscaley(1,j);
    return (unscaley(0,j)==-1.0);
  };

  ///returns 0 if the data is scaled or 1 if the data is unscaled
  int isUnScaled(){
    for(int j=0;j<nvarsr;++j) 
      if(!(((unscalexr(0,j)==1.0)&&(unscalexr(1,j)==0.0))
	   ||(unscalexr(0,j)==-1.0)))
	return 0;
    for(int j=0;j<nout;++j)
      if(!(((unscaley(0,j)==1.0)&&(unscaley(1,j)==0.0))
	   ||(unscaley(0,j)==-1.0)))
	return 0;
    return 1;
  };
  //debating about having these getUnscaleXr and getUnscaleY functions, if they change unscalexr or unscaley then that could be VERY bad, but calling these 2 functions is also very inconvienient because it requires the declaration of a Matrix to hold them

  ///tell me how to unscale the REAL input variables: XRUNSCALED(i,j)=xr(i,j)*fabs(unscalexr(0,j))+unscalexr(1,j) if unscalexr(0,j)<0.0 then all the points in REAL input dimension j share the same value which is stored in unscalexr(1,j), I refer to this as dimension j being "singular" and "singular" dimensions should be excluded from the analysis, unscalexr(0,j) is a "flag" to indicate to exclude dimension j
  inline MtxDbl& getUnscaleXr(MtxDbl& unscale) const {unscale.copy(unscalexr); return unscale;};
  inline MtxDbl& getUnscaleY(MtxDbl& unscale) const {unscale.copy(unscaley); return unscale;};

  ///retrieve an unscaled copy of an output dimension, by default this is the output dimension we want to build the emulator for
  MtxDbl& getYUnScaled(MtxDbl& Y, int jout_want=-99999){
    if(jout_want==-99999) jout_want=jout;    
    assert((0<=jout_want)&&(jout_want<nout));
    double unscale=fabs(unscaley(0,jout_want));
    double unshift=unscaley(1,jout_want);
    if((unscale==1.0)&&(unshift==0.0))
      return (y.getCols(Y,jout_want));

    Y.newSize(npts); 
    for(int i=0; i<npts; i++) 
      Y(i)=y(i,jout_want)*unscale+unshift; 
    return Y;
  };

  ///retrieve an unscaled copy of the real input dimensions, if you want the scaled version you can access them directly with a .xr
  MtxDbl& getXrUnScaled(MtxDbl& XR) {
    XR=xr;
    double unscale, unshift;
    int i;
    for(int j=0;j<nvarsr;j++){
      unscale=fabs(unscalexr(0,j));
      unshift=unscalexr(1,j);
      if(!((unscale==1.0)&&(unshift==0.0)))
	for(i=0; i<npts; i++)
	  XR(i,j)=XR(i,j)*unscale+unshift;
    }
    return XR;
  };

  ///scale the data (output variables and real input variables) Note that this function should only be called by the Model if it is appropriate for the model to do so, it should not be called by anything other than a model. Also note that copies of a SurfData object will use the same default scaling of the original SurfData object, i.e if during construction you told the original SurfData object to use group scaling so will the copy
  void scaleToDefault();

  ///scale real input variables to domain_new(0,j)<=xr(:,j)<=domain_new(1,j) for j=0,1,...,nvarsr;  For a "singular" dimension j xr(:,j)=0.5*(domain_new(0,j)+domain_new(1,j)),unscalexr(0,j)=-1.0 and unscalexr(1,j) = the single value - xr(:,j)
  void scaleXrToDomain(MtxDbl& domain_new);
  
  ///scale the real input variables xr to have an unscale factor of unscale_xr (scale xr so that it is unscaled by unscale_xr and reset unscalexr to unscale_xr) this will change a unscalexr(0,j)<0.0 to a -fabs(unscale_xr(0,j)) 
  void scaleXrToFactor(MtxDbl& unscale_xr);

  ///scale the output variables y to have an unscale factor of unscale_y (scale y so that it is unscaled by unscale_y and reset unscaley to unscale_y) this will change a unscaley(0,j)<0.0 to a -fabs(unscale_y(0,j))   
  void scaleYToFactor(MtxDbl& unscale_y);
  
  ///scale both the output and real input variables to be unscaled by the specified unscale factors in a single function call
  inline void scaleToFactors(MtxDbl& unscale_xr, MtxDbl& unscale_y){
    scaleXrToFactor(unscale_xr);
    scaleYToFactor(unscale_y);
    return;
  };

  ///copy the data to "result" and then unscale "result", return "result" as a reference to a SurfData object, typically only used for displaying output to user, but retain all information needed to scale it if the user later asks us to
  SurfData& unScaleCopy(SurfData& result){
    result=*this; //copy
    return (result.unScale());
  };

  ///copy the data to "result" and then unscale "result", return "result" as a a SurfData object, typically only used for displaying output to user, but retain all information needed to scale it if the user later asks us to
  SurfData unScaleCopy(){
    SurfData result(*this); //copy constructor
    return (result.unScale());
  };

  ///unscale this SurfData object but retain all information needed to scale it if the user later asks us to
  SurfData& unScale();

  MtxDbl& scaleXrOther(MtxDbl& xr_other){
#ifdef __SURFDATA_ERR_CHECK__
    assert(xr_other.getNCols()==nvarsr);
#endif 
    int npts_other=xr_other.getNRows();
    double temp_scale_mult, temp_scale_offset;
    for(int j=0; j<nvarsr; ++j) {
      temp_scale_mult=1.0/unscalexr(0,j);
      temp_scale_offset=unscalexr(1,j);
      for(int i=0; i<npts_other; ++i)
	xr_other(i,j)=(xr_other(i,j)-temp_scale_offset)*temp_scale_mult;
    }
    return xr_other;
  };

  MtxDbl& unScaleXrOther(MtxDbl& xr_other){
#ifdef __SURFDATA_ERR_CHECK__
    assert(xr_other.getNCols()==nvarsr);
#endif 
    int npts_other=xr_other.getNRows();
    double temp_scale_mult, temp_scale_offset;
    for(int j=0; j<nvarsr; ++j) {
      temp_scale_mult=unscalexr(0,j);
      temp_scale_offset=unscalexr(1,j);
      for(int i=0; i<npts_other; ++i)
	xr_other(i,j)=xr_other(i,j)*temp_scale_mult+temp_scale_offset;
    }
    return xr_other;
  };

  inline double scaleYOther(double y_other, int j=-99999){
    if(j==-99999) j=jout;
#ifdef __SURFDATA_ERR_CHECK__
    assert((0<=j)&&(j<nout));
#endif 
    return ((y_other-unscaley(1,j))/fabs(unscaley(0,j)));
  };

  inline double unScaleYOther(double y_other, int j=-99999) {
    if(j==-99999) j=jout;
#ifdef __SURFDATA_ERR_CHECK__
    assert((0<=j)&&(j<nout));
#endif 
    //printf("yunscale(:,%d)=[%g;%g]\n",j,unscaley(0,j),unscaley(1,j));
    return (y_other*fabs(unscaley(0,j))+unscaley(1,j));
  };


  MtxDbl& scaleYOther(MtxDbl& y_other, int j=-99999);

  MtxDbl& unScaleYOther(MtxDbl& y_other, int j=-99999);

  ///for first derivative of y(:,j_y) wrt xr(:,j_xr)
  inline double scaleFactorDerY(int j_xr, int j_y=-99999) {
    if(j_y==-99999) j_y=jout;
#ifdef __SURFDATA_ERR_CHECK__
    assert((0<=j_xr)&&(j_xr<nvarsr)&&(0<=j_y)&&(j_y<nout));
#endif 
    return(fabs(unscalexr(0,j_xr)/unscaley(0,j_y)));
  };

  ///for first derivative of y(:,j_y) wrt xr(:,j_xr)
  inline double unScaleFactorDerY(int j_xr, int j_y=-99999) {
    if(j_y==-99999) j_y=jout;
#ifdef __SURFDATA_ERR_CHECK__
    assert((0<=j_xr)&&(j_xr<nvarsr)&&(0<=j_y)&&(j_y<nout));
#endif 
    return(fabs(unscaley(0,j_y)/unscalexr(0,j_xr)));
  };
    
  ///for arbitrary order mixed partial derivative of y(:,j_y) 
  inline double scaleFactorDerY(const MtxInt& der, int j_y=-99999) {
    if(j_y==-99999) j_y=jout;
#ifdef __SURFDATA_ERR_CHECK__
    assert((der.getNRows()==1)&&(der.getNCols()==nvarsr)&&
	   (0<=j_y)&&(j_y<nout));
#endif 
    double ysf=1.0;
    for(int j_xr=0; j_xr<nvarsr; ++j_xr)
      ysf*=pow(unscalexr(0,j_xr),der(0,j_xr));
    return fabs(ysf/unscaley(0,j_y));    
  };

  ///for arbitrary order mixed partial derivative of y(:,j_y) 
  inline double unScaleFactorDerY(const MtxInt& der, int j_y=-99999) {
    if(j_y==-99999) j_y=jout;
#ifdef __SURFDATA_ERR_CHECK__
    assert((der.getNRows()==1)&&(der.getNCols()==nvarsr)&&
	   (0<=j_y)&&(j_y<nout));
#endif 
    double yusf=unscaley(0,j_y);
    for(int j_xr=0; j_xr<nvarsr; ++j_xr)
      yusf/=pow(unscalexr(0,j_xr),der(0,j_xr));
    return fabs(yusf);    
  };

  inline double scaleFactorVarY(int j_y=-99999) {
    if(j_y==-99999) j_y=jout;
#ifdef __SURFDATA_ERR_CHECK__
    assert((0<=j_y)&&(j_y<nout));
#endif 
    return (1.0/(unscaley(0,j_y)*unscaley(0,j_y)));
  };

  inline double unScaleFactorVarY(int j_y=-99999) {
    if(j_y==-99999) j_y=jout;
#ifdef __SURFDATA_ERR_CHECK__
    assert((0<=j_y)&&(j_y<nout));
#endif 
    return (unscaley(0,j_y)*unscaley(0,j_y));
  };

  ///this function should only be called by putPoints() it decides whether to recommend that this SurfData object be rescaled after newpoints2 are added to the current point and returns 0 or 1 to indicate the recommendation; newpoints2 must be non-empty and already have the same scale as this SurfData object.  * 0 means we do not recommend rescaling: 0 is returned if this SurfData object is unscaled, or if the all new points fall within this SurfData object's range of real inputs and outputs. * 1 means that we recommend rescaling: we recommend rescaling if this SurfData object is scaled and at least one of the points in newpoints2 is outside this SurfData object's range of real inputs and outputs.
  int ifRecommendRescale(SurfData& newpoints2);


  //friend void SurfData_scaling_unit_test();
  friend class SurfDataScaler;

public:
  ///real input variables, by default this will be scaled
  MtxDbl xr;

  ///integer input variables, this is never scaled
  MtxInt xi;

  ///output variables, by default this will be scaled
  MtxDbl y;

  //SurfDataScaler *scaler;

  /* //design decision: do we keep labels locked inside SurfData?
  inline const string& getXrLabel(int j) { 
    if((nvarsr<=j) ||  (xrLabel.size()<=j)) 
      assert(0); //need to replace this with a throw

    return xrLabel[j];
  };

  inline const string& getXiLabel(int j) { 
    if((nvarsi<=j) ||  (xiLabel.size()<=j)) 
      assert(0); //need to replace this with a throw

    return xiLabel[j];
  };

  inline const string& getYLabel(int j) { 
    if((nout<=j) ||  (yLabel.size()<=j)) 
      assert(0); //need to replace this with a throw

    return yLabel[j];
  };
  */


  ///tell me how many points there are
  inline int getNPts() const {return npts;};

  ///tell me how many REAL input variables there are
  inline int getNVarsr() const {return nvarsr;};

  ///tell me how many INTEGER input variables there are 
  inline int getNVarsi() const {return nvarsi;};

  ///tell me how many output variables there are
  inline int getNOut() const {return nout;};

  ///tell me the index of the column of y that you want to build the emulator for
  inline int getJOut() const {return jout;};

  ///set the index of the column of y that you want to build the emulator for
  void setJOut(int jout_new){
    jout=jout_new;
    assert((0<=jout)&&(jout<nout));
    return;};
  
  ///retrieve a (as currently scaled) copy of an output dimension, by default this is the output dimension we want to build the emulator for
  inline MtxDbl& getY(MtxDbl& Y, int jout_want=-99999) const {
    if(jout_want==-99999) jout_want=jout;
    if(!((0<=jout_want)&&(jout_want<nout))){
      printf("jout_want=%d nout=%d\n",jout_want,nout);
      assert((0<=jout_want)&&(jout_want<nout));}
    return (y.getCols(Y,jout_want));};

  inline void setUnscaledDomainSize(const MtxDbl& min_max_xr) {
    minMaxXr.copy(min_max_xr);
    ifHaveMinMaxXr=true;
    return;
  };

  inline void unSetUnscaledDomainSize(const MtxDbl& min_max_xr) {
    ifHaveMinMaxXr=false;
    return;
  };

  //returns true if the unscaled domain size has been specified and
  //false if it has not been specified
  inline bool getUnscaledDomainSize(MtxDbl& min_max_xr) const {
    min_max_xr.copy(minMaxXr);
    return (ifHaveMinMaxXr);
  };

  //specify groups of real inputs dimensions within which we want to lock 
  //the aspect ratios of dimensions when using scaleToDefault() (scaling to
  //a hyper rectangle with volume 1)
  inline void setDimGroups(const MtxInt& dim_groups) {
#ifdef __SURFDATA_ERR_CHECK__
    assert((dim_groups.getNRows()==1)&&(dim_groups.getNCols()==nvarsr));
#endif
    lockxr.newSize(2,nvarsr);
    lockxr.putRows(dim_groups,0);
    for(int ivarsr=0; ivarsr<nvarsr; ++ivarsr)
      lockxr(1,ivarsr)=ivarsr; //need to retain original order of dimensions, 
      //because the columns of lockxr are going to be sorted by elements in 
      //row 0 of lockxr
    lockxr.sortCols();
    return;
  };

  //this lets you change your mind if you've previously said you want to
  //lock the aspect ratio of groups of dimensions
  inline void unSetDimGroups(){
    lockxr.clear();
    return;
  };
  
  ~SurfData() {
    clear();
  };


  ///the default constructor, it produces an empty SurfData object
  SurfData(): npts(0), nvarsr(0), nvarsi(0), nout(0), jout(0) {};

  ///a constructor for when there are real input variables that we don't want to group scale and there are no integer input variables, models that use the default scaling will scale the real input variables and output variable(s) to a hypercube of volume 1
  SurfData(const MtxDbl& XR, const MtxDbl& Y, int jout_set=0);

  ///a constructor for when there are real input variables that we want to group scale and there are no integer input variables, models that use the default scaling will scale the real input variables to a hyper-rectangle of volume 1 and the output variable(s) to a hypercube of volume 1
  SurfData(const MtxInt& LOCKXR, const MtxDbl& XR, const MtxDbl& Y, int jout_set=0);

  ///a constructor for when there is real input variables (that we don't want to group scale) and integer input variables, models that use the default scaling will scale the real input variables and output variable(s) to a hypercube of volume 1, it doesn't scale the integer input variables 
  SurfData(const MtxDbl& XR, const MtxInt& XI, const MtxDbl& Y, int jout_set=0);

  ///a constructor for when there is real input variables (that we DO want the model to group scale if it is appropriate for the model) and integer input variables. models that use the default scaling will scale the real input variables to a hyper-rectangle of volume 1 and the output variable(s) to a hypercube of volume 1, it doesn't scale the integer input variables 
  SurfData(const MtxInt& LOCKXR, const MtxDbl& XR, const MtxInt& XI, const MtxDbl& Y, int jout_set=0);

  ///a constructor that reads data from a file for when you don't want the model to group scale the real input variables, models that use the default scaling will scale the real input variables and output variable(s) to a hypercube of volume 1, it doesn't scale the integer input variables  
  SurfData(const string& filename, int nvarsr_in, int nvarsi_in, int nout_in, int jout_in, int skip_columns);

  ///a constructor that reads data from a file for when you DO want the model to group scale the real input variables if it is appropriate for the model to do so. models that use the default scaling will scale the real input variables to a hyper-rectangle of volume 1 and the output variable(s) to a hypercube of volume 1, it doesn't scale the integer input variables
  SurfData(const string& filename, int nvarsr_in, int nvarsi_in, int nout_in, int jout_in, int skip_columns, const MtxInt& LOCKXR);

  ///copy constructor, performs a deep copy
  SurfData(const SurfData& other);

  ///deep copy constructor that keeps only one column of output
  SurfData(const SurfData& other, int jout_keep);

  ///deep copy constructor that keeps only the multiple columns specified in jout_keep
  SurfData(const SurfData& other, const MtxInt& jout_keep);

  //empty the SurfData object
  void clear();

  //perform a deep copy
  SurfData& copy(const SurfData& other);

  ///assignment operator, performs a deep copy
  inline SurfData& operator=(const SurfData& other){return (copy(other));};

  ///int putPoints(SurfData& newpoints, int ipt) puts the single point in newpoints into this SurfData objects ipt-th point, ipt must be greater than or equal to zero, if ipt is less than the current number of points the previous ipt-th point is replaced with the new one, if ipt is greater than the current number of points this function will enlarge the matrices to hold ipt+1 points before inserting this point, if ipt is unspecified it appends the newpoint to the end of this SurfData object's current list of points (after enlarging the matrices), if newpoints contains multiple points and ipt is unspecified it appends all the newpoints to the end of this SurfData obeject's current list of points (if ipt is specified newpoints must contain exactly one point, otherwise you should call int putPoints(SurfData& newpoints, MtxInt& ipts) instead).  If this SurfData object contained any points before this function was called then the newpoints are automatically rescaled to match the current points before they are inserted.  putPoints returns an integer to recommend whether or not the user/programmer should scale this SurfData object after calling this function: * 0 means we do not recommend rescaling: 0 is returned if newpoints was empty, this SurfData object was unscaled before this function was called, or if the all new points fell within the range of real inputs and outputs that this SurfData object had before this function was called. * -1 means that this SurfData object was empty before this function was called so we kept the scaling in newpoints, * 1 means that we recommend rescaling: we recommend rescaling if this SurfData object was scaled before this function was called and at least one of the newpoints real inputs or outputs was outside the range of range of the old  points.  This function does not actually do the rescaling of this SurfData object itself in case the user/programmer specified had scaled to a factor or scaled to a domain, or if for the sake of consistency/comparison-to-the-values-before-the-new-datapoints-were-added he/she doesn't want to rescale.
  int putPoints(SurfData& newpoints, int ipt=-99999);

  ///int putPoints(SurfData& newpoints, MtxInt& ipts) puts the points in newpoints into the points in this SurfData object listed in ipts, all values contained in ipts must be greater than or equal to zero, if any point listed in ipts is less than the current number of points the previous point is replaced with the new one, if any point listed in ipts is greater than the current number of points this function will enlarge the matrices to hold ipts.max()+1 points before inserting these points, if ipts is unspecified then int putPoints(SurfData& newpoints, int ipt) is called instead (with ipt defaulting to -99999 which is the flag to append the new points to the current list of points).  If this SurfData object contained any points before this function was called then the newpoints are automatically rescaled to match the current points before they are inserted. putPoints returns an integer to recommend whether or not the user/programmer should scale this SurfData object after calling this function: * 0 means we do not recommend rescaling: 0 is returned if newpoints was empty, this SurfData object was unscaled before this function was called, or if the all new points fell within the range of real inputs and outputs that this SurfData object had before this function was called. * -1 means that this SurfData object was empty before this function was called so we kept the scaling in newpoints, * 1 means that we recommend rescaling: we recommend rescaling if this SurfData object was scaled before this function was called and at least one of the newpoints real inputs or outputs was outside the range of range of the old  points.  This function does not actually do the rescaling of this SurfData object itself in case the user/programmer specified had scaled to a factor or scaled to a domain, or if for the sake of consistency/comparison-to-the-values-before-the-new-datapoints-were-added he/she doesn't want to rescale.
  int putPoints(SurfData& newpoints, MtxInt& ipts);

  ///retrieve one point (the one with index ipt), pass it back as SurfData rather than a reference to SurfData; this wraps the next function
  inline SurfData getPoints(int ipt){
    SurfData result;
    return getPoints(result,ipt);
  };

  ///retrieve one point (the one with index ipt), return it as a reference to SurfData
  SurfData& getPoints(SurfData& result, int ipt);

  ///retrieve multiple points (the ones whose indices are stored in ipts), pass them back as a SurfData object, this function wraps the next one, but requires an extra copy
  inline SurfData getPoints(MtxInt& ipts){
    SurfData result;
    return getPoints(result,ipts);
  };

  ///retrieve multiple points (the ones whose indices are stored in ipts), pass them back as a reference to a SurfData object
  SurfData& getPoints(SurfData& result, MtxInt& ipts);

  ///retrieve all points except one (the one with index ipt), return them as a SurfData object, this wraps the next function, but requires an extra copy
  inline SurfData excludePoints(int ipt){
    SurfData result;
    return excludePoints(result, ipt);
  };

  ///retrieve all points except one (the one with index ipt), return them as a reference to a SurfData object
  SurfData& excludePoints(SurfData& result, int ipt);

  ///retrieve all points except the ones whose indices are listed in ipts, return them as a SurfData object, this function wraps the next one but requires an extra copy
  inline SurfData excludePoints(MtxInt& ipts){
    SurfData result;
    return excludePoints(result, ipts);
  };

  ///retrieve all points except the ones whose indices are listed in ipts, return them as a reference to a SurfData object
  SurfData& excludePoints(SurfData& result, MtxInt& ipts);

  ///create a split copy of the points stored in this surfdata, one point, the one whose index is listed in ipt, will be placed in "extracted", the rest will be placed in "rest", this is very useful for the PRESS metric
  void extractPoints(SurfData& rest, SurfData& extracted, int ipt);
    
  ///create a split copy of the points stored in this surfdata, the ones whose indices are listed in ipts will be placed in "extracted", the rest will be placed in "rest", this is very useful for the Leave N Out cross Validation metric
  void extractPoints(SurfData& rest, SurfData& extracted, MtxInt& ipts);

  /// Returns true if file has .bspd extension, false if it has .spd extension. Otherwise, an exception is thrown.
  bool hasBinaryFileExtension(const std::string& filename) const;

  /// Read a set of points into SurfData from either a text or binary file.  Opens file and calls either readText() or readBinary(), currently it only reads text
  void read(const string& filename);

  /// Write a set of points in SurfData to either a text or binary file.  Opens file and calls either writeText() or writeBinary(), currently it only writes text
  void write(const string& filename) const;

  ///Set real input variable, xr, labels to xr0 xr1, etc.; integer input variable, xi, labels to xi0, xi1, etc.; output variable, y, labels to y0 y1, etc.
  void defaultLabels();
  
  bool readLabelsIfPresent(string single_line, int skip_columns=0);

  ///read one point into SurfData from a single_line taken from the text .spd file, skip_columns defaults to zero
  void readPointText(int ipt, const string& single_line, int skip_columns=0);
 
  ///Read a set of points into SurfData from an input stream of a text file "*.spd". At this point nvarsr, nvarsi, nout, and optionally lockxr and skip_columns must have already been set by the user. skip_columns defaults to zero
  void readText(istream& is, int skip_columns=0);

  /// Write a set of SurfPoints to an output stream
  void writeText(ostream& os, bool write_labels) const;

  ///Read a set of points into SurfData from an input stream of a binary file "*.bspd". At this point nvarsr, nvarsi, nout, and optionally lockxr and skip_columns must have already been set by the user. skip_columns defaults to zero (need to do something better for variable labels (there should be a way to read them from the binary file)
  void readBinary(istream& is, int skip_columns=0);

  ///read one point into SurfData from a binary .bspd file, skip_columns defaults to zero
  void readPointBinary(int ipt, istream& is, int skip_columns=0);

};

  ///The purpose of the SurfDataScaler class is to provide SurfPackModels and ONLY SurfPackModels with access to some of SurfData's scaling functions, because we don't want anything but a SurfPackModel scaling the SurfData.
class SurfDataScaler{
public:
  SurfDataScaler(SurfData& sd): mySd(sd) {};
private:
  SurfData& mySd;
  inline int isYSingular(int j, double& singular_y) const {return (mySd.isYSingular(j,singular_y));};
  inline int isUnScaled() const {return (mySd.isUnScaled());};
  inline void scaleToDefault(){mySd.scaleToDefault();};
  inline void scaleXrToDomain(MtxDbl& domain_new){mySd.scaleXrToDomain(domain_new);};
  inline void scaleXrToFactor(MtxDbl& unscale_xr){mySd.scaleXrToFactor(unscale_xr);};
  inline void scaleYToFactor(MtxDbl& unscale_y){mySd.scaleYToFactor(unscale_y);};
  inline void scaleToFactors(MtxDbl& unscale_xr, MtxDbl& unscale_y){mySd.scaleToFactors(unscale_xr,unscale_y);};
  inline SurfData& unScaleCopy(SurfData& result){return (mySd.unScaleCopy(result));};
  inline SurfData& unScale(){return (mySd.unScale());};
  inline MtxDbl& scaleXrOther(MtxDbl& xr_other) const {return (mySd.scaleXrOther(xr_other));};
  inline MtxDbl& unScaleXrOther(MtxDbl& xr_other) const {return (mySd.unScaleXrOther(xr_other));};
  inline double scaleYOther(double y, int j=-99999) const {return (mySd.scaleYOther(y,j));};
  inline double unScaleYOther(double y, int j=-99999) const {return (mySd.unScaleYOther(y));};
  inline MtxDbl& scaleYOther(MtxDbl& y, int j=-99999) const {return (mySd.scaleYOther(y,j));};
  inline MtxDbl& unScaleYOther(MtxDbl& y, int j=-99999) const {return (mySd.unScaleYOther(y,j));};
  inline double scaleFactorDerY(int j_xr, int j_y=-99999) const {return (mySd.scaleFactorDerY(j_xr, j_y));};
  inline double unScaleFactorDerY(int j_xr, int j_y=-99999) const {return (mySd.unScaleFactorDerY(j_xr, j_y));};
  inline double scaleFactorDerY(const MtxInt& der, int j_y=-99999) const {return (mySd.scaleFactorDerY(der, j_y));};
  inline double unScaleFactorDerY(const MtxInt& der, int j_y=-99999)const {return (mySd.unScaleFactorDerY(der, j_y));};
  inline double scaleFactorVarY(int j_y=-99999) const {return (mySd.scaleFactorVarY(j_y));};
  inline double unScaleFactorVarY(int j_y=-99999) const {return (mySd.unScaleFactorVarY(j_y));};

  friend class SurfPackModel;
  friend class KrigingModel;
};

} // end namespace nkm

#endif
