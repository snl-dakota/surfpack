#include "NKM_SurfData.hpp"
#include "NKM_SurfPack.hpp"
#include <iostream>
#include <iomanip>

namespace nkm {

using namespace std;
using std::ios;
using std::istringstream;
using std::istream;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;

/***********************************************************************/
/***********************************************************************/
/**** Unit Test functions for SurfData start here                   ****/
/***********************************************************************/
/***********************************************************************/

//#define __SURFDATA_SCALING_UNIT_TEST__

#ifdef __SURFDATA_SCALING_UNIT_TEST__

void SurfData_scaling_unit_test(){
  int npts=100;
  int nvarsr=7; //this is hard coded below
  int nout=3; //this is hard coded below
  MtxInt lockxr(1,nvarsr);
  MtxInt sortedlockxr(2,nvarsr); //what it should be after being sorted
  //original    sorted into order   original position
  lockxr(0,0)=1; sortedlockxr(0,0)=0; sortedlockxr(1,0)=1; 
  lockxr(0,1)=0; sortedlockxr(0,1)=1; sortedlockxr(1,1)=0;
  lockxr(0,2)=2; sortedlockxr(0,2)=1; sortedlockxr(1,2)=3;
  lockxr(0,3)=1; sortedlockxr(0,3)=1; sortedlockxr(1,3)=5; 
  lockxr(0,4)=2; sortedlockxr(0,4)=2; sortedlockxr(1,4)=2;
  lockxr(0,5)=1; sortedlockxr(0,5)=2; sortedlockxr(1,5)=4;
  lockxr(0,6)=3; sortedlockxr(0,6)=3; sortedlockxr(1,6)=6;
  double SingDimValue=-5.3; //real input dimension 3 is singular
  MtxDbl minmaxd(2,nvarsr);
  MtxDbl minmaxD(2,nvarsr);
  MtxDbl XR(npts,nvarsr);  
  MtxDbl Y(npts,nout);
  MtxDbl domain_new(2,nvarsr);
  domain_new(0,0) = 0.3; domain_new(1,0)= 0.7;
  domain_new(0,1) =-0.3; domain_new(1,1)= 5.0;
  domain_new(0,2) = 5.0; domain_new(1,2)= 8.3;
  domain_new(0,3) =-7.0; domain_new(1,3)=-4.9;
  domain_new(0,4) = 7.0; domain_new(1,4)= 9.4;
  domain_new(0,5) = 4.9; domain_new(1,5)= 7.0;
  domain_new(0,6) =-2.0; domain_new(1,6)= 2.0;

  int maxi, mini, modi;
  double tempd, mind, maxd;
  int i,j;
  for(j=0; j<nvarsr; j++) {
    mini=rand()%10+5;
    maxi=rand()%20+20;
    modi=5*(maxi-mini);
    for(i=0;i<npts;i++)
      XR(i,j)=(rand()%modi)*0.2+mini;
  }
  //make one of the real input dimensions (dimension 3) singular
  for(i=0; i<npts; i++)
    XR(i,3)=SingDimValue;

  for(j=0; j<nout; j++) {
    mini=rand()%10+5;
    maxi=rand()%20+20;
    modi=5*(maxi-mini);
    for(i=0;i<npts;i++)
      Y(i,j)=(rand()%modi)*0.2+mini;
  }


  //test1 is for default scaling WITHOUT locking aspect ratios of dimensions 
  SurfData test1a(XR,Y);
  SurfData test1b=test1a; test1b.unScale();
  SurfData test1c=test1b; test1c.scaleToDefault();
  SurfData test1d(XR,Y,0,'U'); //leave it unscaled
  SurfData test1e=test1a; test1e.scaleXrToDomain(domain_new);
  SurfData test1f=test1d; test1f.scaleXrToDomain(domain_new);
  SurfData test1g=test1d; test1g.scaleToFactors(test1e.unscalexr,test1e.unscaley);
  SurfData test1h=test1a; test1h.scaleXrToFactor(test1e.unscalexr);

  //test the xr scaling
  for(j=0;j<nvarsr;j++) {
    if(j==3){
      //check to make sure that dimension 3 is flagged as being singular
      assert((test1a.unscalexr(0,3)==-1.0)&&
	     (test1a.unscalexr(1,3)==SingDimValue));
      for(i=0;i<npts;i++)
	assert(test1a.xr(i,3)==0.0);
    }
    else{
      assert(test1a.unscalexr(0,j)>0.0);
      mind=maxd=test1a.xr(0,j);
      for(i=1;i<npts;i++){
	if(test1a.xr(i,j)<mind) mind=test1a.xr(i,j);
	if(test1a.xr(i,j)>maxd) maxd=test1a.xr(i,j);
      }
      assert(if_close_enough(mind,-0.5)&&if_close_enough(maxd,0.5));
      mind=maxd=test1e.xr(0,j);
      for(i=1;i<npts;i++){
	if(test1e.xr(i,j)<mind) mind=test1e.xr(i,j);
	if(test1e.xr(i,j)>maxd) maxd=test1e.xr(i,j);
      }
      assert(if_close_enough(mind,domain_new(0,j))&&
	     if_close_enough(maxd,domain_new(1,j)));
    }
    assert(if_close_enough(test1a.unscalexr(0,j),test1c.unscalexr(0,j))&&
	   if_close_enough(test1a.unscalexr(1,j),test1c.unscalexr(1,j))&&
	   if_close_enough(test1b.unscalexr(0,j),test1d.unscalexr(0,j))&&
	   if_close_enough(test1b.unscalexr(1,j),test1d.unscalexr(1,j))&&
	   (test1d.unscalexr(0,j)==1.0)&&(test1d.unscalexr(1,j)==0.0)&&
	   if_close_enough(test1e.unscalexr(0,j),test1f.unscalexr(0,j))&&
	   if_close_enough(test1e.unscalexr(1,j),test1f.unscalexr(1,j))&&
	   if_close_enough(test1e.unscalexr(0,j),test1g.unscalexr(0,j))&&
	   if_close_enough(test1e.unscalexr(1,j),test1g.unscalexr(1,j)));
	   
    for(i=0;i<npts;i++)
      assert(if_close_enough(test1a.xr(i,j),test1c.xr(i,j))&&
	     if_close_enough(XR(i,j),test1b.xr(i,j))&&
	     (XR(i,j)==test1d.xr(i,j))&&
	     if_close_enough(test1e.xr(i,j),test1f.xr(i,j))&&
	     if_close_enough(test1e.xr(i,j),test1g.xr(i,j))&&
	     if_close_enough(test1e.xr(i,j),test1h.xr(i,j)));
  }
  
  //test the y scaling
  for(j=0;j<nout;j++) {
    mind=maxd=test1a.y(0,j);
    for(i=1;i<npts;i++){
      if(test1a.y(i,j)<mind) mind=test1a.y(i,j);
      if(test1a.y(i,j)>maxd) maxd=test1a.y(i,j);
    }
    assert(if_close_enough(mind,-0.5)&&
	   if_close_enough(maxd,0.5)&&
	   if_close_enough(test1a.unscaley(0,j),test1c.unscaley(0,j))&&
	   if_close_enough(test1a.unscaley(1,j),test1c.unscaley(1,j))&&
	   if_close_enough(test1b.unscaley(0,j),test1d.unscaley(0,j))&&
	   if_close_enough(test1b.unscaley(1,j),test1d.unscaley(1,j))&&
	   (test1d.unscaley(0,j)==1.0)&&(test1d.unscaley(1,j)==0.0));

    for(i=0;i<npts;i++){
      assert(if_close_enough(test1a.y(i,j),test1c.y(i,j))&&
	     if_close_enough(Y(i,j),test1b.y(i,j))&&
	     (Y(i,j)==test1d.y(i,j))&&
      	     if_close_enough(test1a.y(i,j),test1g.y(i,j)));
      //      assert(if_close_enough(test1a.y(i,j),test1c.y(i,j)));
      //      assert(if_close_enough(Y(i,j),test1b.y(i,j)));
      //      assert((Y(i,j)==test1d.y(i,j)));
      //      assert(if_close_enough(test1a.y(i,j),test1g.y(i,j)));
    }

  }

  //test2 is for default scaling WITH locking of aspect ratios of groups of dimensions
  SurfData test2a(lockxr,XR,Y);
  SurfData test2b=test2a; test2b.unScale();
  SurfData test2c=test2b; test2c.scaleToDefault();
  SurfData test2d(lockxr,XR,Y,0,'U'); //leave it unscaled

  for(j=0;j<nvarsr;j++)
    for(i=0; i<2; i++)
      assert((sortedlockxr(i,j)==test2a.lockxr(i,j))&&
	     (sortedlockxr(i,j)==test2b.lockxr(i,j))&&
	     (sortedlockxr(i,j)==test2c.lockxr(i,j))&&
	     (sortedlockxr(i,j)==test2d.lockxr(i,j)));

  //test the xr scaling
  for(j=0;j<nvarsr;j++) {
    if(j==3){
      //check to make sure that dimension 3 is flagged as being singular
      assert((test2a.unscalexr(0,3)==-1.0)&&
	     (test2a.unscalexr(1,3)==SingDimValue));
      for(i=0;i<npts;i++)
	assert(test2a.xr(i,3)==0.0);
    }
    else assert(test2a.unscalexr(0,j)>0.0);

    minmaxd(0,j)=minmaxd(1,j)=test2a.xr(0,j);
    minmaxD(0,j)=minmaxD(1,j)=XR(0,j);
    for(i=1;i<npts;i++){
      if(test2a.xr(i,j)<minmaxd(0,j)) minmaxd(0,j)=test2a.xr(i,j);
      if(test2a.xr(i,j)>minmaxd(1,j)) minmaxd(1,j)=test2a.xr(i,j);
      if(XR(i,j)<minmaxD(0,j)) minmaxD(0,j)=XR(i,j);
      if(XR(i,j)>minmaxD(1,j)) minmaxD(1,j)=XR(i,j);
    }
  
    assert(if_close_enough(test2a.unscalexr(0,j),test2c.unscalexr(0,j))&&
	   if_close_enough(test2a.unscalexr(1,j),test2c.unscalexr(1,j)));
    for(i=0;i<npts;i++)
      assert(if_close_enough(test2a.xr(i,j),test2c.xr(i,j))&&
	     if_close_enough(XR(i,j),test2b.xr(i,j))&&
	     (XR(i,j)==test2d.xr(i,j)));

  }
  assert(if_close_enough(1.0,(minmaxd(1,0)-minmaxd(0,0))*
			 (minmaxd(1,5)-minmaxd(0,5)))&&
	 if_close_enough((minmaxD(1,0)-minmaxD(0,0))/(minmaxd(1,0)-minmaxd(0,0)),
			 (minmaxD(1,5)-minmaxD(0,5))/(minmaxd(1,5)-minmaxd(0,5)))&&
	 if_close_enough(1.0,(minmaxd(1,2)-minmaxd(0,2))*
			 (minmaxd(1,4)-minmaxd(0,4)))&&
	 if_close_enough((minmaxD(1,2)-minmaxD(0,2))/(minmaxd(1,2)-minmaxd(0,2)),
			 (minmaxD(1,4)-minmaxD(0,4))/(minmaxd(1,4)-minmaxd(0,4)))&&
	 if_close_enough(1.0,(minmaxd(1,1)-minmaxd(0,1)))&&
	 if_close_enough(1.0,(minmaxd(1,6)-minmaxd(0,6)))&&
	 if_close_enough(0.0,(minmaxd(1,3)-minmaxd(0,3))));
	 

  //test the y scaling
  for(j=0;j<nout;j++) {
    mind=maxd=test2a.y(0,j);
    for(i=1;i<npts;i++){
      if(test2a.y(i,j)<mind) mind=test2a.y(i,j);
      if(test2a.y(i,j)>maxd) maxd=test2a.y(i,j);
    }
    assert(if_close_enough(mind,-0.5)&&
	   if_close_enough(maxd,0.5)&&
	   if_close_enough(test2a.unscaley(0,j),test2c.unscaley(0,j))&&
	   if_close_enough(test2a.unscaley(1,j),test2c.unscaley(1,j)));
    for(i=0;i<npts;i++)
      assert(if_close_enough(test2a.y(i,j),test2c.y(i,j))&&
	     if_close_enough(Y(i,j),test2b.y(i,j)));
  }

  printf("class SurfData; passed the scaling unit tests\n");
  return;
}

int main(){
  SurfData_scaling_unit_test();
  return 0;
};



#endif

/***********************************************************************/
/***********************************************************************/
/**** Unit Test functions for SurfData end here                     ****/
/***********************************************************************/
/***********************************************************************/



/***********************************************************************/
/***********************************************************************/
/**** SurfData member functions start here                          ****/
/***********************************************************************/
/***********************************************************************/
  

///a constructor for when there are real input variables that we don't want the model to group scale and there are no integer input variables, this will result in the model scaling its real input variables and output variable(s) to a hypercube of volume 1
  SurfData::SurfData(const MtxDbl& XR, const MtxDbl& Y, int jout_set) : npts(XR.getNRows()), nvarsr(XR.getNCols()), nvarsi(0), nout(Y.getNCols()), jout(jout_set), ifHaveMinMaxXr(false)
{
  npts  =XR.getNRows();
  nvarsr=XR.getNCols();
  assert(Y.getNRows()==npts);
  nvarsi=0;
  nout  =Y.getNCols();
  jout  =jout_set;
  
  if(0<npts) {
    assert((0<=jout)&&(jout<nout));
    xr=XR;
    y=Y;
    dontScale();
  }
  else{
    jout=0;
    cerr << "Warning: SurfData() constructor was passed empty data matrices!!!" << endl;
  }
  defaultLabels();
  return;
}

  

///a constructor for when there are real input variables that we want the model to group scale (if it is appropriate to the model) and there are no integer input variables.  If it is appropriate to the model to do so, the model will automatically scale the real input variables to a hyper-rectangle of volume 1 and the output variable(s) to a hypercube of volume 1
SurfData::SurfData(const MtxInt& LOCKXR, const MtxDbl& XR, const MtxDbl& Y, int jout_set) : npts(XR.getNRows()), nvarsr(XR.getNCols()), nvarsi(0), nout(Y.getNCols()), jout(jout_set), ifHaveMinMaxXr(false)
{
  //npts  =XR.getNRows();
  //nvarsr=XR.getNCols();
  assert((LOCKXR.getNElems()==nvarsr)&&
	 (Y.getNRows()==npts));
  //nvarsi=0;
  //nout  =Y.getNCols();
  //jout  =jout_set;
  
  if(0<npts) {
    assert((0<=jout)&&(jout<nout));
    xr=XR;
    y=Y;
    lockxr.newSize(2,nvarsr);
    int j;
    for(j=0; j<nvarsr; j++) {
      lockxr(0,j)=LOCKXR(j);
      lockxr(1,j)=j; //need to retain original order of dimensions, because the columns of lockxr are going to be sorted by elements in row 0 of lockxr
    }
    //you want us to scale by groups, and although there are group flags for each dimension (lockxr), those dimensions aren't arranged so that member dimensions of a group are listed sequentionally so we need to sort lockxr by groups.
    lockxr.sortCols();
    dontScale();
  }
  else{
    jout=0;
    cerr << "Warning: SurfData() constructor was passed empty data matrices!!!" <<endl;
  }
  defaultLabels();
  return;
}


///a constructor for when there is real input variables (that we don't want to group scale) and integer input variables.  Models that use the default scaling will automatically scale the real input variables and output variable(s) to a hypercube of volume 1, they doesn't scale the integer input variables 
SurfData::SurfData(const MtxDbl& XR, const MtxInt& XI, const MtxDbl& Y, int jout_set) : npts(XR.getNRows()), nvarsr(XR.getNCols()), nvarsi(XI.getNCols()), nout(Y.getNCols()), jout(jout_set), ifHaveMinMaxXr(false)
{
  //npts  =XR.getNRows();
  //nvarsr=XR.getNCols();
  assert((XI.getNRows()==npts)&&(Y.getNRows()==npts));
  //nvarsi=XI.getNCols();
  //nout  =Y.getNCols();
  //jout  =jout_set;
  
  if(0<npts) {
    assert((0<=jout)&&(jout<nout));
    xr=XR;
    y=Y;
    dontScale();
    xi=XI; //note that integer input variables will never be scaled
  }
  else{
    jout=0;
    cerr << "Warning: SurfData() constructor was passed empty data matrices!!!" <<endl;
  }
  defaultLabels();
  return;
}


///a constructor for when there are real input variables (that we DO want the model to group scale if it is appropriate for the model) and integer input variables.  Models that use the default scaling will automatically scale the real input variables to a hyper-rectangle of volume 1 and the output variable(s) to a hypercube of volume 1, the integer input variables won't be scaled
  SurfData::SurfData(const MtxInt& LOCKXR, const MtxDbl& XR, const MtxInt& XI, const MtxDbl& Y, int jout_set) : npts(XR.getNRows()), nvarsr(XR.getNCols()), nvarsi(XI.getNCols()), nout(Y.getNCols()), jout(jout_set), ifHaveMinMaxXr(false)
{
  //npts  =XR.getNRows();
  //nvarsr=XR.getNCols();
  assert((LOCKXR.getNElems()==nvarsr)&&
	 (XI.getNRows()==npts)&&(Y.getNRows()==npts));
  //nvarsi=XI.getNCols();
  //nout  =Y.getNCols();
  //jout  =jout_set;
  
  if(0<npts) {
    assert((0<=jout)&&(jout<nout));
    xr=XR;
    y=Y;
    lockxr.newSize(2,nvarsr);
    int j;
    for(j=0; j<nvarsr; j++) {
      lockxr(0,j)=LOCKXR(j);
      lockxr(1,j)=j; //need to retain original order of dimensions, because the columns of lockxr are going to be sorted by elements in row 0 of lockxr
    }
    //you want us to scale by groups, and although there are group flags for each dimension (lockxr), those dimensions aren't arranged so that member dimensions of a group are listed sequentionally so we need to sort lockxr by groups.
    lockxr.sortCols();
    dontScale();
    xi=XI; //note that integer input variables will never be scaled
  }
  else{
    jout=0;
    cerr << "Warning: SurfData() constructor was passed empty data matrices!!!" << endl;
  }
  defaultLabels();
  return;
}

///a constructor that reads data from a file for when you don't want the model to group scale the real input variables
  SurfData::SurfData(const string& filename, int nvarsr_in, int nvarsi_in, int nout_in, int jout_in, int skip_columns) : nvarsr(nvarsr_in), nvarsi(nvarsi_in), nout(nout_in), jout(jout_in), ifHaveMinMaxXr(false)
{
  //cout << "filename=[" << filename <<"]\nnvarsr=" << nvarsr 
  //   <<"\nnvarsi=" << nvarsi <<"\nnout="<< nout << "\njout=" << jout 
  //   <<"\nskip_columns=" << skip_columns << "\nifscale=[" << ifscale <<"]\n";

  assert((nvarsr>0)&&(nvarsi>=0)&&(nout>0)&&(jout>=0));

  //lockxr is an empy matrix
  dontScale(); //set the initial scaling to identity
  read(filename);

  return;
}

///a constructor that reads data from a file for when you DO want the model to group scale the real input variables if it is appropriate for the model
  SurfData::SurfData(const string& filename, int nvarsr_in, int nvarsi_in, int nout_in, int jout_in, int skip_columns, const MtxInt& LOCKXR) : nvarsr(nvarsr_in), nvarsi(nvarsi_in), nout(nout_in), jout(jout_in), ifHaveMinMaxXr(false)
{
  assert((nvarsr>0)&&(nvarsi>=0)&&(nout>0)&&(jout>=0)&&
	 (LOCKXR.getNElems()==nvarsr));
  lockxr.newSize(2,nvarsr);
  int j;
  for(j=0; j<nvarsr; ++j) {
    lockxr(0,j)=LOCKXR(j);
    lockxr(1,j)=j; //need to retain original order of dimensions, because the columns of lockxr are going to be sorted by elements in row 0 of lockxr
  }
  //you want us to scale by groups, and although there are group flags for each dimension (lockxr), those dimensions aren't arranged so that member dimensions of a group are listed sequentionally so we need to sort lockxr by groups.
  lockxr.sortCols();

  dontScale(); //set the initial scaling to identity  
  read(filename);

  return;
}


///copy constructor performs a deep copy
SurfData::SurfData(const SurfData& other) : npts(other.npts), nvarsr(other.nvarsr), nvarsi(other.nvarsi), nout(other.nout), jout(other.jout), xr(other.xr), xi(other.xi), y(other.y), unscalexr(other.unscalexr), unscaley(other.unscaley), lockxr(other.lockxr), ifHaveMinMaxXr(false)
 //effective c++ says to initialize rather than assign 
{
  //I don't know if vectors have copy constructors so use the assignment operator which I know does the right thing.
  xrLabels=other.xrLabels;
  xiLabels=other.xiLabels;
  yLabels =other.yLabels;
  //copy(other);
  return;
}

///deep copy constructor that keeps only one column of output
SurfData::SurfData(const SurfData& other, int jout_keep) : npts(other.npts), nvarsr(other.nvarsr), nvarsi(other.nvarsi), nout(1), jout(0), ifHaveMinMaxXr(false), xr(other.xr), xi(other.xi), unscalexr(other.unscalexr), lockxr(other.lockxr) //effective c++ says to initialize rather than assign 
{
  //printf("inside SurfData::SurfData(const SurfData& other, int jout_keep)\n");  fflush(stdout);

  assert((0<=jout_keep)&&(jout_keep<other.nout)); //we shouldn't "throw" exceptions during construction but there is no prohibition of asserts during construction
  if(jout_keep==-1)
    jout_keep=other.jout;
  y.newSize(nout);
  //printf("other.getNOut()=%d other.unscaley.getNCols()=%d\n", other.getNOut(),other.unscaley.getNCols());  fflush(stdout);
  other.y.getCols(y,jout_keep);
  other.unscaley.getCols(unscaley,jout_keep);
  //printf("about to get labels\n"); fflush(stdout);
  xrLabels=other.xrLabels;
  xiLabels=other.xiLabels;
  //printf("got xrLabels and xiLabels, about to get yLabels\n"); fflush(stdout);
  yLabels.resize(1);
  yLabels[0]=other.yLabels[jout_keep];
  //printf("leaving SurfData::SurfData(const SurfData& other, int jout_keep)\n"); fflush(stdout);
  return;
}
/*
///deep copy constructor that keeps only the multiple columns specified in jout_keep
SurfData::SurfData(const SurfData& other, const MtxInt& jout_keep) : npts(other.npts), nvarsr(other.nvarsr), nvarsi(other.nvarsi), nout(jout_keep.getNElems()), jout(0), ifHaveMinMaxXr(false), xr(other.xr), xi(other.xi), unscalexr(other.unscalexr), lockxr(other.lockxr) //effective c++ says to initialize rather than assign 
{
  assert((0<=jout_keep.minElem())&&(jout_keep.maxElem()<other.nout)); //we shouldn't "throw" exceptions during construction but there is no prohibition of asserts during construction

  y.newSize(npts,nout);
  other.y.getCols(y,jout_keep);
  other.unscaley.getCols(unscaley,jout_keep);
  xrLabels=other.xrLabels;
  xiLabels=other.xiLabels;
  yLabels.resize(nout); //yLabels is a vector of strings not a SurfMat object
  for(int j=0; j<nout; ++j) {
    yLabels[j]=other.yLabels[jout_keep(j)];
    if(other.jout==jout_keep(j))
      jout=j; //keep the other's current output as this copies output
  }
  return;
}
*/

void SurfData::clear(){
  npts=nvarsr=nvarsi=nout=jout=0;
  xr.clear();
  xi.clear();
  y.clear();
  unscalexr.clear();
  unscaley.clear();
  lockxr.clear();
  xrLabels.clear();
  xiLabels.clear();
  yLabels.clear();
}


///make a deep copy
SurfData& SurfData::copy(const SurfData& other) {
  npts     =other.npts;
  nvarsr   =other.nvarsr;
  nvarsi   =other.nvarsi;
  nout     =other.nout;
  jout     =other.jout;
  xr       =other.xr;
  xi       =other.xi;
  y        =other.y;
  unscalexr=other.unscalexr;
  unscaley =other.unscaley;
  lockxr   =other.lockxr;
  xrLabels =other.xrLabels;
  xiLabels =other.xiLabels;
  yLabels  =other.yLabels;
  return *this;
}



/// Returns true if file has .bspd extension, false if it has .spd extension. Otherwise, an exception is thrown.
bool SurfData::hasBinaryFileExtension(const string& filename) const
{
  if (surfpack::hasExtension(filename,".bspd")) {
    return true;
  } else if (surfpack::hasExtension(filename,".spd")) {
    return false;
  } else if (surfpack::hasExtension(filename,".dat")) {
    return false;
  } else {
    throw surfpack::io_exception(
      "Unrecognized filename extension.  Use .bspd or .spd"
    );
  }
}

/// Read a set of points into SurfData from either a text or binary file.  Opens file and calls either readText() or readBinary()
void SurfData::read(const string& filename)
{
  // Open file in binary or text mode based on filename extension (.bspd or .spd)
  bool binary = hasBinaryFileExtension(filename);
  //cout << "binary=" << binary << endl;
  ifstream infile(filename.c_str(), (binary ? ios::in|ios::binary : ios::in));
  if (!infile) {
    //cout << "couldn't open file" << endl;
    throw surfpack::file_open_failure(filename);
  } else if (binary) {
    cout << "attempting to open a binary file" << endl;
    assert(0);
    readBinary(infile);
  } else {
    //cout << "attempting to open a text file" << endl;
    readText(infile);
  }
  //cout << "done reading from file\n";
  // Object may have already been created
  infile.close();
}

void SurfData::write(const string& filename) const
{
  //if (mapping.empty()) {
  //ostringstream errormsg;
  //errormsg << "Cannot write SurfData object to stream."
  //     << "  No active data points." << endl;
  //throw bad_surf_data(errormsg.str());
  //}
  bool binary = hasBinaryFileExtension(filename);
  ofstream outfile(filename.c_str(), 
    (binary ? ios::out|ios::binary : ios::out));
  if (!outfile) {
    throw surfpack::file_open_failure(filename);
  } else if (binary) {
    cout << "attempting to write a binary file" <<endl;
    assert(0);
    //writeBinary(outfile);
  } else {
    // Write the header and label info for .spd, not for .dat
    bool write_labels = surfpack::hasExtension(filename,".spd");
    writeText(outfile, write_labels);
  }
  outfile.close();
}


/// Set real input variable, xr, labels to xr0 xr1, etc.; integer input variable, xi, labels to xi0, xi1, etc.; output variable, y, labels to y0 y1, etc.
void SurfData::defaultLabels()
{
  int j;
  //real input variable labels
  xrLabels.resize(nvarsr);
  for(j=0; j<nvarsr; ++j) {
    ostringstream os;
    os << "xr" << j ;
    xrLabels[j] = os.str();
  }

  //integer input variable labels
  xiLabels.resize(nvarsi);
  for(j=0; j<nvarsi; ++j) {
    ostringstream os;
    os << "xi" << j;
    xiLabels[j] = os.str();
  }

  //output variable labels
  yLabels.resize(nout);
  for(j=0; j<nout; ++j) {
    ostringstream os;
    os << "y" << j ;
    yLabels[j] = os.str();
  }
}

bool SurfData::readLabelsIfPresent(string single_line, int skip_columns)
{
  if(! ((single_line[0] == '%')||(single_line[0] == '#'))) {
    defaultLabels();
    return false;
  } else { // use custom labels

    single_line[0] = ' ';
    string dummy;
    xrLabels.resize(nvarsr); //this is a vector of strings
    xiLabels.resize(nvarsi); //this is a vector of strings
    yLabels.resize(nout); //this is a vector of strings
    istringstream is(single_line);
    int j;

    //labels for columns the user told us to skip
    for(j=0; j<skip_columns; ++j) {
      is >> dummy;
      if(dummy == "") {
	// Not enough heading names in the line of column headings
	// Use the default headings and return
	defaultLabels();
	return false;
      }
    }

    //real input variable labels
    for(j=0; j<nvarsr; ++j) {
      is >> xrLabels[j];
      if(xrLabels[j] == "") { 
	// Not enough heading names in the line of column headings
	// Use the default headings and return
        defaultLabels();
        return false;
      }
    } 

    //integer input variable labels
    for(j=0; j<nvarsi; ++j) {
      is >> xiLabels[j];
      if(xiLabels[j] =="") {
	// Not enough heading names in the line of column headings
	// Use the default headings and return
	defaultLabels();
	return false;
      }
    }

    //output variable labels
    for(j=0; j<nout; ++j) {
      is >> yLabels[j];
      if (yLabels[j] == "") { 
	// Not enough heading names in the line of column headings
	// Use the default headings and return
        defaultLabels();
        return false;
      }
    } 

  } // use custom labels
  return true;
}

///read one point into SurfData from a single_line taken from the text .spd file, skip_columns defaults to zero
void SurfData::readPointText(int ipt, const string& single_line, 
			     int skip_columns)
{
  int nvarsr_read=0, nvarsi_read=0, nout_read=0, nskip_read=0;
  string dummy;
  try {
    // read the point as text
    istringstream streamline(single_line);

    //skip leading columns
    for(nskip_read=0; nskip_read<skip_columns; ++nskip_read) {
      // Throw an exception if there are fewer values on this line that
      // expected.
      surfpack::checkForEOF(streamline);
      streamline >> dummy;
    }

    //read in real input variables "xr"
    for(nvarsr_read=0; nvarsr_read<nvarsr; ++nvarsr_read) {
      // Throw an exception if there are fewer values on this line that
      // expected.
      surfpack::checkForEOF(streamline);
      streamline >> xr(ipt,nvarsr_read);
    }

    //read in integer input variables "xi"
    for(nvarsi_read=0; nvarsi_read<nvarsi; ++nvarsi_read) {
      // Throw an exception if there are fewer values on this line that
      // expected.
      surfpack::checkForEOF(streamline);
      streamline >> xi(ipt,nvarsi_read);
    } 

    //read in output variables "y"
    for(nout_read=0; nout_read<nout; ++nout_read) {
      // Throw an exception if there are fewer values on this line that
      // expected.
      surfpack::checkForEOF(streamline);
      streamline >> y(ipt,nout_read);
    } 
  } catch(surfpack::io_exception&) {
    cerr << "Bad SurfPoint: " << single_line 
	 << "\nExpected on this line: " 
         << "\n  " << skip_columns << " leading columns to skip and"
	 << "\n  " << nvarsr << " real input variable(s) and"
	 << "\n  " << nvarsi << " integer input variable(s) and"
         << "\n  " << nout << " output variables(s) " << endl
         << "Found: " 
	 << "\n  " << nskip_read << " leading columns to skip and"
	 << "\n  " << nvarsr_read << " real input variable(s) and"
	 << "\n  " << nvarsi_read << " integer input variable(s) and"
         << "\n  " << nout_read << " output variables(s) " << endl;
    throw;
  } catch (...) {
    cerr << "Exception caught and rethrown in SurfData::readPointText(...)" << endl;
    throw;
  }
  return;
}
    

///read one point into SurfData from a binary .bspd file, skip_columns defaults to zero
void SurfData::readPointBinary(int ipt, istream& is, int skip_columns)
{
  cout << "SurfData: reading from a binary file has not yet been implemented NEEDS MUCH WORK must deal with cross platform endian-ness variation\n";
  assert(0);

  int nvarsr_read=0, nvarsi_read=0, nout_read=0, nskip_read=0;
  double dummy;

  try {
    // read the point in binary format
    for (nvarsr_read=0; nvarsr_read<nvarsr; ++nvarsr_read) {
       // Throw an exception if there are fewer values on this line that
       // expected.
       surfpack::checkForEOF(is);
       is.read(reinterpret_cast<char*>(xr.ptr(ipt,nvarsr_read)),sizeof(xr(ipt,nvarsr_read)));
    }
    for (nout_read=0; nout_read<nout; ++nout_read) {
       // Throw an exception if there are fewer values on this line that
       // expected.
       surfpack::checkForEOF(is);
       is.read(reinterpret_cast<char*>(y.ptr(ipt,nout_read)),sizeof(y(ipt,nout_read)));
    }
  } catch (surfpack::io_exception&) {
    cerr << "Bad SurfPoint: binary file"  
	 << "\nExpected on this line: " 
         << "\n  " << skip_columns << " leading columns to skip and"
	 << "\n  " << nvarsr << " real input variable(s) and"
	 << "\n  " << nvarsi << " integer input variable(s) and"
         << "\n  " << nout << " output variables(s) " << endl
         << "Found: " 
	 << "\n  " << nskip_read << " leading columns to skip and"
	 << "\n  " << nvarsr_read << " real input variable(s) and"
	 << "\n  " << nvarsi_read << " integer input variable(s) and"
         << "\n  " << nout_read << " output variables(s) " << endl;
    throw;
  } catch (...) {
    cerr << "Exception rethrown in SurfData::readPointBinary(...)" 
         << endl;
    throw;
  }
}







/** Read a set of points into SurfData from an input stream of a text file 
    "*.spd". At this point nvarsr, nvarsi, nout, and optionally lockxr and 
    skip_columns must have already been set by the user. skip_columns 
    defaults to zero
*/
void SurfData::readText(istream& is, int skip_columns) 
{
  string single_line;
  int nlines=0;
  npts=0;
  try {
    //determine how many lines are in the file, so we can allocate matrices if any lines are blank or commented (such as a variable label line) then we will need to resize() matrices (not reshape() not newSize()) at the end
    while(!(is.eof())) {
      nlines++;
      getline(is,single_line);
    }
    //cout << "nlines=" << nlines << endl;
    assert(nlines&&nvarsr&&nout); //replace with a
    //if(nlines==0) throw;
    xr.newSize(nlines,nvarsr);
    xi.newSize(nlines,nvarsi);
    y.newSize(nlines,nout);
    //cout << "just newsized arrays\n";

    cout << "TODO in SurfData.cpp: void SurfData::readText(istream&is, int skip_columns)  need to check for \"failbit\" and \"badbit\" before doing \"is.clear()\"\n";
    is.clear();
    is.seekg(0, ios::beg);
    
    //cout << "I think I just reset file position to beginning of file\n";
    getline(is,single_line);
    //cout << "firstline is=[" << single_line << "]\n";
    istringstream streamline(single_line);
    if (!readLabelsIfPresent(single_line, skip_columns)) {
      if ((single_line != "") && 
	  (single_line != "\n") && 
	  (single_line[0] != '%')&&
	  (single_line[0] != '#')) {
	readPointText(npts, single_line, skip_columns);
        npts = 1;
      }
      //cout << "I didn't find labels; nout=" << nout << endl;
    }

    //we don't know whether npts=0 or 1 will occur at this line during 
    //execution (Could be either depending on input file) so we have to 
    //program generically => separate ilines and npts
    for(int iline=1; iline<nlines; ++iline) {
      //cout << "iline=" << iline << "  npts=" << npts << endl;
      getline(is,single_line);
      if ((single_line != "") && 
	  (single_line != "\n") && 
	  (single_line[0] != '%')&&
	  (single_line[0] != '#')) {
	readPointText(npts, single_line, skip_columns);
	++npts; //only increase npts when we add a new point
      }
    }
  } catch(surfpack::io_exception& exception) {
    //cout << "encountered an exception in readText()\n";
    cerr << exception.what() << endl;
    throw;
  } 
  if(npts < nlines){
    //cout << "downsizing arrays because npts=" <<npts << " < nlines=" << nlines <<endl; 
    //there were one or more blank or commented lines so we will remove the estra lines from the input/output variable matrices with the resize() command
    xr.resize(npts,nvarsr);
    xi.resize(npts,nvarsi);
    y.resize(npts,nout);
  }else if(npts > nlines) {
    assert(0);  //replace with a throw, the only way to get here should be that whatever called this function passed it an istream that was not at the beginning of the .spd file, in which case we've read values into memory locations beyond the bounds of the matrices, i.e. it's a segfault waiting to happen
  }
}

/// Write a set of SurfPoints to an output stream
void SurfData::writeText(ostream& os, bool write_labels) const
{
  assert((nvarsr>=0)&&(nvarsi>=0)&&(nvarsr+nvarsi>0)&&(nout>=1));
  int ivarsr, ivarsi, iout;
  stringstream ss; //output to file is slow because of overhead of each 
  //output call, so "output" each line to a stringstream, and when done 
  //with the line output the stringstream to the file so we have fewer 
  //write to file calls, also since I declared stringstream ss; in this 
  //function it disappears when the function exits so I can format ss 
  //instead of os so I don't have to worry about restoring the previous
  //format settings for os
    if (write_labels) {
      ss.setf(ios_base::left, ios_base::adjustfield);

      ss << '#';
      int correction = 1;      
      for(ivarsr=0; ivarsr<nvarsr; ++ivarsr) {
        ss << setw(surfpack::field_width - correction) << xrLabels[ivarsr] 
	   << " ";
        correction=0;
      }
      for(ivarsi=0; ivarsi<nvarsi; ++ivarsi) {
        ss << setw(surfpack::field_width - correction) << xiLabels[ivarsi]
	   << " ";
        correction=0;
      }
      for(iout=0; iout<nout-1; ++iout) {
        ss << setw(surfpack::field_width) << yLabels[iout] << " ";
      }
      ss << yLabels[nout-1];
      os << ss.str() << endl;
    }

  // ios_base::flags returns ios::fmtflags object, but OSF compiler doesn't 
  // like that.
  // Save the stream flags.  The output precision may be modified, but it 
  // will be restored to its old value before the method exits. 
    ss.precision(surfpack::output_precision);
    ss.setf(ios::scientific);
    //ss.width(surfpack::field_width);
    //ss.setw(surfpack::field_width);
    for(int ipt=0; ipt<npts; ++ipt) {
      ss.str("");
      for(ivarsr=0; ivarsr<nvarsr; ++ivarsr) 
	ss << setw(surfpack::field_width) << xr(ipt,ivarsr) << " ";
      for(ivarsi=0; ivarsi<nvarsi; ++ivarsi)
	ss << setw(surfpack::field_width) << xi(ipt,ivarsi) << " ";
      for(iout=0; iout<nout-1; ++iout)
	ss << setw(surfpack::field_width)<< y(ipt,iout) << " ";
      ss << y(ipt,nout-1);
      os << ss.str() << endl;
    }
}

/** Read a set of points into SurfData from an input stream of a binary file 
    "*.bspd". At this point nvarsr, nvarsi, nout, and optionally lockxr and 
    skip_columns must have already been set by the user. skip_columns 
    defaults to zero, need to do something better for variable labels (there should be a way to read them from the binary file)
*/
void SurfData::readBinary(istream& is, int skip_columns) 
{
  int npts_read = 0;
  try {
    is.read((char*)&npts  ,sizeof(npts));
    is.read((char*)&nvarsr,sizeof(nvarsr));
    is.read((char*)&nvarsi,sizeof(nvarsi));
    is.read((char*)&nout  ,sizeof(nout));
    
    xr.newSize(npts,nvarsr);
    xi.newSize(npts,nvarsi);
    y.newSize(npts,nout);

    defaultLabels(); //need to replace this with something better (there should be a way to store labels in the binary file)

    for (npts_read=0; npts_read < npts; ++npts_read) {
      // Throw an exception if we hit the end-of-file before we've
      // read the number of points that were supposed to be there.
      surfpack::checkForEOF(is);
      // True for fourth argument signals a binary read
      readPointBinary(npts_read, is, skip_columns);
    }

  } catch(surfpack::io_exception&) {
    cerr << "Expected: " << npts << " points.  "
         << "Read: " << npts_read << " points." << endl;
    throw;
  } 
}



///scale the data (output variable and real input variables) Note that this function should only be called by the Model if it is appropriate for the model to do so, it should not be called by anything other than a model. Also note that copies of a SurfData object will use the same default scaling of the original SurfData object, i.e if during construction you told the original SurfData object to use group scaling so will the copy
void SurfData::scaleToDefault()
{

  //WE WANT TO CHECK IF IT IS OK TO SCALE BEFORE SCALING
  int i, j, ifalreadyscaled=0;
  MtxDbl minmaxgroup(0,0);
  if(unscalexr.getNRows()!=0) {
    //if this is a non empty matrix
    assert((unscalexr.getNRows()==2)&&(unscalexr.getNCols()==nvarsr)&&
	   (unscaley.getNRows()==2)&&(unscaley.getNCols()==nout ));

    for(j=0; j<nvarsr; j++)
      if((fabs(unscalexr(0,j))!=1.0)||(unscalexr(1,j)!=0.0)) {
	ifalreadyscaled=1;
	break;
      }

    if(!ifalreadyscaled)
      for(j=0; j<nout; j++)
	if((fabs(unscaley(0,j))!=1.0)||(unscaley(1,j)!=0.0)) {
	  ifalreadyscaled=1;
	  break;
	}
  }
  else{
    printf("Warning: You just tried to scale() an empty surfdata object; ignoring request to scale!!!\n");
    assert(0);
    return;
  }

  if(ifalreadyscaled) unScale();
  
  //scale each output dimensions/variables individually
  //printf("scaleToDefault: before indivScale: y(0)=%g y(1)=%g",y(0),y(1));
  indivScale(y,unscaley,minmaxgroup,false);
  //printf(" after indivScale y(0)=%g y(1)=%g\n",y(0),y(1));
  if(lockxr.getNElems()) {
    //at this point we know that the user wants us to scale by groups and that the dimensions (in lockxr) have been sorted by groups, note that some groups may have only one dimension in them so individually scale the only-one-dimension-groups and only group-scale the groups with more than one dimension
    int igroupstart=0, igroupstop=0, ngroup,k;
    MtxInt igroup;    
    MtxDbl group;
    MtxDbl unscalegroup;

    for(j=1;j<nvarsr;j++) {
      if(lockxr(0,j)==lockxr(0,j-1))
	igroupstop=j;
      
      if((igroupstop==nvarsr-1)&&(igroupstart==0)) {
	groupScale(xr,unscalexr,minMaxXr,ifHaveMinMaxXr);
	break;
      }

      if(((lockxr(0,j)!=lockxr(0,j-1))||
	  (igroupstop==nvarsr-1))&&
	 (igroupstop!=igroupstart)) {
	//group scale if statement, we've found the end of a group is lock(0,j)!==lock(0,j-1) or we've reached the last dimension and the start of the group is not the same element as the end of the group (i.e. don't "group scale" a single dimension)
	ngroup=igroupstop-igroupstart+1;
	igroup.newSize(ngroup);
	for(k=0,i=igroupstart; i<=igroupstop; i++,k++)
	  igroup(k)=lockxr(1,i);	  
	xr.getCols(group,igroup);
	if(ifHaveMinMaxXr==true)
	  minMaxXr.getCols(minmaxgroup,igroup);
	groupScale(group,unscalegroup,minmaxgroup,ifHaveMinMaxXr);
	xr.putCols(group,igroup);
	unscalexr.putCols(unscalegroup,igroup);
	igroupstop=igroupstart=igroupstop+1;
      }

      if((igroupstop<j)||(igroupstop==nvarsr-1)) {
	//if igroup<j then we're at the end of a group and I didn't enter the "group scale" if statement so there's only a single dimension in this group so I want to individual scale; however, if the final dimension is the only one in it's group, and there was a "real" (more than 1 dimension) group that ended on the dimension before the final dimension, i.e. if igroupstop was nvarsr-2 and was just increased (i.e. igroupstop=igroupstart=igroupstop+1) then igroupstart=igroupstop=nvarsr-1 and in that case we want to individual scale the final (i.e. nvarsr-1) dimension even though the preceding (group scale) if statement was just executed.
	i=lockxr(1,igroupstart);
	xr.getCols(group,i);
	if(ifHaveMinMaxXr==true)
	  minMaxXr.getCols(minmaxgroup,i);
	indivScale(group,unscalegroup,minmaxgroup,ifHaveMinMaxXr);
	//indivScale(group,unscalegroup);
	xr.putCols(group,i);
	unscalexr.putCols(unscalegroup,i);
	igroupstop=igroupstart=igroupstop+1;
      }
	  
    } //yay!!! we're done scaling all the groups.
  }
  else{
    //printf("scaleToDefault::calling indivScale(xr,unscalexr)\n");
    //the user wanted us to individually scale (i.e not group scale) all the real input variables/dimensions
    indivScale(xr,unscalexr,minMaxXr,ifHaveMinMaxXr);
  }
  return;
}

  ///individually scale each dimension to length 1 centered at zero, (a "singular" dimension, i.e. a dimension in which all points share the same value, receives unscalea(0)=-1.0 and unscalea(1) is set to the single constant value, a negative unscalea(0) is a flag that says to exclude this dimension from the analysis) at the time of this writing this function is intended to only be called by the scaleToDefault() function
  void SurfData::indivScale(MtxDbl& a, MtxDbl& unscalea, const MtxDbl& minmaxa, bool if_have_minmaxa)
{
  int nrowsa=a.getNRows();
  int ncolsa=a.getNCols();
  unscalea.newSize(2,ncolsa);
  
  //printf("indivScale:: nrowsa=%d ncolsa=%d unscalea.getNRows()=%d unscalea.getNCols()=%d\n",nrowsa,ncolsa,unscalea.getNRows(),unscalea.getNCols()); fflush(stdout);
  

  //  unscalea.newsize(2,ncolsa);

  //printf("indivScale:: nrowsa=%d ncolsa=%d\n unscalea.getNRows()=%d unscalea.getNCols()=%d",nrowsa,ncolsa,unscalea.getNRows(),unscalea.getNCols());
  


  int i, j;
  double mina, maxa, scalea;

  for(j=0; j<ncolsa; j++) {
    if(if_have_minmaxa==true) {
      mina=minmaxa(0,j);
      maxa=minmaxa(1,j);
    }
    else 
      mina=maxa=a(0,j);
    
    for(i=0; i<nrowsa; i++) 
      if(a(i,j)<mina) mina=a(i,j);
      else if(a(i,j)>maxa) maxa=a(i,j);
  
    
    unscalea(0,j)=maxa-mina;
    unscalea(1,j)=0.5*(maxa+mina);

    if(unscalea(0,j)==0.0){
      //I refer to this as a "singular" dimension, all of points in this 
      //dimension are the same so exclude this dimension from anaylsis  
      unscalea(0,j)=-1.0; //a negative length scale is a flag to exclude this
      //dimension from the analysis; but we retain unscalea(1,j) to tell us 
      //what that constant value is
      for(i=0;i<nrowsa;i++)
	a(i,j)=0.0;
    }
    else{
      scalea=1.0/unscalea(0,j);
      for(i=0; i<npts; i++)
	a(i,j)=(a(i,j)-unscalea(1,j))*scalea;
    }
  }
  return;
};

  ///scale _this_ group of dimensions to a hyperrectangle of volume 1 centered at zero while preserving/locking the aspect ratio of the group of dimensions (note that this should only be called for __a__ __group__ (as in 1 group at a time) of real input variables, (a "singular" dimension j, i.e. a dimension in which all points share the same value, is not counted as part of the group, instead it is mapped to 0.0, has unscalea(0,j)=-1.0 and unscalea(1,j) is set to the single constant value) at the time of this writing this function is intended to only be called by the scaleToDefault() function
void SurfData::groupScale(MtxDbl& a, MtxDbl& unscalea, const MtxDbl& minmaxa, bool if_have_minmaxa)
{  
  int nrowsa=a.getNRows();
  int ncolsa=a.getNCols();
  unscalea.newSize(2,ncolsa);

  int i, j;
  double mina, maxa, scalea, unscale_all;
  int numdiffa=0;
  double proddiffa=1.0;
  
  for(j=0; j<ncolsa; j++) {
    if(if_have_minmaxa==true) {
      mina=minmaxa(0,j);
      maxa=minmaxa(1,j);
    }
    else 
      mina=maxa=a(0,j);
    
    for(i=0; i<nrowsa; i++) 
      if(a(i,j)<mina) mina=a(i,j);
      else if(a(i,j)>maxa) maxa=a(i,j);
    

    unscalea(0,j)=maxa-mina;
    unscalea(1,j)=0.5*(maxa+mina);

    if(unscalea(0,j)==0.0){
      //I refer to this as a "singular" dimension, all of points in this 
      //dimension are the same so exclude this dimension from anaylsis  
      unscalea(0,j)=-1.0; //a negative length scale is a flag to exclude this
      //dimension from the analysis; but we retain unscalea(1,j) to tell us 
      //what that constant value is
      for(i=0;i<nrowsa;i++)
	a(i,j)=0.0;
    }
    else{
      numdiffa++;
      proddiffa*=unscalea(0,j);
    }
  }

  unscale_all=pow(proddiffa,1.0/numdiffa);
  scalea=1.0/unscale_all;

  for(j=0; j<ncolsa; j++) 
    if(unscalea(0,j)!=-1.0) {
      //exclude dimensions in which all points are the same from analysis
      unscalea(0,j)=unscale_all;
      for(i=0; i<nrowsa; i++)
	a(i,j)=(a(i,j)-unscalea(1,j))*scalea;
    }
  return;
}

///scale real input variables to domain_new(0,j)<=xr(:,j)<=domain_new(1,j) for j=0,1,...,nvarsr;  For "singular" dimension j xr(:,j)=0.5*(domain_new(0,j)+domain_new(1,j)),unscalexr(0,j)=-1.0 and unscalexr(1,j) = the single value - xr(:,j)
void SurfData::scaleXrToDomain(MtxDbl& domain_new){
  int i;
  double scaleby, shiftby, minxr, maxxr, temp;  
  for(int j=0;j<nvarsr;j++) {
      temp=domain_new(1,j)-domain_new(0,j);
      assert(temp!=0.0);	
      //scaleby=unscalexr(0,j); scaleby=fabs(scaleby);
      scaleby=fabs(this->unscalexr(0,j));
      shiftby=unscalexr(1,j);
      //unscale xr and determine the unscaled xr's min and max values
      if(ifHaveMinMaxXr==true) {
	minxr=minMaxXr(0,j)*scaleby+shiftby;
	maxxr=minMaxXr(1,j)*scaleby+shiftby;
	if(!((scaleby==1.0)&&(shiftby=0.0)))
	  for(i=1; i<npts; ++i)
	    xr(i,j)=xr(i,j)*scaleby+shiftby;
      }
      else {
	minxr=maxxr=xr(0,j)=xr(0,j)*scaleby+shiftby;
	for(i=1; i<npts; i++) {
	  xr(i,j)=xr(i,j)*scaleby+shiftby;
	  if(xr(i,j)<minxr) minxr=xr(i,j);
	  if(xr(i,j)>maxxr) maxxr=xr(i,j);
	}
      }
      if(maxxr==minxr){
	unscalexr(0,j)=-1.0; //can't define any length scale 
	shiftby=0.5*(domain_new(1,j)+domain_new(0,j));
	unscalexr(1,j)=maxxr-shiftby;	
	for(i=0;i<npts;i++)
	  xr(i,j)=shiftby;
      }
      else{
	unscalexr(0,j)=(maxxr-minxr)/temp;
	unscalexr(1,j)=minxr-unscalexr(0,j)*domain_new(0,j);
	scaleby=1.0/unscalexr(0,j);
	shiftby=-unscalexr(1,j);
	for(i=0; i<npts; i++) 
	  xr(i,j)=(xr(i,j)+shiftby)*scaleby;
      }
  }
  return;
}

void SurfData::scaleXrToFactor(MtxDbl& unscale_xr) {
  assert((unscale_xr.getNRows()==2)&&(unscale_xr.getNCols()==nvarsr));
  int i, j;
  double scaleby, shiftby;
  for(j=0; j<nvarsr; j++) {
    assert(unscale_xr(0,j));
    scaleby=fabs(unscalexr(0,j))/fabs(unscale_xr(0,j));
    shiftby=(unscalexr(1,j)-unscale_xr(1,j))/fabs(unscale_xr(0,j));
    unscalexr(0,j)=unscale_xr(0,j);
    unscalexr(1,j)=unscale_xr(1,j);
    for(i=0; i<npts; i++)
      xr(i,j)=xr(i,j)*scaleby+shiftby;
  }
  return;
}

void SurfData::scaleYToFactor(MtxDbl& unscale_y) {
  assert((unscale_y.getNRows()==2)&&(unscale_y.getNCols()==nout));
  int i, j;
  double scaleby, shiftby;
  for(j=0; j<nout; j++) {
    assert(unscale_y(0,j));
    scaleby=fabs(unscaley(0,j))/fabs(unscale_y(0,j));
    shiftby=(unscaley(1,j)-unscale_y(1,j))/fabs(unscale_y(0,j));
    unscaley(0,j)=unscale_y(0,j);
    unscaley(1,j)=unscale_y(1,j);
    for(i=0; i<npts; i++)
      y(i,j)=y(i,j)*scaleby+shiftby;
  }
  return;
}

///unscale this surfdata object 
SurfData& SurfData::unScale()
{
  //unscale the real input variables
  int i, j;
  double scaleby, shiftby;
  for(j=0; j<nvarsr; j++) {
    scaleby=fabs(unscalexr(0,j)); //in case this dimension was singular
    shiftby=unscalexr(1,j);
    unscalexr(0,j)=1.0;
    unscalexr(1,j)=0.0;
    for(i=0; i<npts; i++)
      xr(i,j)=xr(i,j)*scaleby+shiftby;
  }
  
  //unscale the output variables  
  for(j=0; j<nout; j++) {
    scaleby=fabs(unscaley(0,j)); //in case this dimension was singular
    shiftby=unscaley(1,j);
    unscaley(0,j)=1.0;
    unscaley(1,j)=0.0;
    for(i=0; i<npts; i++)
      y(i,j)=y(i,j)*scaleby+shiftby;
  }
  
  return *this;
}

MtxDbl& SurfData::scaleYOther(MtxDbl& y_other, int j)
{
  if(j==-99999) j=jout;
#ifdef __SURFDATA_ERR_CHECK__
  assert((0<=j)&&(j<nout));
#endif
  int npts_other=y_other.getNRows();
  int nout_other=y_other.getNCols();
  double scaleby, shiftby;
  if(nout_other==1) {
    scaleby=1.0/fabs(unscaley(0,j));
    shiftby=unscaley(1,j);
    for(int i=0; i<npts_other; ++i)
      y_other(i)=(y_other(i)-shiftby)*scaleby;
  }
  else if(nout_other==nout) {
    for(j=0; j<nout; ++j) {
      scaleby=1.0/fabs(unscaley(0,j));
      shiftby=unscaley(1,j);
      for(int i=0; i<npts_other; ++i)
	y_other(i,j)=(y_other(i,j)-shiftby)*scaleby;
    }
  }
  else {
    printf("MtxDbl& SurfData::scaleYOther(MtxDbl& y_other, int j=jout)... nout=%d & nout_other=%d but must equal 1 or nout\n",nout,nout_other);
#ifdef __SURFDATA_ERR_CHECK__
    fflush(stdin);
    assert((nout_other==1)||(nout_other==nout));
#endif
  }
  return y_other;
}
  
MtxDbl& SurfData::unScaleYOther(MtxDbl& y_other, int j)
{
  if(j==-99999) j=jout;
#ifdef __SURFDATA_ERR_CHECK__
  assert((0<=j)&&(j<nout));
#endif
  int npts_other=y_other.getNRows();
  int nout_other=y_other.getNCols();
  double scaleby, shiftby;
  if(nout_other==1) {
    scaleby=fabs(unscaley(0,j));
    shiftby=unscaley(1,j);
    for(int i=0; i<npts_other; ++i)
      y_other(i)=y_other(i)*scaleby+shiftby;
  }
  else if(nout_other==nout) {
    for(j=0; j<nout; ++j) {
      scaleby=fabs(unscaley(0,j));
      shiftby=unscaley(1,j);
      for(int i=0; i<npts_other; ++i)
	y_other(i,j)=y_other(i,j)*scaleby+shiftby;
    }
  }
  else{
    printf("MtxDbl& SurfData::unScaleYOther(MtxDbl& y_other, int j=jout)... nout=%d & nout_other=%d but must equal 1 or nout\n",nout,nout_other);      
#ifdef __SURFDATA_ERR_CHECK__
    fflush(stdin);
    assert((nout_other==1)||(nout_other==nout));
#endif 
  }
  return y_other;
}




///this function should only be called by putPoints() it decides whether to recommend that this SurfData object be rescaled after newpoints2 are added to the current point and returns 0 or 1 to indicate the recommendation; newpoints2 must be non-empty and already have the same scale as this SurfData object.  * 0 means we do not recommend rescaling: 0 is returned if this SurfData object is unscaled, or if the all new points fall within this SurfData object's range of real inputs and outputs. * 1 means that we recommend rescaling: we recommend rescaling if this SurfData object is scaled and at least one of the points in newpoints2 is outside this SurfData object's range of real inputs and outputs.
int SurfData::ifRecommendRescale(SurfData& newpoints2) {
  int nnew=newpoints2.npts;
  int i,j;
  //we need to determine if we should recommend rescaling to the 
  //user/programmer
  int ifrecommendrescale;
  //First thing to check: if the current data is unscaled we will 
  //recommend to NOT rescale the data after adding the new points
  
  int ifunscaled=1;
  for(i=0;i<nvarsr;i++)
    if((unscalexr(0,i)!=1.0)||(unscalexr(1,i)!=0.0)){
      ifunscaled=0;
      break;
    }
  if(ifunscaled)
    for(i=0;i<nout;i++)
      if((unscalexr(0,i)!=1.0)||(unscalexr(1,i)!=0.0)){
	ifunscaled=0;
	break;
      }
  if(ifunscaled) {
    //the current data is unscaled so we recommend to NOT rescale the data
    ifrecommendrescale=0;
  }
  else{
    ifrecommendrescale=0; //by default we recommend not to rescale
    //we need to check if any of the new points are outside the range of the 
    //current points, if so we will recommend to rescale if not we will 
    //recommend to NOT rescale
    MtxDbl minmaxold(2), minmaxnew(2);
    //first we check xr
    for(j=0; j<nvarsr; j++) {
      minmaxold(0)=minmaxold(1)=xr(0,j);
      for(i=1;i<npts;i++) {
	if(xr(i,j)<minmaxold(0)) minmaxold(0)=xr(i,j);
	if(xr(i,j)>minmaxold(1)) minmaxold(1)=xr(i,j);	       
      }
      minmaxnew(0)=minmaxnew(1)=newpoints2.xr(0,j);
      for(i=1;i<nnew;i++){
	if(newpoints2.xr(i,j)<minmaxnew(0)) minmaxnew(0)=newpoints2.xr(i,j);
	if(newpoints2.xr(i,j)>minmaxnew(1)) minmaxnew(1)=newpoints2.xr(i,j);	       
      }
      if((minmaxnew(0)<minmaxold(0))||(minmaxold(1)<minmaxnew(1))){
	ifrecommendrescale=1;
	break;
      }
    }
    if(!ifrecommendrescale) {
      //now we check y
      for(j=0; j<nout; j++) {
	minmaxold(0)=minmaxold(1)=y(0,j);
	for(i=1;i<npts;i++) {
	  if(y(i,j)<minmaxold(0)) minmaxold(0)=y(i,j);
	  if(y(i,j)>minmaxold(1)) minmaxold(1)=y(i,j);	       
	}
	minmaxnew(0)=minmaxnew(1)=newpoints2.y(0,j);
	for(i=1;i<nnew;i++){
	  if(newpoints2.y(i,j)<minmaxnew(0)) minmaxnew(0)=newpoints2.y(i,j);
	  if(newpoints2.y(i,j)>minmaxnew(1)) minmaxnew(1)=newpoints2.y(i,j);	       
	}
	if((minmaxnew(0)<minmaxold(0))||(minmaxold(1)<minmaxnew(1))){
	  ifrecommendrescale=1;
	  break;
	}
      } 
    }
  }
  return ifrecommendrescale;
}

///int putPoints(SurfData& newpoints, int ipt) puts the single point in newpoints into this SurfData objects ipt-th point, ipt must be greater than or equal to zero, if ipt is less than the current number of points the previous ipt-th point is replaced with the new one, if ipt is greater than the current number of points this function will enlarge the matrices to hold ipt+1 points before inserting this point, if ipt is unspecified it appends the newpoint to the end of this SurfData object's current list of points (after enlarging the matrices), if newpoints contains multiple points and ipt is unspecified it appends all the newpoints to the end of this SurfData obeject's current list of points (if ipt is specified newpoints must contain exactly one point, otherwise you should call int putPoints(SurfData& newpoints, MtxInt& ipts) instead).  If this SurfData object contained any points before this function was called then the newpoints are automatically rescaled to match the current points before they are inserted.  putPoints returns an integer to recommend whether or not the user/programmer should scale this SurfData object after calling this function: * 0 means we do not recommend rescaling: 0 is returned if newpoints was empty, this SurfData object was unscaled before this function was called, or if the all new points fell within the range of real inputs and outputs that this SurfData object had before this function was called. * -1 means that this SurfData object was empty before this function was called so we kept the scaling in newpoints, * 1 means that we recommend rescaling: we recommend rescaling if this SurfData object was scaled before this function was called and at least one of the newpoints real inputs or outputs was outside the range of range of the old  points.  This function does not actually do the rescaling of this SurfData object itself in case the user/programmer specified had scaled to a factor or scaled to a domain, or if for the sake of consistency/comparison-to-the-values-before-the-new-datapoints-were-added he/she doesn't want to rescale.
int SurfData::putPoints(SurfData& newpoints, int ipt) {
  int nnew=newpoints.npts;
  if(nnew==0) {
    cerr << "Warning!!! in: 'int SurfData::putPoints(SurfData& newpoints, int ipt)' newpoints is empty." << endl;
    return 0; //recommend to NOT rescale this SurfData object because we didn't change it
  }

  if(npts>0)
    assert((nvarsr==newpoints.nvarsr)&&
	   (nvarsi==newpoints.nvarsi)&&
	   (nout  ==newpoints.nout  ));
  else{
    if(ipt==-99999){
      *this=newpoints;
      return -1; //we kept the scaling in newpoints because we didn't have 
      //any points in this SurfData object previously
    }
    assert((nnew==1)&&(0<=ipt));
    npts=ipt+1;
    nvarsr   =newpoints.nvarsr;
    nvarsi   =newpoints.nvarsi;
    nout     =newpoints.nout;
    jout     =newpoints.jout;
    unscalexr=newpoints.unscalexr;
    unscaley =newpoints.unscaley;
    lockxr   =newpoints.lockxr;
    xrLabels =newpoints.xrLabels;
    xiLabels =newpoints.xiLabels;
    yLabels  =newpoints.yLabels;
    xr.newSize(npts,nvarsr); xr.putRows(newpoints.xr,ipt);
    xi.newSize(npts,nvarsi); xi.putRows(newpoints.xi,ipt);
    y.newSize( npts,nout  ); y.putRows( newpoints.y ,ipt);
    return -1;  //we kept the scaling in newpoints because we didn't have 
    //any points in this SurfData object previously
  }

  //we already have some points so we need to make sure the newpoints
  //have the same number of each type of variable as the current points
  assert((nvarsr==newpoints.nvarsr)&&
	 (nvarsi==newpoints.nvarsi)&&
	 (nout  ==newpoints.nout  ));
  //convert the newpoints to the same scale as the current points
  SurfData newpoints2=newpoints; //don't want to change this in case the
  //user is still interested in using it for something else
  newpoints2.scaleToFactors(unscalexr,unscaley);
  
  //decide whether or not to recommend that this SurfData object after 
  //the newpoints have been added to the current points
  int ifrecommendrescale=ifRecommendRescale(newpoints2);

  if(ipt==-99999) {
    //the flag to append at the end so we need to enlarge our matrices
    //this is triggered in by NOT passing in a valued for ipt, -99999
    //is a default value
    ipt=npts; //if it's a single point, set ipt to put it at the end
    npts=npts+nnew;
    xr.resize(npts,nvarsr);
    xi.resize(npts,nvarsi);
    y.resize( npts,nout  );
    if(nnew==1){
      //there's just one point so don't do anything special
      xr.putRows(newpoints2.xr,ipt);
      xi.putRows(newpoints2.xi,ipt);
      y.putRows( newpoints2.y ,ipt);
    }
    else{
      //there's multiple new points, so we need to create an integer 
      //vector of indices to specify the rows of the matrices we want
      //to fill in
      MtxInt ipts(nnew);  
      for(int i=0; i<nnew; i++)
	ipts(i)=ipt+i;
      xr.putRows(newpoints2.xr,ipts);
      xi.putRows(newpoints2.xi,ipts);
      y.putRows( newpoints2.y ,ipts);
    }
  }
  else{
    //the user is specifying a particular row to fill in
    
    //make sure there is only one point, otherwise he would have needed 
    //to pass in a MtxInt& ipts containing the indices, and that would have
    //called the other (overloaded) function with the same name, i.e. 
    //void putPoints(SurfData& newpoints, MtxInt& ipts)
    assert((nnew==1)&&(0<=ipt));
    if(npts<=ipt) {
      //the user is specifying a point beyond the edge of our matrices so 
      //we need to enlarge them before we fill in the new point.
      //if the specified point to "replace" is "greater than" (or equal to) 
      //the number of points (i.e. the size of the matrices) we currently 
      //have, enlarge the matrices
      npts=npts+nnew;
      xr.resize(npts,nvarsr);
      xi.resize(npts,nvarsi);
      y.resize( npts,nout  );
    }
    //fill in the new point
    xr.putRows(newpoints2.xr,ipt);
    xi.putRows(newpoints2.xi,ipt);
    y.putRows( newpoints2.y ,ipt);
  }

  return ifrecommendrescale; //make the recommendation about whether or
  //not to rescale
}


///int putPoints(SurfData& newpoints, MtxInt& ipts) puts the points in newpoints into the points in this SurfData object listed in ipts, all values contained in ipts must be greater than or equal to zero, if any point listed in ipts is less than the current number of points the previous point is replaced with the new one, if any point listed in ipts is greater than the current number of points this function will enlarge the matrices to hold ipts.max()+1 points before inserting these points, if ipts is unspecified then int putPoints(SurfData& newpoints, int ipt) is called instead (with ipt defaulting to -99999 which is the flag to append the new points to the current list of points).  If this SurfData object contained any points before this function was called then the newpoints are automatically rescaled to match the current points before they are inserted. putPoints returns an integer to recommend whether or not the user/programmer should scale this SurfData object after calling this function: * 0 means we do not recommend rescaling: 0 is returned if newpoints was empty, this SurfData object was unscaled before this function was called, or if the all new points fell within the range of real inputs and outputs that this SurfData object had before this function was called. * -1 means that this SurfData object was empty before this function was called so we kept the scaling in newpoints, * 1 means that we recommend rescaling: we recommend rescaling if this SurfData object was scaled before this function was called and at least one of the newpoints real inputs or outputs was outside the range of range of the old  points.  This function does not actually do the rescaling of this SurfData object itself in case the user/programmer specified had scaled to a factor or scaled to a domain, or if for the sake of consistency/comparison-to-the-values-before-the-new-datapoints-were-added he/she doesn't want to rescale.
int SurfData::putPoints(SurfData& newpoints, MtxInt& ipts) {
  int nnew=newpoints.npts;
  if(nnew==0) {
    cerr << "Warning!!! in: 'int SurfData::putPoints(SurfData& newpoints, MtxInt& ipts)' newpoints is empty." << endl;
    return 0; //recommend to NOT rescale this SurfData object because we didn't change it
  }
  assert((0<=ipts.minElem())&&
	 (newpoints.npts==ipts.getNElems()));
  int ifrecommendrescale;
  SurfData newpoints2=newpoints; 
  if(npts==0) {
    nvarsr   =newpoints.nvarsr;
    nvarsi   =newpoints.nvarsi;
    nout     =newpoints.nout;
    jout     =newpoints.jout;
    unscalexr=newpoints.unscalexr;
    unscaley =newpoints.unscaley;
    lockxr   =newpoints.lockxr;
    xrLabels =newpoints.xrLabels;
    xiLabels =newpoints.xiLabels;
    yLabels  =newpoints.yLabels;
    ifrecommendrescale=-1;
  }
  else{
    //check to make sure that we have valid input
    assert((nvarsr==newpoints.nvarsr)&&
	   (nvarsi==newpoints.nvarsi)&&
	   (nout  ==newpoints.nout  ));

    newpoints2.scaleToFactors(unscalexr,unscaley);

    //decide whether or not to recommend that this SurfData object after 
    //the newpoints have been added to the current points
    int ifrecommendrescale=ifRecommendRescale(newpoints2);
  }

  int iptsmax=ipts.maxElem();
  if(npts<=iptsmax){
    npts=iptsmax+1;
    xr.resize(npts,nvarsr);
    xi.resize(npts,nvarsi);
    y.resize( npts,nout  );
  }
  xr.putRows(newpoints2.xr,ipts);
  xi.putRows(newpoints2.xi,ipts);
  y.putRows( newpoints2.y ,ipts);

  return ifrecommendrescale;
}

///retrieve one point (the one with index ipt), return it as a reference to SurfData
SurfData& SurfData::getPoints(SurfData& result, int ipt) {
  assert((0<=ipt)&&(ipt<npts));
  result.npts     =1;
  result.nvarsr   =nvarsr;
  result.nvarsi   =nvarsi;
  result.nout     =nout;
  result.jout     =jout;
  result.unscalexr=unscalexr;
  result.unscaley =unscaley;
  result.lockxr   =lockxr;
  result.xrLabels =xrLabels;
  result.xiLabels =xiLabels;
  result.yLabels  =yLabels;

  xr.getRows(result.xr,ipt);
  xi.getRows(result.xi,ipt);
  y.getRows( result.y ,ipt);

  return result;
}

///retrieve multiple points (the ones whose indices are stored in ipts), pass them back as a reference to a SurfData object
SurfData& SurfData::getPoints(SurfData& result, MtxInt& ipts) {
  int nget=ipts.getNElems();
  ipts.uniqueElems();
  assert(nget==ipts.getNElems());
  assert((0<=ipts(0))&&(ipts(nget-1)<npts));
  result.npts     =nget;
  result.nvarsr   =nvarsr;
  result.nvarsi   =nvarsi;
  result.nout     =nout;
  result.jout     =jout;
  result.unscalexr=unscalexr;
  result.unscaley =unscaley;
  result.lockxr   =lockxr;
  result.xrLabels =xrLabels;
  result.xiLabels =xiLabels;
  result.yLabels  =yLabels;

  xr.getRows(result.xr,ipts);
  xi.getRows(result.xi,ipts);
  y.getRows( result.y ,ipts);

  return result;
}

///retrieve all points except one (the one with index ipt), return them as a reference to a SurfData object
SurfData& SurfData::excludePoints(SurfData& result, int ipt) {
  assert((0<=ipt)&&(ipt<npts));
  result.npts     =npts-1;
  result.nvarsr   =nvarsr;
  result.nvarsi   =nvarsi;
  result.nout     =nout;
  result.jout     =jout;
  result.unscalexr=unscalexr;
  result.unscaley =unscaley;
  result.lockxr   =lockxr;
  result.xrLabels =xrLabels;
  result.xiLabels =xiLabels;
  result.yLabels  =yLabels;

  xr.excludeRows(result.xr,ipt);
  xi.excludeRows(result.xi,ipt);
  y.excludeRows( result.y ,ipt);

  return result; 
}


///retrieve all points except the ones whose indices are listed in ipts, return them as a reference to a SurfData object
SurfData& SurfData::excludePoints(SurfData& result, MtxInt& ipts) {
  int nexclude=ipts.getNElems();
  ipts.uniqueElems();
  assert(nexclude==ipts.getNElems());
  assert((0<=ipts(0))&&(ipts(nexclude-1)<npts));
  result.npts     =npts-nexclude;
  result.nvarsr   =nvarsr;
  result.nvarsi   =nvarsi;
  result.nout     =nout;
  result.jout     =jout;
  result.unscalexr=unscalexr;
  result.unscaley =unscaley;
  result.lockxr   =lockxr;
  result.xrLabels =xrLabels;
  result.xiLabels =xiLabels;
  result.yLabels  =yLabels;

  xr.excludeRows(result.xr,ipts);
  xi.excludeRows(result.xi,ipts);
  y.excludeRows( result.y ,ipts);

  return result;
}

///create a split copy of the points stored in this surfdata, one point, the one whose index is listed in ipt, will be placed in "extracted", the rest will be placed in "rest", this is very useful for the PRESS metric
void SurfData::extractPoints(SurfData& rest, SurfData& extracted, int ipt) {
  assert((0<=ipt)&&(ipt<npts));
  /*
  assert((npts  ==xr.getNRows())&&
	 (npts  ==y.getNRows())&&
	 (nvarsr==xr.getNCols())&&
	 (nout  ==y.getNCols()));  
  printf("extractPoints(): npts=%d nvarsr=%d nout=%d\n",npts,nvarsr,nout); fflush(stdout);
  */

  extracted.npts     =1;
  extracted.nvarsr   =nvarsr;
  extracted.nvarsi   =nvarsi;
  extracted.nout     =nout;
  extracted.jout     =jout;
  extracted.unscalexr=unscalexr;
  //for(int i=0; i<nvarsr; ++i)
  //printf("unscalexr(:,%d)=[%g;%g] ",i,unscalexr(0,i), unscalexr(1,i));
  //printf("\n");
  extracted.unscaley =unscaley;
  //for(int i=0; i<nout; ++i)
  //printf("unscaley(:,%d)=[%g;%g] ",i,unscaley(0,i), unscaley(1,i));
  //printf("\n");
  extracted.lockxr   =lockxr;
  extracted.xrLabels =xrLabels;
  extracted.xiLabels =xiLabels;
  extracted.yLabels  =yLabels;

  /*
  assert((npts ==xr.getNRows())&&
	 (npts ==y.getNRows())&&
	 (nvarsr==xr.getNCols())&&
	 (nout ==y.getNCols()));  
  printf("extractPoints(): npts=%d nvarsr=%d nout=%d\n",npts,nvarsr,nout); fflush(stdout);
  */

  xr.getRows(extracted.xr,ipt);
  if(nvarsi) xi.getRows(extracted.xi,ipt);
  y.getRows(extracted.y,ipt);

  rest.npts     =npts-1;
  rest.nvarsr   =nvarsr;
  rest.nvarsi   =nvarsi;
  rest.nout     =nout;
  rest.jout     =jout;
  rest.unscalexr=unscalexr;
  rest.unscaley =unscaley;
  rest.lockxr   =lockxr;
  rest.xrLabels =xrLabels;
  rest.xiLabels =xiLabels;
  rest.yLabels  =yLabels;

  xr.excludeRows(rest.xr,ipt);
  if(nvarsi) xi.excludeRows(rest.xi,ipt);
  y.excludeRows(rest.y,ipt);

  return;
}

///create a split copy of the points stored in this surfdata, the ones whose indices are listed in ipts will be placed in "extracted", the rest will be placed in "rest", this is very useful for the Leave N Out cross Validation metric
void SurfData::extractPoints(SurfData& rest, SurfData& extracted, MtxInt& ipts) {
  int nextract=ipts.getNElems();
  ipts.uniqueElems();
  assert(nextract==ipts.getNElems());
  assert((0<=ipts(0))&&(ipts(nextract-1)<npts));
  extracted.npts     =nextract;
  extracted.nvarsr   =nvarsr;
  extracted.nvarsi   =nvarsi;
  extracted.nout     =nout;
  extracted.jout     =jout;
  extracted.unscalexr=unscalexr;
  extracted.unscaley =unscaley;
  extracted.lockxr   =lockxr;
  extracted.xrLabels =xrLabels;
  extracted.xiLabels =xiLabels;
  extracted.yLabels  =yLabels;

  xr.getRows(extracted.xr,ipts);
  if(nvarsi) xi.getRows(extracted.xi,ipts);
  y.getRows( extracted.y ,ipts);

  rest.npts     =npts-nextract;
  rest.nvarsr   =nvarsr;
  rest.nvarsi   =nvarsi;
  rest.nout     =nout;
  rest.jout     =jout;
  rest.unscalexr=unscalexr;
  rest.unscaley =unscaley;
  rest.lockxr   =lockxr;
  rest.xrLabels =xrLabels;
  rest.xiLabels =xiLabels;
  rest.yLabels  =yLabels;

  xr.excludeRows(rest.xr,ipts);
  if(nvarsi) xi.excludeRows(rest.xi,ipts);
  y.excludeRows( rest.y ,ipts);

  return;
}


/***********************************************************************/
/***********************************************************************/
/**** SurfData member functions end here                            ****/
/***********************************************************************/
/***********************************************************************/

} // end namespace nkm
