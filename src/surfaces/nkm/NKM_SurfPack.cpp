#include <math.h>
#include "NKM_SurfPack.hpp"

//the purpose of this file is to contain generic functions usable by everything

//#define __SURFPACK_ERR_CHECK__
//#define __SURFPACK_UNIT_TEST__
//#define __SURFPACK_DER_TEST__

namespace nkm {

//not sure which of these usings are necessary, achieve compilation by random guessing
using std::cerr;
using std::endl;
using std::ifstream;
using std::istream;
using std::ios;
//using std::numeric_limits;
using std::ofstream;
using std::ostream;
using std::setw;
using std::string;
//using std::vector;

#ifdef __SURFPACK_UNIT_TEST__
int main(){
  int i, j, k;
  for(int n=0; n<5; n++)
    for(k=0; k<=n; k++)
      printf("nchoosek(%d,%d)=%d\n",n,k,nchoosek(n,k));

  MtxInt poly; 
  multi_dim_poly_power(poly,4,4);
  int npoly=poly.getNRows();
  int nvarsr=poly.getNCols();
  printf("poly.getNRows()=%d\npoly=...\n",poly.getNRows());
  
  for(i=0; i<npoly; i++) {
    printf("%8d [",i);
    for(j=0; j<nvarsr; j++)
      printf(" %d",poly(i,j));
    printf(" ]\n");
  }
  
  double pi=acos(0.0)*2.0;
  MtxDbl EulAng(6), Rot(4,4), RotShouldBe(4,4); 
  EulAng(0)=pi/6.0;
  EulAng(1)=pi/4.0;
  EulAng(2)=pi/3.0;
  EulAng(3)=pi/6.0;
  EulAng(4)=pi/4.0;
  EulAng(5)=pi/3.0;
  //This RotShouldBe was evaluated for these EulAng by matlab code that I know works
  RotShouldBe(0,0)= -0.057800215120219;
  RotShouldBe(1,0)=  0.353765877365274;
  RotShouldBe(2,0)=  0.768283046242747;
  RotShouldBe(3,0)=  0.530330085889911;
  RotShouldBe(0,1)= -0.695272228311384;
  RotShouldBe(1,1)= -0.649306566066328;
  RotShouldBe(2,1)=  0.035320133098213;
  RotShouldBe(3,1)=  0.306186217847897;
  RotShouldBe(0,2)=  0.647692568794007;
  RotShouldBe(1,2)= -0.414729655649473;
  RotShouldBe(2,2)= -0.183012701892219;
  RotShouldBe(3,2)=  0.612372435695795;
  RotShouldBe(0,3)= -0.306186217847897;
  RotShouldBe(1,3)=  0.530330085889911;
  RotShouldBe(2,3)= -0.612372435695795;
  RotShouldBe(3,3)=  0.500000000000000;
  gen_rot_mat(Rot,EulAng,4);
  printf("Rot=...\n");
  for(i=0; i<nvarsr; i++){
    printf("    [ ");
    for(j=0; j<nvarsr; j++) {
      printf("%9f ",Rot(i,j));
      if(fabs(Rot(i,j)-RotShouldBe(i,j))>1.0e-15) {
	printf("\n(%d,%d):Rot=%.16f RotShouldBe=%.16f\n",i,j,Rot(i,j),RotShouldBe(i,j));
	fflush(stdout);
	assert(Rot(i,j)==RotShouldBe(i,j));
      }
    }
    printf("]\n");
  }
  printf("We know Rot is ok because it didn't assert(false)\n");

  MtxDbl should_be_identity(4,4);
  matrix_mult(should_be_identity,Rot,Rot,0.0,1.0,'N','T');
  printf("should_be_identity=Rot*Rot^T=...\n");
  for(i=0; i<nvarsr; i++){
    printf("    [ ");
    for(j=0; j<nvarsr; j++) {
      printf("%9f ",should_be_identity(i,j));
    }
    printf("]\n");
  }

  matrix_mult(should_be_identity,Rot,Rot,0.0,1.0,'T','N');
  printf("should_be_identity=Rot^T*Rot=...\n");
  for(i=0; i<nvarsr; i++){
    printf("    [ ");
    for(j=0; j<nvarsr; j++) {
      printf("%9f ",should_be_identity(i,j));
    }
    printf("]\n");
  }


  int rotnvarsr=5;
  MtxDbl BaseAxis(2*rotnvarsr,rotnvarsr); BaseAxis.zero();
  for(i=0; i<rotnvarsr; i++) {
    BaseAxis(i,i)=1.0;
    BaseAxis(i+rotnvarsr,i)=-1.0;
  }
  MtxDbl Axis;

  int ioct;
  int noct=(int) pow(2,rotnvarsr);
  MtxInt InOct(noct); InOct.zero();
  int NDV=(rotnvarsr*(rotnvarsr-1))/2; //nchoosek(nvarsr,2);
  MtxDbl lowerBounds(1,NDV); lowerBounds.zero();
  MtxDbl upperBounds(1,NDV);
  int mymod=1048576; //2^20 instead of 10^6 to be kind to the computer
  for(j=0;j<NDV; j++) upperBounds(j)=pi;

  int nguess=2*rotnvarsr*noct*100;

  EulAng.newSize(NDV);
  for(i=0; i<nguess; i++) {
    for(j=0;j<NDV;j++)
      EulAng(j)=(rand()%mymod)*upperBounds(j)/mymod;
    gen_rot_mat(Rot,EulAng,rotnvarsr);
    matrix_mult(Axis,BaseAxis,Rot,0.0,1.0);
    for(k=0;k<2*rotnvarsr;k++){
      //{k=0;
      ioct=(Axis(k,0)>=0.0);
      for(j=1;j<rotnvarsr;j++)
	ioct+=(Axis(k,j)>=0)*((int) pow(2.0,j));
      InOct(ioct)++;
    }
  }
  printf("relative # of axis endpoints per octant\n");
  for(i=0;i<noct;i++)
    printf("  InOct(%d/%d)=%g\n",i,noct,InOct(i)/(2.0*(nguess*rotnvarsr/noct)));

  return 0;


}
#endif

int if_close_enough(double a, double b)
{
  if(fabs(a-b)>1.0e-5){
    printf("a=%20.14f b=%20.14f\n",a,b);  fflush(stdout);
    assert((fabs(a-b)<=1.0e-5));
  }
  return (fabs(a-b)<=1.0e-5);
};

///this function should be moved to surfpack.cpp
int nchoosek(int n, int k){
  int nck=1;

  int nloop=(k<n-k)?k:n-k;
  if(nloop>0) {
    nck=n;
    for(int iloop=1; iloop<nloop; iloop++)
      nck=(nck*(n-iloop))/(iloop+1);
  }
    
  return nck;
}

///this function should be moved to surfpack.cpp
MtxDbl& gen_rot_mat(MtxDbl& Rot, const MtxDbl& EulAng, int nvarsr){
#ifdef __SURFPACK_ERR_CHECK__
  assert(EulAng.getNElems()==(nvarsr*(nvarsr-1)/2));
#endif
  MtxDbl I(nvarsr,nvarsr), R(nvarsr,nvarsr), RotTemp(nvarsr,nvarsr);
  I.zero();
  int ivarr;
  for(ivarr=0; ivarr<nvarsr; ivarr++) 
    I(ivarr,ivarr)=1.0;
  Rot=I;
  int nang=nvarsr, iang, Iang=0;
  double c, s;
  for(ivarr=0;ivarr<nvarsr-1;ivarr++) {
    nang--;
    for(iang=0; iang<nang; iang++) {
      c=cos(EulAng(Iang));
      s=sin(EulAng(Iang));
      R=I;
      R(iang  ,iang  )= c;
      R(iang  ,iang+1)=-s;
      R(iang+1,iang  )= s;
      R(iang+1,iang+1)= c;
      matrix_mult(RotTemp,Rot,R,0,1.0);
      Rot=RotTemp;
      Iang++;
    }
  }
  return Rot;
}

//all coordinates are between -1 and 1
MtxDbl& gen_rand_rot_mat(MtxDbl& rot,int nvarsr) 
{
  int n_eul_ang=nchoosek(nvarsr, 2);
  //printf("n_eul_ang=%d\n",n_eul_ang);
  MtxDbl eul_ang(n_eul_ang);
  double pi=2.0*acos(0.0);
  int mymod = 1048576; //2^20 instead of 10^6 to be kind to the computer
  for(int i=0; i<n_eul_ang; ++i)
    eul_ang(i)=(std::rand() % mymod)*pi/mymod;
  rot.newSize(nvarsr,nvarsr);
  gen_rot_mat(rot, eul_ang, nvarsr);
  return rot;
}


///generates 2*nvarsr random samples between 0 and 1, the sample design is binning optimal with the bins being chosen as the end points of a randomly rotated set of axes (so the BINS are maximin spaced) the design is NOT a latin hypercube and it is not symmetric, opposite octants are sequential
MtxDbl& gen_rand_axis_bin_opt_samples_0to1(MtxDbl& xr, int nvarsr) 
{
  gen_rand_rot_mat(xr,nvarsr);
  xr.resize(2*nvarsr,nvarsr);
  int mymod = 1048576; //2^20 instead of 10^6 to be kind to the computer
  for(int i=nvarsr-1; i>=0; --i) {
    //printf("surfpack.cpp: i=%d",i);
    for(int j=0; j<nvarsr; ++j) {
      //printf(" j=%d",j); fflush(stdout);
      xr(2*i,j)=2.0*floor(1.0+xr(i,j))-1.0;
      xr(2*i+1,j)=0.5*((-xr(2*i,j)*(std::rand() % mymod))/mymod+1.0);
      xr(2*i,j)=0.5*((xr(2*i,j)*(std::rand() % mymod))/mymod+1.0);
    }
    //printf("\n");
  }
  return xr;
}

/***********************************************************************/
/// num_multi_dim_poly_coef(int Nvarsr, int Ndeg) is for use with multi_dim_poly_power(), says how many coefficients are needed for a Nvarsr-dimensional polynomial "of degree abs(Ndeg)" if Ndeg>0 the polynomial contains all multidimensional monomials UPTO (inclusive) degree Ndeg, if Ndeg<0 it contains only those multidimensional monomials OF EXACT DEGREE abs(Ndeg)
int num_multi_dim_poly_coef(int Nvarsr, int Ndeg){
  int Npoly;
  if(Ndeg<0)
    Npoly=nchoosek(-Ndeg-1+Nvarsr,-Ndeg);
  else
    Npoly=nchoosek(Ndeg+Nvarsr,Ndeg);
  return Npoly;
}

/***********************************************************************/
/**** multi_dim_poly_power() is a rather long function so I gave it ****/
/**** this special "boundary" comment so you can easily tell where  ****/
/**** it begins and ends                                            ****/
/***********************************************************************/
///poly=multi_dim_poly_power(poly,Nvarsr,Ndeg,istart,jstart) returns a matrix of size={Npoly,Nvarsr}; if Ndeg>=0 this function returns the mixed partial powers of all Nvarsr-dimensional polynomials of degree less than or equal to zero (there are Npoly=nchoosek(Ndeg+Nvarsr,Ndeg) of these), if Ndeg<=0 this function returns the mixed partial powers of all Nvarsr-dimensional polynomials of exactly degree abs(Ndeg) (there are Npoly=nchoosek(abs(Ndeg)-1+Nvarsr,abs(Ndeg)) of these); istart and jstart are offsets from the beginning of the matrix that say where to start the next "group" of powers (istart and jstart are there to avoid needing to RECURSIVELY allocate, fill, copy contents of this "poly" to the next larger "poly" and then deallocate this poly (it's a performance/speed thing), when any function other than multi_dim_poly_power() calls multi_dim_poly_power(), istart and jstart should be left unspecified, which makes them default to zero; if the "user" specifies non-zero istart and jstart then it "voids the warranty" on multi_dim_poly_power(), specifically you could end up going beyond the bounds of poly (a memory error) 
MtxInt& multi_dim_poly_power(MtxInt& poly, int Nvarsr, int Ndeg, int istart, int jstart, int iffirst){
  //printf("istart=%d jstart=%d iffirst=%d Nvarsr=%d Ndeg=%d poly.NRows()=%d\n",istart,jstart,iffirst,Nvarsr,Ndeg,poly.getNRows());
  int Npoly, npoly, istartorig=istart;
  //determine the total number of polynomials that this call of multi_dim_poly_power() is supposed to find mixed partial powers for
  if(Ndeg<0)
    Npoly=nchoosek(-Ndeg-1+Nvarsr,-Ndeg);
  else
    Npoly=nchoosek(Ndeg+Nvarsr,Ndeg);

  //istart=jstart=0 (the default when istart and jstart are not specified) is the flag for this function being called by the "user" (as opposed to a bigger multi_dim_poly_power()) and we don't want to require the user to know how big poly should be so we will right size it for him/her, if the user specified nonzero istart and/or jstart then it is his/her responsibility for making sure that poly is already big enough but since I'm a nice guy I put in an assert for when the user runs this in debug mode
  if((istart==0)&&(jstart==0)&&(iffirst==1)) {
    //printf("newSizing: Npoly=%d\n",Npoly);
    poly.newSize(Npoly,Nvarsr);
  }
  else{
#ifdef __SURFPACK_ERR_CHECK__
    if(!((istart+Npoly<=poly.getNRows())&&(jstart+Nvarsr<=poly.getNCols()))) {
      printf("Error in multi_dim_poly_power(): you asked me to fill in mixed partial polynomial powers beyond the bounds of the MtxInt& poly that you gave me.  If you want me to resize poly don't specify a nonzero istart or jstart\n, istart=%d Npoly=%d poly.NRows=%d jstart=%d Nvarsr=%d poly.NCols=%d\n",istart,Npoly,poly.getNRows(),jstart,Nvarsr,poly.getNCols());
      assert((istart+Npoly<=poly.getNRows())&&(jstart+Nvarsr<=poly.getNCols()));
    }
#endif
  }

  int i, j;
  int ndeg, nvarsr; //these are for recursive calls to smaller multi_dim_poly_power()'s that's why they start with little "n" instead of big "N"
  if(Ndeg==0) {
    // a Nvarsr-dimensional polynomial of total degree 0 has all mixed partial powers equal to zero, istart and jstart should be 0 but the user could have done something strange
    for(j=0;j<Nvarsr; j++)
      poly(istart,jstart+j)=0;
  }
  else if(Nvarsr==1) {
    //this is a failsafe for direct user input, it isn't necessary for recursive calls
    if(Ndeg>0){
      //all powers less than or equal to Ndeg for 1 variable are
      // [   0; ...
      //     1; ...
      //     2; ...
      //     3; ...
      //     .
      //     : 
      //  Ndeg];    
      // make sure to offset by istart and jstart (it shouldn't be necessary but who knows the user might have actually specified istart and jstart)
      for(i=0;i<=Ndeg;i++)
	poly(istart+i,jstart)=i;      
    }
    else{
      //Ndeg<0 says you only want powers of exactly abs(Ndeg) for 1 variable this is [-Ndeg]; istart and jstart "should be" zero but who knows what the user gave us as input
      poly(istart,jstart)=-Ndeg;
    }
  }
  else if(Ndeg==-1) {
    //Ndeg = -1 says you want all polynomials of EXACTLY degree 1, for which the "mixed" partial powers is just an identity matrix, but make sure to offset it by istart and jstart  
    for(j=0; j<Nvarsr; j++) {
      for(i=0; i<Nvarsr; i++)
	poly(istart+i,jstart+j)=0;    
      poly(istart+j,jstart+j)=1;
    }
  }
  else if(Ndeg>0) {
    //this logic should only be executed when this multi_dim_poly_power() is called by the "user" (it should not be executed when this multi_dim_poly_power() is called by a bigger multi_dim_poly_power()); Ndeg > 0 says "I want all polynomials of degree less than or equal to Ndeg" so we will do this in batches of exactly degree 0, exactly degree 1, exactly degree 2, ..., exactly degree Ndeg
    for(ndeg=0; ndeg<=Ndeg; ndeg++) {
      npoly=nchoosek(ndeg-1+Nvarsr,ndeg);
      multi_dim_poly_power(poly, Nvarsr, -ndeg, istart, jstart, 0);
      istart+=npoly;
    }
#ifdef __SURFPACK_ERR_CHECK__
    assert(istart-istartorig==Npoly); //istartorig should be zero but who knows what the user gave us
#endif      
  }
  else if(Nvarsr==2) { 
    //we know that Ndeg < -1, so we want all 2-dimensional polynomials of exaclty degree abs(Ndeg) and the answer for 2 dimensions is easy
    // [ abs(Ndeg)        0      ;...
    //   abs(Ndeg)-1      1      ;...
    //   abs(Ndeg)-2      2      ;...
    //      .             .      ;...
    //      :             :      ;...
    //      2        abs(Ndeg)-2 ;...
    //      1        abs(Ndeg)-1 ;...
    //      0        abs(Ndeg)   ]
    //done in 2 loops to fill down columns (access memory sequentially)
    for(i=0; i<=-Ndeg; i++)
      poly(istart+i,jstart  )=-Ndeg-i;
    for(i=0; i<=-Ndeg; i++)
      poly(istart+i,jstart+1)=i;
  }
  else{ //we know that Ndeg < -1 and Nvarsr > 2
    //we want all Nvarsr-dimensional polynomials of exaclty degree abs(Ndeg) so we are going to start with the first column being ideg=abs(Ndeg) decrimenting that and pairing it with ALL (Nvarsr-1)-dimensional polynomials of exactly degree ndeg=abs(Ndeg)-ideg
    int nvarsr=Nvarsr-1;
    for(int ideg=-Ndeg; ideg>=0; ideg--){
      ndeg=-Ndeg-ideg;
      npoly=nchoosek(ndeg-1+nvarsr,ndeg);
      for(i=0; i<npoly; i++)
	poly(istart+i,jstart)=ideg;
      multi_dim_poly_power(poly,nvarsr,-ndeg,istart,jstart+1,0);
      istart+=npoly;
    }
#ifdef __SURFPACK_ERR_CHECK__
    assert(istart-istartorig==Npoly); //don't expect istartorig to be zero 
#endif
  }
  
  return poly;  
}

/***********************************************************************/
/**** multi_dim_poly_power() ends here                              ****/
/***********************************************************************/

/// if ndeg>=0 generates a matrix "poly" of nvarsr-dimensional polynomial powers that for all polynomials WITHOUT MIXED POWERS up to (and including) degree ndeg.  If ndeg<0 it generates a nvarsr by nvarsr diagonal matrix "poly" of non-mixed polynomials of exact degree abs(ndeg), this diagonal matrix has abs(ndeg) for all its diagonal elements.
MtxInt& main_effects_poly_power(MtxInt& poly, int nvarsr, int ndeg) {



#ifdef __SURFPACK_ERR_CHECK__
  assert(nvarsr>0);
#endif

  if(ndeg<0) {
    int abs_ndeg=abs(ndeg);
    poly.newSize(nvarsr,nvarsr);
    poly.zero();
    for(int ivarsr=0; ivarsr<nvarsr; ++ivarsr)
      poly(ivarsr,ivarsr)=abs_ndeg;
    return poly;
  }
  else if(ndeg==0) {
    poly.newSize(1,nvarsr);
    poly.zero();
    return poly;
  }
  
  poly.newSize(1+nvarsr*ndeg,nvarsr);
  poly.zero();
  int ipoly=0;
  for(int ideg=1; ideg<=ndeg; ++ideg)
    for(int ivarsr=0; ivarsr<nvarsr; ++ivarsr)
      poly(++ipoly,ivarsr)=ideg;

  return poly;
}



//evalBasis evaluates every polynomial basis function in "poly" at every point in xr (the r is for Real), and returns them in matrix g in a reasonably efficient way; it works for arbitrary number of dimensions and polynomial order
MtxDbl& evaluate_poly_basis(MtxDbl& g, const MtxInt& poly, const MtxDbl& xr) {
  int npts=xr.getNRows();
  int nvarsr=xr.getNCols(); //number of real variables
  int npoly=poly.getNRows();
#ifdef __SURFPACK_ERR_CHECK__
  if(!((nvarsr>0)&&(nvarsr==poly.getNCols())))
    assert((nvarsr>0)&&(nvarsr==poly.getNCols()));
#endif
  g.newSize(npts,npoly);
  
  int ivarr, ipoly, ipt;
  for(ipoly=0; ipoly<npoly; ipoly++) {
#ifdef __SURFPACK_ERR_CHECK__
    for(ivarr=0;ivarr<nvarsr;ivarr++)
      assert(poly(ipoly,ivarr)>=0);
#endif
    
    //save multiplications by not doing anything for polynomial powers of zero
    for(ivarr=0;ivarr<nvarsr;ivarr++)
      if(poly(ipoly,ivarr)!=0) 
	break;

    if(ivarr==nvarsr) {
      //if all polynomial powers are zero then prod(anything^0)=1.0
      for(ipt=0; ipt<npts; ipt++)
	g(ipt,ipoly)=1.0;
    }
    else {

      //do simple assignment for the first non-zero power dimension (save assignment of 1.0 and the product of 1.0 times the power of xr)
      switch(poly(ipoly,ivarr)) {
      case 0:
	//we shouldn't need this but better safe than sorry
	break;
      case 1:
	//linear functions are so simple and common that we don't want to use the "expensive" pow function
	for(ipt=0; ipt<npts; ipt++)
	  g(ipt,ipoly)=xr(ipt,ivarr);
	break;
      case 2:
	//quadratic functions are so simple and common that we don't want to use the "expensive" pow function
	for(ipt=0; ipt<npts; ipt++)
	  g(ipt,ipoly)=xr(ipt,ivarr)*xr(ipt,ivarr);
	break;
      default:
	//we pulled out cubics as being simple and common too, but they're not as simple or common, typical use case for Kriging is power<=2, and we have to stop somewhere
	for(ipt=0; ipt<npts; ipt++)
	  g(ipt,ipoly)=pow(xr(ipt,ivarr),poly(ipoly,ivarr));
      }
      ivarr++;
     
      //do every dimension after the first non-zero power dimension in the general fashion... i.e. we have to multiply these by the non-unity product we already have
      for(; ivarr<nvarsr; ivarr++) 
	switch(poly(ipoly,ivarr)) {
	case 0:
	  //we DO need this
	  break;
	case 1:
	  //linear functions are so simple and common that we don't want to use the "expensive" pow function
	  for(ipt=0; ipt<npts; ipt++)
	    g(ipt,ipoly)*=xr(ipt,ivarr);
	  break;
	case 2:
	  //quadratic functions are so simple and common that we don't want to use the "expensive" pow function
	  for(ipt=0; ipt<npts; ipt++)
	    g(ipt,ipoly)*=xr(ipt,ivarr)*xr(ipt,ivarr);
	  break;
	default:
	  //we could have pulled out cubics as being simple and common too, but they're not as simple or common, the typical use case for Kriging is power<=2, and we have to stop somewhere
	  for(ipt=0; ipt<npts; ipt++)
	    g(ipt,ipoly)*=pow(xr(ipt,ivarr),poly(ipoly,ivarr));
	}
    }
  }
      
  return g;
}


///evaluate_poly_der_basis produces dg(ipt+ider*npts,ipoly)=d^sum(der(ider,:)) (g(ipt,ipoly))/prod(dx(ipt,ivar)^der(ider,ivar)) where g is created by evaluate_poly_basis(g,poly,xr); and npts=xr.getNRows(); is the number of points in xr. matrix_mult(dfunc,evaluate_poly_der_basis(dg,poly,der,xr),coef); dfunc.reshape(xr.getNRows(),der.getNRows()); should produce the same dfunc as evaluate_poly_der(dfunc,poly,der,coef,xr); but take more memory to do it... the primary reason for this function's existence is for use in constructing polynomial fits from gradients (and higher order derivatives if desired)
MtxDbl& evaluate_poly_der_basis(MtxDbl& dg, const MtxInt& poly, const MtxInt& der, const MtxDbl& xr){
  int nder=der.getNRows();
  int npts=xr.getNRows();
  int npoly=poly.getNRows();
  int nvars=xr.getNCols();
  assert((poly.getNCols()==nvars)&&(der.getNCols()==nvars));
  dg.newSize(npts*nder,npoly);
  dg.zero(); //initialize whole matrix to zero
  
  double tempcoef; 
  int temppolypow;
  int ivar;

  for(int ipoly=0; ipoly<npoly; ++ipoly) 
    for(int ider=0; ider<nder; ++ider) {
      int ideroffset=ider*npts;
      
      //temp coef is double so that the poly power will converts from int to 
      //double "once" instead of for every point, i.e. it's faster
      tempcoef=1.0;
      for(ivar=0; ivar<nvars; ++ivar) {
	//determine the coefficient of the derivative of this 
	//multidimensional monomial
	if(der(ider,ivar)>poly(ipoly,ivar)) {
	  //if the coefficient of the derivative polynomial is zero (because
	  //we took the derivative of a constant) we don't need to do any more
	  tempcoef=0.0;
	  break;
	}
	else{
	  //since the derivative polynomial's coefficient isn't zero (so far)
	  //we need to know what it is
	  int thisder=der(ider,ivar);
	  for(int jder=0; jder<thisder; ++jder) 
	    tempcoef*=(poly(ipoly,ivar)-jder);	  
	}
      }
      
#ifdef __SURFPACK_DER_TEST__
      printf("poly(ipoly=%3d,:)={%d",ipoly,poly(ipoly,0));
      for(int jvar=1; jvar<nvars; ++jvar)
	printf(",%d",poly(ipoly,jvar));
      printf("} der(ider=%3d,:)={%d",ider,der(ider,0));
      for(int jvar=1; jvar<nvars; ++jvar)
	printf(",%d",der(ider,jvar));
      printf("} poly(ipoly=%3d,:)-der(ider=%3d,:)={%2d",ipoly,ider,poly(ipoly,0)-der(ider,0));
      for(int jvar=1; jvar<nvars; ++jvar)
	printf(",%2d",poly(ipoly,jvar)-der(ider,jvar));
      printf("} tempcoef=%g\n",tempcoef);
#endif
      


      if(tempcoef==0.0) {
	//do nothing because the corresponding entries in dg need to be zero 
	//and we initialized the whole matrix to zero	
      }
      else{
	for(ivar=0; ivar<nvars; ++ivar) {
	  temppolypow=poly(ipoly,ivar)-der(ider,ivar);
	  if(temppolypow>0) 
	    break; //the value of ivar is retained
	}
	
	switch(temppolypow) {
	case 0:
	  //to get here ivar=nvars which means all powers are zero and 
	  //anything to the zeroth power is 1, and tempcoef*1 is tempcoef
	  //so just assign tempcoef
	  for(int ipt=0; ipt<npts; ++ipt)
	    dg(ipt+ideroffset,ipoly)=tempcoef;
	  break;
	case 1:
	  //pull out the linear term as a special case because it is 
	  //common and can be done much faster than pow(xr(ipt,ivar),1);
	  for(int ipt=0; ipt<npts; ++ipt)
	    dg(ipt+ideroffset,ipoly)=tempcoef*xr(ipt,ivar);
	  break;
	case 2:
	  //pull out the quadratic term as a special case because it is 
	  //common and can be done much faster than pow(xr(ipt,ivar),2);
	  for(int ipt=0; ipt<npts; ++ipt)
	    dg(ipt+ideroffset,ipoly)=tempcoef*xr(ipt,ivar)*xr(ipt,ivar);
	  break;
	default:
	  //do the generic "pow" case for polynomial powers>=3
	  for(int ipt=0; ipt<npts; ++ipt)
	    dg(ipt+ideroffset,ipoly)=tempcoef*pow(xr(ipt,ivar),temppolypow);
	}
	++ivar;//increase the retained value of ivar by 1
	
	//starting from the retained value of ivar multiply in the
	//contributions from the other variables
	for(; ivar<nvars; ++ivar) { 
	  temppolypow=poly(ipoly,ivar)-der(ider,ivar);
	  switch(temppolypow) {
	  case 0:
	    //anything to the zeroth power is one so save the work of
	    //multiplying by one by doing nothing
	    break;
	  case 1:
	    //pull out the linear term as a special case because it is 
	    //common and can be done much faster than pow(xr(ipt,ivar),1);
	    for(int ipt=0; ipt<npts; ++ipt)
	      dg(ipt+ideroffset,ipoly)*=xr(ipt,ivar);
	    break;
	  case 2:
	    //pull out the quadratic term as a special case because it is 
	    //common and can be done much faster than pow(xr(ipt,ivar),2);
	    for(int ipt=0; ipt<npts; ++ipt)
	      dg(ipt+ideroffset,ipoly)*=xr(ipt,ivar)*xr(ipt,ivar);
	    break;
	  default:
	    //do the generic "pow" case for polynomial powers>=3
	    for(int ipt=0; ipt<npts; ++ipt)
	      dg(ipt+ideroffset,ipoly)*=pow(xr(ipt,ivar),temppolypow);
	  } //switch(temppolypow>0)
	} //for(; ivar<nvar; ++ivar) 	
      } //if(tempcoef==0.0) {...} else
    } //for(int ider=0; ider<nder; ++ider)
  //for(int ipoly=0; ipoly<npoly; ++ipoly) doesn't have a "{" to match

  return dg;  
}


///evaluate_poly_der(dfunc,poly,der,coef,xr); should produce the same dfunc as  matrix_mult(dfunc,evaluate_poly_der_basis(dg,poly,der,xr),coef); dfunc.reshape(xr.getNRows(),der.getNRows()); but use less memory to do it
MtxDbl& evaluate_poly_der(MtxDbl& result, const MtxInt& poly, const MtxInt& der, const MtxDbl& coef, const MtxDbl& xr){
  int nder=der.getNRows();
  int npts=xr.getNRows();
  int npoly=poly.getNRows();
  int nvars=xr.getNCols();
  assert((poly.getNCols()==nvars)&&(der.getNCols()==nvars)&&
	 (coef.getNRows()==npoly)&&(coef.getNCols()==1));
  result.newSize(npts,nder);
  result.zero();
  double tempcoef;
  int temppolypow;
  MtxDbl tempres(npts);
  int ivar;


  for(int ider=0; ider<nder; ++ider) {

    for(int ipoly=0; ipoly<npoly; ++ipoly) {

      double tempcoef=coef(ipoly);
      if(tempcoef) {
	//only proceed with this monomial term if its coefficient isn't zero
	for(ivar=0; ivar<nvars; ++ivar) {
	  //determine the coefficient of the derivative of this 
	  //multidimensional monomial
	  if(der(ider,ivar)>poly(ipoly,ivar)) {
	    //if the coefficient of the derivative polynomial is zero (because
	    //we took the derivative of a constant) we don't need to do any 
	    //more
	    tempcoef=0.0;
	    break;
	  }
	  else{
	    //since the derivative polynomial's coefficient isn't zero (so far)
	    //we need to know what it is
	    int thisder=der(ider,ivar);
	    for(int jder=0; jder<thisder; ++jder) 
	      tempcoef*=(poly(ipoly,ivar)-jder);	  
	  }
	}
	
	if(tempcoef==0.0) {
	  //taking the derivative of a constant term in the multidimensional 
	  //monomial made the coefficient of derivative monomial zero so 
	  //rather than calculating the term and then multiplying by zero 
	  //(voiding the term's "contribution") we're going to take a short
	  //cut and just not calculate it. i.e. we "do nothing" for this
	  //derivative monomial term if it's coefficient is zero
	}
	else{

	  //find the first variable with a non-zero poly power in the 
	  //derivative
	  for(ivar=0; ivar<nvars; ++ivar) {
	    temppolypow=poly(ipoly,ivar)-der(ider,ivar);
	    if(temppolypow>0) 
	      break; //the value of ivar is retained
	  }

	  if(ivar==nvars) { //using the retained value of ivar
	    //all poly powers were zero and anything to the zeroth power
	    //is one and one times tempcoef is tempcoef so just add temp 
	    //coef (skip the multiplication)
	    for(int ipt=0; ipt<npts; ipt++)
	      result(ipt,ider)+=tempcoef;
	  }
	  else if(ivar==nvars-1) {
	    //the only non-zero poly power in the derivative is for the 
	    //last variable so directly add it's contribution to the result
	    //don't store it to tempres since we don't need to multiply in
	    //any contributions from other variables
	    switch(temppolypow) {
	    case 1:
	      for(int ipt=0; ipt<npts; ++ipt)
		result(ipt,ider)+=tempcoef*xr(ipt,ivar);
	      break;
	    case 2:
	      for(int ipt=0; ipt<npts; ++ipt)
		result(ipt,ider)+=tempcoef*xr(ipt,ivar)*xr(ipt,ivar);
	      break;
	    default:
	      for(int ipt=0; ipt<npts; ++ipt)
		result(ipt,ider)+=tempcoef*pow(xr(ipt,ivar),temppolypow);
	    }
	  }
	  else{
	    //there could be more variables with non-zero poly powers in
	    //the derivative polynomial, so we will store the first 
	    //variable's contribution to tempres, THEN multiply
	    //the other variables' contributions, THEN add tempres to result
	    switch(temppolypow) {
	    case 1:
	      for(int ipt=0; ipt<npts; ++ipt)
		tempres(ipt)=tempcoef*xr(ipt,ivar);
	      break;
	    case 2:
	      for(int ipt=0; ipt<npts; ++ipt)
		tempres(ipt)=tempcoef*xr(ipt,ivar)*xr(ipt,ivar);
	      break;
	    default:
	      for(int ipt=0; ipt<npts; ++ipt)
		tempres(ipt)=tempcoef*pow(xr(ipt,ivar),temppolypow);
	    }
	    ++ivar; //increase the retained value of ivar by 1

	    //starting from the retained value of ivar multiply in the
	    //contributions from the other variables
	    for(; ivar<nvars; ++ivar) { 
	      temppolypow=poly(ipoly,ivar)-der(ider,ivar);
	      switch(temppolypow) {
	      case 0:
		break;
	      case 1:
		for(int ipt=0; ipt<npts; ++ipt)
		  tempres(ipt)*=xr(ipt,ivar);
		break;
	      case 2:
		for(int ipt=0; ipt<npts; ++ipt)
		  tempres(ipt)*=xr(ipt,ivar)*xr(ipt,ivar);
		break;
	      default:
		for(int ipt=0; ipt<npts; ++ipt)
		  tempres(ipt)*=pow(xr(ipt,ivar),temppolypow);
	      } //switch(temppolypow)
	    } //for(; ivar<nvar; ++ivar)

	    //add tempres (the contribution from all variables for this
	    //multi-dimensional monomial term in the derivative polynomial)
	    //to result
	    for(int ipt=0; ipt<npts; ipt++)
	      result(ipt,ider)+=tempres(ipt);
	      

	  } //if(ivar==nvar){...}else 
	} //the 2nd if(tempcoef)
      } //the 1st if(tempcoef)
    } //for(int ipoly=0; ipoly<npoly; ++ipoly)
  } //for(int ider=0; ider<nder; ++ider)

  return result;
}











/// Throw an exception if end-of-file has been reached 
void surfpack::checkForEOF(istream& is)
{
  if (is.eof()) {
    throw surfpack::io_exception("End of file reached unexpectedly.");
  }
}

/// Return true if the file specified by parameter file name has the extension specified by parameter extension
bool surfpack::hasExtension(const string& filename, const string extension)
{
  return (filename.find(extension) == filename.size() - extension.size());
}

} // end namespace nkm
