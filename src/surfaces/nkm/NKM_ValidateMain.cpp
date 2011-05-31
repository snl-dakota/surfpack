#include "NKM_SurfMat.hpp"
#include "NKM_SurfData.hpp"
#include "NKM_KrigingModel.hpp"
#include "NKM_GradKrigingModel.hpp"
#include <iostream>
#include <cstdlib>

//#define __PROFILING_TEST__ //not iplemented yet
//#define __TIMING_BENCH__
#define __FAST_TEST__
//#define __WITH_PAV_500__
//#define __FASTER_TEST__
//#define __EVEN_FASTER_TEST__
//#define __VALGRIND_TEST__
//#define __GKM_USE_KM_CORR_LEN__
using std::cout;
using std::endl;
using std::string;

void validate();
void validate_grad();
void hack();
void check_matrix();

int main(int argc, char* argv[])
{
  //hack();
  //validate();
  validate_grad();
  //check_matrix();
  return 0;
}

void check_matrix() 
{
  nkm::MtxDbl I10(10,10), I5;

  I10.identity();
  I5.copy(I10);
  I5.resize(5,5);

  nkm::MtxDbl a(10,10),  A(10,10), AChol(10,10), Ainv(10,10), IA(10,10), IA2(10,10), IA3(10,10);
  
  int nrows=10, ncols=10;

  for(int j=0; j<nrows; ++j)
    for(int i=0; i<ncols; ++i) 
      a(i,j)=static_cast<double>(std::rand());

  nkm::matrix_mult(A,a,a,0.0,1.0,'N','T');
  AChol.copy(A);
  int info;
  double rcond;
  nkm::Chol_fact(AChol,info,rcond);
  Ainv.copy(AChol);
  nkm::inverse_after_Chol_fact(Ainv);
  nkm::matrix_mult(IA,A,Ainv);
  nkm::matrix_mult(IA2,Ainv,A);
  nkm::solve_after_Chol_fact(IA3,AChol,A);


  nkm::MtxDbl B(15,15), BChol(15,15), Binv(15,15), Binv2(15,15), IB(15,15), IB2(15,15), IB3(15,15), IB4(14,14);

  B.copy(A);
  BChol.copy(B);
  nkm::Chol_fact(BChol,info,rcond);  
  Binv.copy(BChol);
  nkm::inverse_after_Chol_fact(Binv);
  nkm::matrix_mult(IB,B,Binv);
  nkm::matrix_mult(IB2,Binv,B);
  nkm::solve_after_Chol_fact(IB3,BChol,B);
  nkm::solve_after_Chol_fact(Binv2,BChol,I10);
  nkm::matrix_mult(IB4,Binv2,B,0.0,1.0,'N','T');



  FILE *fp=fopen("test_mat.txt","w");
  nrows=I10.getNRows();
  ncols=I10.getNCols();
  fprintf(fp,"I10 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,I10.getNRowsAct(),I10.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",I10(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",I10(i,j));
  }

  nrows=I5.getNRows();
  ncols=I5.getNCols();
  fprintf(fp,"\n\nI5 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,I5.getNRowsAct(),I5.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",I5(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",I5(i,j));
  }


  nrows=a.getNRows();
  ncols=a.getNCols();
  fprintf(fp,"\n\na (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,a.getNRowsAct(),a.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",a(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",a(i,j));
  }

  nrows=A.getNRows();
  ncols=A.getNCols();
  fprintf(fp,"\n\nA (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,A.getNRowsAct(),A.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",A(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",A(i,j));
  }

  
  int j;
  nrows=AChol.getNRows();
  ncols=AChol.getNCols();
  fprintf(fp,"\n\nAChol (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,AChol.getNRowsAct(),AChol.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",AChol(i,0));
    for(j=1; j<=i; ++j)
      fprintf(fp," %12.6g",AChol(i,j));
    for(; j<ncols; ++j)
      fprintf(fp," %12.6g",0.0);
  }
  

  
  nrows=Ainv.getNRows();
  ncols=Ainv.getNCols();
  fprintf(fp,"\n\nAinv (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,Ainv.getNRowsAct(),Ainv.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",Ainv(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",Ainv(i,j));
  }

  nrows=IA.getNRows();
  ncols=IA.getNCols();
  fprintf(fp,"\n\nIA (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IA.getNRowsAct(),IA.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IA(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IA(i,j));
  }


  nrows=IA2.getNRows();
  ncols=IA2.getNCols();
  fprintf(fp,"\n\nIA2 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IA2.getNRowsAct(),IA2.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IA2(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IA2(i,j));
  }


  nrows=IA3.getNRows();
  ncols=IA3.getNCols();
  fprintf(fp,"\n\nIA3 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IA3.getNRowsAct(),IA3.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IA3(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IA3(i,j));
  }


  nrows=B.getNRows();
  ncols=B.getNCols();
  fprintf(fp,"\n\nB (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,B.getNRowsAct(),B.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",B(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",B(i,j));
  }

  nrows=BChol.getNRows();
  ncols=BChol.getNCols();
  fprintf(fp,"\n\nBChol (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,BChol.getNRowsAct(),BChol.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",BChol(i,0));
    for(j=1; j<=i; ++j)
      fprintf(fp," %12.6g",BChol(i,j));
    for(; j<ncols; ++j)
      fprintf(fp," %12.6g",0.0);
  }
  
  nrows=Binv.getNRows();
  ncols=Binv.getNCols();
  fprintf(fp,"\n\nBinv (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,Binv.getNRowsAct(),Binv.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",Binv(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",Binv(i,j));
  }

  nrows=Binv2.getNRows();
  ncols=Binv2.getNCols();
  fprintf(fp,"\n\nBinv2 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,Binv2.getNRowsAct(),Binv2.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",Binv2(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",Binv2(i,j));
  }

  nrows=IB.getNRows();
  ncols=IB.getNCols();
  fprintf(fp,"\n\nIB (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IB.getNRowsAct(),IB.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IB(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IB(i,j));
  }

  nrows=IB2.getNRows();
  ncols=IB2.getNCols();
  fprintf(fp,"\n\nIB2 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IB2.getNRowsAct(),IB2.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IB2(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IB2(i,j));
  }

  nrows=IB3.getNRows();
  ncols=IB3.getNCols();
  fprintf(fp,"\n\nIB3 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IB3.getNRowsAct(),IB3.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IB3(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IB3(i,j));
  }


  nrows=IB4.getNRows();
  ncols=IB4.getNCols();
  fprintf(fp,"\n\nIB4 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IB4.getNRowsAct(),IB4.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IB4(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IB4(i,j));
  }


  nkm::MtxDbl C(5,5),CChol(5,5), Cinv(5,5), IC(5,5);

  nrows=ncols=5;
  for(int j=0; j<ncols; ++j) 
    for(int i=0; i<nrows; ++i)
      C(i,j)=A(i,j);
  //for(int k=0; k<nrows*ncols; ++k)
  //C(k)=A(k);

  CChol.copy(C);
  nkm::Chol_fact(CChol,info,rcond);  
  Cinv.copy(CChol);
  nkm::inverse_after_Chol_fact(Cinv);
  nkm::matrix_mult(IC,C,Cinv);

  nrows=C.getNRows();
  ncols=C.getNCols();
  fprintf(fp,"\n\nC (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,C.getNRowsAct(),C.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",C(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",C(i,j));
  }

  
  nrows=CChol.getNRows();
  ncols=CChol.getNCols();
  fprintf(fp,"\n\nCChol (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,CChol.getNRowsAct(),CChol.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",CChol(i,0));
    for(j=1; j<=i; ++j)
      fprintf(fp," %12.6g",CChol(i,j));
    for(; j<ncols; ++j)
      fprintf(fp," %12.6g",0.0);
  }
    
  nrows=Cinv.getNRows();
  ncols=Cinv.getNCols();
  fprintf(fp,"\n\nCinv (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,Cinv.getNRowsAct(),Cinv.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",Cinv(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",Cinv(i,j));
  }

  nrows=IC.getNRows();
  ncols=IC.getNCols();
  fprintf(fp,"\n\nIC (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,IC.getNRowsAct(),IC.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",IC(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",IC(i,j));
  }



  nkm::MtxDbl D(15,15), DChol(15,15), Dinv(15,15), ID(15,15);

  D.copy(B); D.resize(5,5);

  DChol.copy(D);
  nkm::Chol_fact(DChol,info,rcond);  
  Dinv.copy(DChol);
  nkm::inverse_after_Chol_fact(Dinv);
  nkm::matrix_mult(ID,D,Dinv);


  nrows=D.getNRows();
  ncols=D.getNCols();
  fprintf(fp,"\n\nD (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,D.getNRowsAct(),D.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",D(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",D(i,j));
  }

  
  nrows=DChol.getNRows();
  ncols=DChol.getNCols();
  fprintf(fp,"\n\nDChol (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,DChol.getNRowsAct(),DChol.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",DChol(i,0));
    for(j=1; j<=i; ++j)
      fprintf(fp," %12.6g",DChol(i,j));
    for(; j<ncols; ++j)
      fprintf(fp," %12.6g",0.0);
  }
    
  nrows=Dinv.getNRows();
  ncols=Dinv.getNCols();
  fprintf(fp,"\n\nDinv (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,Dinv.getNRowsAct(),Dinv.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",Dinv(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",Dinv(i,j));
  }

  nrows=ID.getNRows();
  ncols=ID.getNCols();
  fprintf(fp,"\n\nID (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,ID.getNRowsAct(),ID.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",ID(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",ID(i,j));
  }

  nkm::MtxDbl E1(13,13), E2(13,13), E3(15,25);
  E1.newSize(10,10);
  nrows=ncols=10;
  for(int ij=0; ij<nrows*ncols; ++ij)
    E1(ij)=B(ij);

  nrows=E1.getNRows();
  ncols=E1.getNCols();
  fprintf(fp,"\n\nE1 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,E1.getNRowsAct(),E1.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",E1(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",E1(i,j));
  }

  E1.reshape(20,5);
  nrows=E1.getNRows();
  ncols=E1.getNCols();
  fprintf(fp,"\n\nE1 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,E1.getNRowsAct(),E1.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",E1(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",E1(i,j));
  }
  
  E2.copy(E1);
  E2.reshape(10,10);

  nrows=E2.getNRows();
  ncols=E2.getNCols();
  fprintf(fp,"\n\nE2 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,E2.getNRowsAct(),E2.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",E2(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",E2(i,j));
  }

  E3.copy(B);
  E3.reshape(5,20);

  nrows=E3.getNRows();
  ncols=E3.getNCols();
  fprintf(fp,"\n\nE3 (nrows=%d ncols=%d nrows_act=%d ncols_act=%d) =",nrows,ncols,E3.getNRowsAct(),E3.getNColsAct());
  for(int i=0; i<nrows; ++i) {
    fprintf(fp,"\n%12.6g",E3(i,0));
    for(int j=1; j<ncols; ++j)
      fprintf(fp," %12.6g",E3(i,j));
  }



  fclose(fp);

  return;
}



void hack()
{
  string filename="ORIG_DATA.spd";
  nkm::SurfData orig_data( filename , 6, 0, 10, 0, 1);
  
  std::map< std::string, std::string> km_params;
  km_params["order"] = "2";


  for(int jout=0; jout<10; ++jout) {
    orig_data.setJOut(jout);
    nkm::KrigingModel km(orig_data , km_params); km.create();
  }
  return;
}

/*
void validate_grad() {
  printf("validating Gradient Enhanced Kriging Model\n");
  
  string grad_validate2d_10 ="grad_validate2d_10.spd";
  nkm::SurfData sd2d10(grad_validate2d_10 , 2, 0, 3, 0, 1, 0);
  string grad_validate2d_100 ="grad_validate2d_100.spd";
  nkm::SurfData sd2d100(grad_validate2d_100 , 2, 0, 3, 0, 1, 0);
  string grad_validate2d_500 ="grad_validate2d_500.spd";
  nkm::SurfData sd2d500(grad_validate2d_500 , 2, 0, 3, 0, 1, 0);
  string validate2d_10K="grad_validate2d_10K.spd";
  nkm::SurfData sd2d10K(    validate2d_10K, 2, 0, 3, 0, 1, 0);
  
  nkm::MtxDbl yeval10(    10);
  nkm::MtxDbl yeval100(  100);
  nkm::MtxDbl yeval500(  500);
  nkm::MtxDbl yeval10K(10000);
  int jout=0; //the 0th output column is Rosenbrock  
  
  nkm::MtxDbl roserror(3,4); roserror.zero();
  nkm::MtxDbl roserror_grad(3,4); roserror_grad.zero();
  sd2d10.setJOut( jout);
  sd2d100.setJOut(jout);
  sd2d500.setJOut(jout);
  sd2d10K.setJOut(jout);
  
  std::map< std::string, std::string> km_params;
  km_params["lower_bounds"]="-2.0 -2.0";
  km_params["upper_bounds"]="2.0 2.0";
  //km_params["correlation_lengths"]="0.419055 2.93523";
  
  //km_params["correlation_lengths"]="1.0 2.0";
  //km_params["optimization_method"]="none";
  //km_params["optimization_method"]="local";
  km_params["order"] = "2";
  km_params["reduced_polynomial"]=nkm::toString<bool>(true);
  
  //nkm::KrigingModel kmros10( sd2d10 , km_params); kmros10.create();
  //nkm::KrigingModel kmros100( sd2d100, km_params); kmros100.create();
  //nkm::GradKrigingModel gkmros10( sd2d10 , km_params); gkmros10.create();
  nkm::GradKrigingModel gkmros100( sd2d100, km_params); gkmros100.create();
  //nkm::GradKrigingModel gkmros500( sd2d500, km_params); gkmros500.create();
  
  //km_params["optimization_method"]="local";
							  
  //km_params["correlation_lengths"]="0.53066 2.11213";
  //km_params["optimization_method"]="none";

  
  //evaluate error the 10 pt rosenbrock kriging model at 10K points
  //kmros10.evaluate(yeval10K,sd2d10K.xr);
  //for(int i=0; i<10000; ++i)
  //roserror(0,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  //roserror(0,3)=sqrt(roserror(0,2)/10000.0);

  //evaluate error the 10 pt rosenbrock kriging model at build points  
  //kmros10.evaluate(yeval10,sd2d10.xr);
  //for(int i=0; i<10; ++i)
  //roserror(0,0)+=pow(yeval10(i)-sd2d10.y(i,jout),2);
  //roserror(0,1)=sqrt(roserror(0,0)/10.0);

  //evaluate error the 10 pt rosenbrock Grad kriging model at 10K points
  //gkmros10.evaluate(yeval10K,sd2d10K.xr);
  //for(int i=0; i<10000; ++i)
  //roserror_grad(0,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  //roserror_grad(0,3)=sqrt(roserror_grad(0,2)/10000.0);

  //evaluate error the 10 pt rosenbrock Grad kriging model at build points  
  //gkmros10.evaluate(yeval10,sd2d10.xr);
  //for(int i=0; i<10; ++i)
  //roserror_grad(0,0)+=pow(yeval10(i)-sd2d10.y(i,jout),2);
  //roserror_grad(0,1)=sqrt(roserror_grad(0,0)/10.0);


  //evaluate error the 100 pt rosenbrock kriging model at 10K points
  //kmros100.evaluate(yeval10K,sd2d10K.xr);
  //for(int i=0; i<10000; ++i)
  //roserror(1,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  //roserror(1,3)=sqrt(roserror(1,2)/10000.0);

  //evaluate error the 100 pt rosenbrock kriging model at build points  
  //kmros100.evaluate(yeval100,sd2d100.xr);
  //for(int i=0; i<100; ++i)
  //roserror(1,0)+=pow(yeval100(i)-sd2d100.y(i,jout),2);
  //roserror(1,1)=sqrt(roserror(1,0)/100.0);
  
  //evaluate error the 100 pt rosenbrock Grad kriging model at 10K points
  gkmros100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
  roserror_grad(1,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  roserror_grad(1,3)=sqrt(roserror_grad(1,2)/10000.0);

  //evaluate error the 100 pt rosenbrock Grad kriging model at build points  
  gkmros100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
  roserror_grad(1,0)+=pow(yeval100(i)-sd2d100.y(i,jout),2);
  roserror_grad(1,1)=sqrt(roserror_grad(1,0)/100.0);

  
  //evaluate error the 500 pt rosenbrock Grad kriging model at 10K points
  //gkmros500.evaluate(yeval10K,sd2d10K.xr);
  //for(int i=0; i<10000; ++i)
    //roserror_grad(2,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  //roserror_grad(2,3)=sqrt(roserror_grad(2,2)/10000.0);

  //evaluate error the 500 pt rosenbrock Grad kriging model at build points  
  //gkmros500.evaluate(yeval500,sd2d500.xr);
  //for(int i=0; i<500; ++i)
    //roserror_grad(2,0)+=pow(yeval500(i)-sd2d500.y(i,jout),2);
  //roserror_grad(2,1)=sqrt(roserror_grad(2,0)/500.0);
  
  FILE *fpout=fopen("grad_Kriging.validate","w");

  
  //fprintf(fpout,"rosenbrock\n");
  //fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
  //fprintf(fpout,"Kriging:\n");
  //fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",10,roserror(0,0),roserror(0,1),roserror(0,2),roserror(0,3));  
  //fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",100,roserror(1,0),roserror(1,1),roserror(1,2),roserror(1,3));  
  
  fprintf(fpout,"Grad Kriging:\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",10,roserror_grad(0,0),roserror_grad(0,1),roserror_grad(0,2),roserror_grad(0,3));
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",100,roserror_grad(1,0),roserror_grad(1,1),roserror_grad(1,2),roserror_grad(1,3));
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",500,roserror_grad(2,0),roserror_grad(2,1),roserror_grad(2,2),roserror_grad(2,3));
  
  fclose(fpout);


  return;
}
*/

/*
void validate_grad()
{
  printf("validating Gradient Enhanced Kriging Model\n");

  //filenames

  string validate2d_100 ="grad_validate2d_100.spd";
  string validate2d_10K="validate2d_10K.spd";

  nkm::SurfData sd2d100( validate2d_100 , 2, 0, 3, 0, 1, 0);
  nkm::SurfData sd2d10K(validate2d_10K, 2, 0, 3, 0, 0, 0);

//string paviani10d_50  ="grad_paviani10d_50.spd";
//string paviani10d_10K ="paviani10d_10K.spd";



  nkm::MtxDbl yeval100(    10);
  nkm::MtxDbl yeval10K(10000);

  int jout;

  std::map< std::string, std::string> km_params;
  km_params["constraint_type"] = "r";
  //km_params["order"] = "linear";
  km_params["order"] = "2";
  km_params["reduced_polynomial"]=nkm::toString<bool>(true);

  std::map< std::string, std::string> gkm_params;
  gkm_params=km_params;
  //gkm_params["optimization_method"]="none";

  
  km_params["lower_bounds"]="-2.0 -2.0";
  km_params["upper_bounds"]="2.0 2.0";


  printf("*****************************************************************\n");
  printf("*** running shubert 2D tests ************************************\n");
  printf("*****************************************************************\n");

  nkm::MtxDbl shuerror(3,4); shuerror.zero();
  jout=1;
  sd2d100.setJOut( jout);
  sd2d10K.setJOut(jout);

  gkm_params=km_params;
  nkm::GradKrigingModel gkmshu100( sd2d100 , gkm_params); gkmshu100.create();

  //evaluate error the 10 pt shubert kriging model at 10K points
  gkmshu100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    shuerror(0,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  shuerror(0,3)=sqrt(shuerror(0,2)/10000.0);

  //evaluate error the 10 pt shubert kriging model at build points  
  gkmshu100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
    shuerror(0,0)+=pow(yeval100(i)-sd2d100.y(i,jout),2);
  shuerror(0,1)=sqrt(shuerror(0,0)/100.0);

#ifdef NONO
  printf("*****************************************************************\n");
  printf("*** running paviani 10D tests ***********************************\n");
  printf("*****************************************************************\n");

  km_params["lower_bounds"]=" 2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0";
  km_params["upper_bounds"]="10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0";

  gkm_params=km_params;

  nkm::MtxDbl paverror(3,4); paverror.zero();
  nkm::SurfData sdpav50(  paviani10d_50  , 10, 0, 1, 0, 1, 0);

  nkm::SurfData sdpav10K( paviani10d_10K , 10, 0, 1, 0, 0, 0);
  nkm::GradKrigingModel gkmpav50( sdpav50 , gkm_params); gkmpav50.create();


  nkm::MtxDbl yeval50(50);

  //evaluate error the 10 pt paviani10d kriging model at 10K points
  gkmpav50.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(0,2)+=pow(yeval10K(i)-sdpav10K.y(i),2);
  paverror(0,3)=sqrt(paverror(0,2)/10000.0);

  sdpav10K.clear();

  //evaluate error the 50 pt paviani10d kriging model at build points  
  gkmpav50.evaluate(yeval50,sdpav50.xr);
  for(int i=0; i<50; ++i)
    paverror(0,0)+=pow(yeval50(i)-sdpav50.y(i),2);
  paverror(0,1)=sqrt(paverror(0,0)/50.0);

  sdpav50.clear();
  yeval50.clear();
#endif

  printf("*****************************************************************\n");
  printf("*** writing output **********************************************\n");
  printf("*****************************************************************\n");


  FILE *fpout=fopen("grad_Kriging.validate","w");

  fprintf(fpout,"shubert\n");
  fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",10,shuerror(0,0),shuerror(0,1),shuerror(0,2),shuerror(0,3));


#ifdef NONO
  fprintf(fpout,"paviani\n");
  fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",50,paverror(0,0),paverror(0,1),paverror(0,2),paverror(0,3));
#endif

  fclose(fpout);
  return;
}
*/


std::string mtxdbl_2_string(nkm::MtxDbl& md) { 
  std::ostringstream oss;
  oss.precision(16);
  oss << md(0);
  for(int i=1; i<md.getNElems(); ++i)
    oss << " " << md(i);
  return oss.str();
}



void validate_grad()
{
  printf("validating Gradient Enhanced Kriging Model\n");

  //filenames

#ifndef __TIMING_BENCH__  
#ifndef __PROFILING_TEST__
  string validate2d_10 ="grad_validate2d_10.spd";
#endif
  string validate2d_100="grad_validate2d_100.spd";
#ifndef __PROFILING_TEST__
  string validate2d_500="grad_validate2d_500.spd";
#endif
  string validate2d_10K="validate2d_10K.spd";

#ifndef __PROFILING_TEST__
  nkm::SurfData sd2d10( validate2d_10 , 2, 0, 3, 0, 1, 0);
#endif
  nkm::SurfData sd2d100(validate2d_100, 2, 0, 3, 0, 1, 0);
#ifndef __PROFILING_TEST__
  nkm::SurfData sd2d500(validate2d_500, 2, 0, 3, 0, 1, 0);
#endif
  nkm::SurfData sd2d10K(validate2d_10K, 2, 0, 3, 0, 0, 0);
#endif

#ifndef __PROFILING_TEST__
  string paviani10d_50  ="grad_paviani10d_50.spd";
  string paviani10d_500 ="grad_paviani10d_500.spd";
  string paviani10d_2500="grad_paviani10d_2500.spd";
  string paviani10d_10K ="paviani10d_10K.spd";
#endif

#ifndef __PROFILING_TEST__
  nkm::MtxDbl yeval10(    10);
#endif
#ifndef __TIMING_BENCH__
  nkm::MtxDbl yeval100(  100);
#ifndef __PROFILING_TEST__
  nkm::MtxDbl yeval500(  500);
#endif
#endif
  nkm::MtxDbl yeval10K(10000);

  int jout;

  std::map< std::string, std::string> km_params;
  km_params["constraint_type"] = "r";
  //km_params["order"] = "linear";
  km_params["order"] = "2";
  km_params["reduced_polynomial"]=nkm::toString<bool>(true);

  std::map< std::string, std::string> gkm_params;
  gkm_params=km_params;
  //gkm_params["optimization_method"]="none";

#ifndef __TIMING_BENCH__  
  printf("*****************************************************************\n");
  printf("*** running rosenbrock 2D tests *********************************\n");
  printf("*****************************************************************\n");

  nkm::MtxDbl corr_lengths;
  nkm::MtxDbl roserror(3,4); roserror.zero();
  
  jout=0; //the 0th output column is Rosenbrock  
#ifndef __PROFILING_TEST__
  sd2d10.setJOut( jout);
#endif
  sd2d100.setJOut(jout);
#ifndef __PROFILING_TEST__
  sd2d500.setJOut(jout);
#endif
  sd2d10K.setJOut(jout);

  km_params["lower_bounds"]="-2.0 -2.0";
  km_params["upper_bounds"]="2.0 2.0";
  //km_params["optimization_method"]="local";
  //km_params["optimization_method"]="none";
  //km_params["nugget_formula"]="2";

#ifndef __PROFILING_TEST__
#ifndef __GKM_USE_KM_CORR_LEN__
  gkm_params=km_params;
#endif
  
#ifdef __GKM_USE_KM_CORR_LEN__
  nkm::KrigingModel kmros10(sd2d10,km_params); kmros10.create();
  gkm_params["correlation_lengths"]=mtxdbl_2_string(kmros10.get_correlation_lengths(corr_lengths));
#endif
  nkm::GradKrigingModel gkmros10( sd2d10 , gkm_params); gkmros10.create();
#endif

#ifndef __VALGRIND_TEST__
#ifndef __EVEN_FASTER_TEST__
#ifdef __GKM_USE_KM_CORR_LEN__
  nkm::KrigingModel kmros100(sd2d100,km_params); kmros100.create();
  gkm_params["correlation_lengths"]=mtxdbl_2_string(kmros100.get_correlation_lengths(corr_lengths));
#endif
  nkm::GradKrigingModel gkmros100(sd2d100, gkm_params); gkmros100.create();


#ifndef __FASTER_TEST__
#ifndef __PROFILING_TEST__
#ifdef __GKM_USE_KM_CORR_LEN__
  nkm::KrigingModel kmros500(sd2d500,km_params); kmros500.create();
  gkm_params["correlation_lengths"]=mtxdbl_2_string(kmros500.get_correlation_lengths(corr_lengths));  
#endif
  nkm::GradKrigingModel gkmros500(sd2d500, gkm_params); gkmros500.create();
#endif //__PROFILING_TEST__
#endif //__FASTER_TEST__
#endif //__EVEN_FASTER_TEST__
#endif
  //exit(0);

#ifndef __PROFILING_TEST__  
  //evaluate error the 10 pt rosenbrock kriging model at 10K points
  gkmros10.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    roserror(0,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  roserror(0,3)=sqrt(roserror(0,2)/10000.0);
#endif  //__PROFILING_TEST__

#ifndef __VALGRIND_TEST__
#ifndef __EVEN_FASTER_TEST__
  //evaluate error the 100 pt rosenbrock kriging model at 10K points
  gkmros100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    roserror(1,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  roserror(1,3)=sqrt(roserror(1,2)/10000.0);

#ifndef __FASTER_TEST__  
#ifndef __PROFILING_TEST__
  //evaluate error the 500 pt rosenbrock kriging model at 10K points
  gkmros500.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    roserror(2,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  roserror(2,3)=sqrt(roserror(2,2)/10000.0);
  
  //sd2d10K.clear();

  //evaluate error the 500 pt rosenbrock kriging model at build points
  gkmros500.evaluate(yeval500,sd2d500.xr);
  for(int i=0; i<500; ++i)
    roserror(2,0)+=pow(yeval500(i)-sd2d500.y(i,jout),2);
  roserror(2,1)=sqrt(roserror(2,0)/500.0);
#endif //__PROFILING_TEST__
#endif //__FASTER_TEST__


  //evaluate error the 100 pt rosenbrock kriging model at build points
  gkmros100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
    roserror(1,0)+=pow(yeval100(i)-sd2d100.y(i,jout),2);
  roserror(1,1)=sqrt(roserror(1,0)/100.0);
#endif //__EVEN_FASTER_TEST__
#endif //__VALGRIND_TEST__

#ifndef __PROFILING_TEST__
  //evaluate error the 10 pt rosenbrock kriging model at build points  
  gkmros10.evaluate(yeval10,sd2d10.xr);
  for(int i=0; i<10; ++i)
    roserror(0,0)+=pow(yeval10(i)-sd2d10.y(i,jout),2);
  roserror(0,1)=sqrt(roserror(0,0)/10.0);
  
#ifndef __VALGRIND_TEST__
  printf("*****************************************************************\n");
  printf("*** running shubert 2D tests ************************************\n");
  printf("*****************************************************************\n");

  nkm::MtxDbl shuerror(3,4); shuerror.zero();
  jout=1;
  sd2d10.setJOut( jout);
  sd2d100.setJOut(jout);
  sd2d500.setJOut(jout);
  sd2d10K.setJOut(jout);

#ifndef __GKM_USE_KM_CORR_LEN__
  gkm_params=km_params;
#endif

#ifdef __GKM_USE_KM_CORR_LEN__
  nkm::KrigingModel kmshu10(sd2d10,km_params); kmshu10.create();
  gkm_params["correlation_lengths"]=mtxdbl_2_string(kmshu10.get_correlation_lengths(corr_lengths));  
#endif
  nkm::GradKrigingModel gkmshu10( sd2d10 , gkm_params); gkmshu10.create();

#ifndef __EVEN_FASTER_TEST__
#ifdef __GKM_USE_KM_CORR_LEN__
  nkm::KrigingModel kmshu100(sd2d100,km_params); kmshu100.create();
  gkm_params["correlation_lengths"]=mtxdbl_2_string(kmshu100.get_correlation_lengths(corr_lengths));  
#endif
  nkm::GradKrigingModel gkmshu100(sd2d100, gkm_params); gkmshu100.create();

#ifndef __FASTER_TEST__
#ifdef __GKM_USE_KM_CORR_LEN__
  nkm::KrigingModel kmshu500(sd2d500,km_params); kmshu500.create();
  gkm_params["correlation_lengths"]=mtxdbl_2_string(kmshu500.get_correlation_lengths(corr_lengths));  
#endif
  nkm::GradKrigingModel gkmshu500(sd2d500, gkm_params); gkmshu500.create();
#endif //FASTER_TEST
#endif //EVEN_FASTER_TEST

  //evaluate error the 10 pt shubert kriging model at 10K points
  gkmshu10.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    shuerror(0,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  shuerror(0,3)=sqrt(shuerror(0,2)/10000.0);

#ifndef __EVEN_FASTER_TEST__
  //evaluate error the 100 pt shubert kriging model at 10K points
  gkmshu100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    shuerror(1,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  shuerror(1,3)=sqrt(shuerror(1,2)/10000.0);

#ifndef __FASTER_TEST__   
  //evaluate error the 500 pt shubert kriging model at 10K points
  gkmshu500.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    shuerror(2,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  shuerror(2,3)=sqrt(shuerror(2,2)/10000.0);
  

  //evaluate error the 500 pt shubert kriging model at build points
  gkmshu500.evaluate(yeval500,sd2d500.xr);
  for(int i=0; i<500; ++i)
    shuerror(2,0)+=pow(yeval500(i)-sd2d500.y(i,jout),2);
  shuerror(2,1)=sqrt(shuerror(2,0)/500.0);
#endif //FASTER_TEST

  //evaluate error the 100 pt shubert kriging model at build points
  gkmshu100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
    shuerror(1,0)+=pow(yeval100(i)-sd2d100.y(i,jout),2);
  shuerror(1,1)=sqrt(shuerror(1,0)/100.0);
#endif //EVEN_FASTER_TEST

  //evaluate error the 10 pt shubert kriging model at build points  
  gkmshu10.evaluate(yeval10,sd2d10.xr);
  for(int i=0; i<10; ++i)
    shuerror(0,0)+=pow(yeval10(i)-sd2d10.y(i,jout),2);
  shuerror(0,1)=sqrt(shuerror(0,0)/10.0);

  printf("*****************************************************************\n");
  printf("*** running herbie 2D tests *************************************\n");
  printf("*****************************************************************\n");

  nkm::MtxDbl herberror(3,4); herberror.zero();

  jout=2;
  sd2d10.setJOut( jout);
  sd2d100.setJOut(jout);
  sd2d500.setJOut(jout);
  sd2d10K.setJOut(jout);

#ifndef __GKM_USE_KM_CORR_LEN__
  gkm_params=km_params;
#endif

#ifdef __GKM_USE_KM_CORR_LEN__
  nkm::KrigingModel kmherb10(sd2d10,km_params); kmherb10.create();
  gkm_params["correlation_lengths"]=mtxdbl_2_string(kmherb10.get_correlation_lengths(corr_lengths));
#endif
  nkm::GradKrigingModel gkmherb10( sd2d10 , gkm_params); gkmherb10.create();

#ifndef __EVEN_FASTER_TEST__
#ifdef __GKM_USE_KM_CORR_LEN__
  nkm::KrigingModel kmherb100(sd2d100,km_params); kmherb100.create();
  gkm_params["correlation_lengths"]=mtxdbl_2_string(kmherb100.get_correlation_lengths(corr_lengths));
#endif
  nkm::GradKrigingModel gkmherb100(sd2d100, gkm_params); gkmherb100.create();

#ifndef __FASTER_TEST__
#ifdef __GKM_USE_KM_CORR_LEN__
  nkm::KrigingModel kmherb500(sd2d500,km_params); kmherb500.create();
  gkm_params["correlation_lengths"]=mtxdbl_2_string(kmherb500.get_correlation_lengths(corr_lengths));
#endif
  nkm::GradKrigingModel gkmherb500(sd2d500, gkm_params); gkmherb500.create();
#endif //FASTER_TEST
#endif //EVEN_FASTER_TEST

  //evaluate error the 10 pt herbie kriging model at 10K points
  gkmherb10.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    herberror(0,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  herberror(0,3)=sqrt(herberror(0,2)/10000.0);

#ifndef __EVEN_FASTER_TEST__
  //evaluate error the 100 pt herbie kriging model at 10K points
  gkmherb100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    herberror(1,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  herberror(1,3)=sqrt(herberror(1,2)/10000.0);

#ifndef __FASTER_TEST__   
  //evaluate error the 500 pt herbie kriging model at 10K points
  gkmherb500.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    herberror(2,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  herberror(2,3)=sqrt(herberror(2,2)/10000.0);
  

  //evaluate error the 500 pt herbie kriging model at build points
  gkmherb500.evaluate(yeval500,sd2d500.xr);
  for(int i=0; i<500; ++i)
    herberror(2,0)+=pow(yeval500(i)-sd2d500.y(i,jout),2);
  herberror(2,1)=sqrt(herberror(2,0)/500.0);
#endif //FASTER_TEST

  //evaluate error the 100 pt herbie kriging model at build points
  gkmherb100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
    herberror(1,0)+=pow(yeval100(i)-sd2d100.y(i,jout),2);
  herberror(1,1)=sqrt(herberror(1,0)/100.0);
#endif //EVEN_FASTER_TEST

  //evaluate error the 10 pt herbie kriging model at build points  
  gkmherb10.evaluate(yeval10,sd2d10.xr);
  for(int i=0; i<10; ++i)
    herberror(0,0)+=pow(yeval10(i)-sd2d10.y(i,jout),2);
  herberror(0,1)=sqrt(herberror(0,0)/10.0);
  
  sd2d10.clear();
  sd2d100.clear();
  sd2d500.clear();
  sd2d10K.clear();
  yeval10.clear();
  yeval100.clear();
#endif //__VALGRIND_TEST__
#endif //__PROFILING_TEST__
#endif //__TIMING_BENCH__  

#ifndef __PROFILING_TEST__
#ifndef __VALGRIND_TEST__
#ifndef __EVEN_FASTER_TEST__
  printf("*****************************************************************\n");
  printf("*** running paviani 10D tests ***********************************\n");
  printf("*****************************************************************\n");

  km_params["lower_bounds"]=" 2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0";
  km_params["upper_bounds"]="10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0";

#ifndef __GKM_USE_KM_CORR_LEN__
  gkm_params=km_params;
#endif

  nkm::MtxDbl paverror(3,4); paverror.zero();
#ifndef __TIMING_BENCH__
  nkm::SurfData sdpav50(  paviani10d_50  , 10, 0, 1, 0, 1, 0);
#ifndef __FASTER_TEST__
#ifdef __WITH_PAV_500__
  nkm::SurfData sdpav500( paviani10d_500 , 10, 0, 1, 0, 1, 0);
#endif //__WITH_PAV_500__
#endif //FASTER TEST
#endif //TIMING_BENCH
#ifndef __FASTER_TEST__
#ifndef __FAST_TEST__
  nkm::SurfData sdpav2500(paviani10d_2500, 10, 0, 1, 0, 1, 0);
#endif //FAST
#endif //FASTER


  nkm::SurfData sdpav10K( paviani10d_10K , 10, 0, 1, 0, 0, 0);
#ifndef __TIMING_BENCH__
#ifdef __GKM_USE_KM_CORR_LEN__
  nkm::KrigingModel kmpav50(sdpav50,km_params); kmpav50.create();
  gkm_params["correlation_lengths"]=mtxdbl_2_string(kmpav50.get_correlation_lengths(corr_lengths));
#endif
  nkm::GradKrigingModel gkmpav50( sdpav50 , gkm_params); gkmpav50.create();
#endif //TIMING_BENCH

#ifndef __FASTER_TEST__
#ifndef __TIMING_BENCH__
#ifdef __WITH_PAV_500__
#ifdef __GKM_USE_KM_CORR_LEN__
  nkm::KrigingModel kmpav500(sdpav500,km_params); kmpav500.create();
  gkm_params["correlation_lengths"]=mtxdbl_2_string(kmpav500.get_correlation_lengths(corr_lengths));
#endif
  nkm::GradKrigingModel gkmpav500(sdpav500, gkm_params); gkmpav500.create();
#endif //__WITH_PAV_500__
#endif //TIMING_BENCH

#ifndef __FAST_TEST__
#ifdef __GKM_USE_KM_CORR_LEN__
  nkm::KrigingModel kmpav2500(sdpav2500,km_params); kmpav2500.create();
  gkm_params["correlation_lengths"]=mtxdbl_2_string(kmpav2500.get_correlation_lengths(corr_lengths));
#endif
  nkm::GradKrigingModel gkmpav2500(sdpav2500, gkm_params); gkmpav2500.create();
  //cout << kmpav2500.model_summary_string();
#endif //FAST_TEST
#endif //FASTER_TEST

#ifndef __TIMING_BENCH__
  nkm::MtxDbl yeval50(50);
#endif
#ifndef __FASTER_TEST__
#ifndef __FAST_TEST__
  nkm::MtxDbl yeval2500(2500);
#endif //FAST_TEST
#endif //FASTER_TEST

#ifndef __TIMING_BENCH__
  //evaluate error the 10 pt paviani10d kriging model at 10K points
  gkmpav50.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(0,2)+=pow(yeval10K(i)-sdpav10K.y(i),2);
  paverror(0,3)=sqrt(paverror(0,2)/10000.0);

#ifndef __FASTER_TEST__
#ifdef __WITH_PAV_500__
  //evaluate error the 100 pt paviani10d kriging model at 10K points
  gkmpav500.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(1,2)+=pow(yeval10K(i)-sdpav10K.y(i),2);
  paverror(1,3)=sqrt(paverror(1,2)/10000.0);
#endif //__WITH_PAV_500__
#endif //FASTER_TEST
#endif //TIMING_BENCH

#ifndef __FASTER_TEST__
#ifndef __FAST_TEST__      
  //evaluate error the 2500 pt paviani10d kriging model at 10K points
  gkmpav2500.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(2,2)+=pow(yeval10K(i)-sdpav10K.y(i),2);
  paverror(2,3)=sqrt(paverror(2,2)/10000.0);
#endif //FAST_TEST
#endif //FASTER_TEST

#ifdef __TIMING_BENCH__  
  sdpav10K.y.copy(yeval10K);
  string pav10Kout="grad_paviani10d_10K_nkm_out.spd";
  sdpav10K.write(pav10Kout);
#endif //TIMING_BENCH

  sdpav10K.clear();

#ifndef __FASTER_TEST__
#ifndef __FAST_TEST__
  //evaluate error the 2500 pt paviani10d kriging model at build points
  gkmpav2500.evaluate(yeval2500,sdpav2500.xr);
  for(int i=0; i<2500; ++i)
    paverror(2,0)+=pow(yeval2500(i)-sdpav2500.y(i),2);
  paverror(2,1)=sqrt(paverror(2,0)/2500.0);
  
  sdpav2500.clear();
  yeval2500.clear();
#endif //__FAST_TEST__

#ifndef __TIMING_BENCH__    
#ifdef __WITH_PAV_500__
  //evaluate error the 500 pt paviani10d kriging model at build points
  gkmpav500.evaluate(yeval500,sdpav500.xr);
  for(int i=0; i<500; ++i)
    paverror(1,0)+=pow(yeval500(i)-sdpav500.y(i),2);
  paverror(1,1)=sqrt(paverror(1,0)/500.0);

  sdpav500.clear();
#endif //__WITH_PAV_500__
#endif //__TIMING_BENCH__
#endif //__FASTER_TEST__

#ifndef __TIMING_BENCH__  
  //evaluate error the 50 pt paviani10d kriging model at build points  
  gkmpav50.evaluate(yeval50,sdpav50.xr);
  for(int i=0; i<50; ++i)
    paverror(0,0)+=pow(yeval50(i)-sdpav50.y(i),2);
  paverror(0,1)=sqrt(paverror(0,0)/50.0);

  sdpav50.clear();
  yeval50.clear();

  yeval500.clear();
#endif //__TIMING_BENCH__  
#endif //__EVEN_FASTER_TEST__
#endif //__VALGRIND_TEST__
#endif //__PROFILING_TEST__
  printf("*****************************************************************\n");
  printf("*** writing output **********************************************\n");
  printf("*****************************************************************\n");

  FILE *fpout=fopen("grad_Kriging.validate","w");

#ifndef __TIMING_BENCH__
  fprintf(fpout,"rosenbrock\n");
  fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
#ifndef __PROFILING_TEST__
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",10,roserror(0,0),roserror(0,1),roserror(0,2),roserror(0,3));
#endif //__PROFILING_TEST__

#ifndef __VALGRIND_TEST__
#ifndef __EVEN_FASTER_TEST__
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",100,roserror(1,0),roserror(1,1),roserror(1,2),roserror(1,3));
#ifndef __FASTER_TEST__
#ifndef __PROFILING_TEST__
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",500,roserror(2,0),roserror(2,1),roserror(2,2),roserror(2,3));
#endif //__PROFILING_TEST__
#endif //FASTER_TEST
#endif //EVEN_FASTER_TEST

#ifndef __PROFILING_TEST__  
  fprintf(fpout,"shubert\n");
  fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",10,shuerror(0,0),shuerror(0,1),shuerror(0,2),shuerror(0,3));
#ifndef __EVEN_FASTER_TEST__
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",100,shuerror(1,0),shuerror(1,1),shuerror(1,2),shuerror(1,3));
#ifndef __FASTER_TEST__
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",500,shuerror(2,0),shuerror(2,1),shuerror(2,2),shuerror(2,3));
#endif //FASTER_TEST
#endif //EVEN_FASTER_TEST

  fprintf(fpout,"herbie\n");
  fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",10,herberror(0,0),herberror(0,1),herberror(0,2),herberror(0,3));
#ifndef __EVEN_FASTER_TEST__
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",100,herberror(1,0),herberror(1,1),herberror(1,2),herberror(1,3));
#ifndef __FASTER_TEST__
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",500,herberror(2,0),herberror(2,1),herberror(2,2),herberror(2,3));
#endif //FASTER_TEST

  fprintf(fpout,"paviani\n");
  fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",50,paverror(0,0),paverror(0,1),paverror(0,2),paverror(0,3));
#ifndef __FASTER_TEST__
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",500,paverror(1,0),paverror(1,1),paverror(1,2),paverror(1,3));
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",2500,paverror(2,0),paverror(2,1),paverror(2,2),paverror(2,3));
#endif //__FASTER_TEST__
#endif //__EVEN_FASTER_TEST__
#endif //__PROFILING_TEST__
#endif //__VALGRIND_TEST__
#endif //TIMING_BENCH

#ifndef __PROFILING_TEST__
#ifdef TIMING_BENCH
  fprintf(fpout,"paviani\n");
  fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",2500,paverror(2,0),paverror(2,1),paverror(2,2),paverror(2,3));
#endif //TIMING_BENCH
#endif //__PROFILING_TEST__

  fclose(fpout);
  return;
}


void validate()
{
  printf("validating Kriging Model\n");

  //filenames

#ifndef __TIMING_BENCH__  
  string validate2d_10 ="grad_validate2d_10.spd";
  string validate2d_100="grad_validate2d_100.spd";
  string validate2d_500="grad_validate2d_500.spd";
  string validate2d_10K="grad_validate2d_10K.spd";

  nkm::SurfData sd2d10( validate2d_10 , 2, 0, 3, 0, 1, 0);
  nkm::SurfData sd2d100(validate2d_100, 2, 0, 3, 0, 1, 0);
  nkm::SurfData sd2d500(validate2d_500, 2, 0, 3, 0, 1, 0);
  nkm::SurfData sd2d10K(validate2d_10K, 2, 0, 3, 0, 1, 0);
#endif

  string paviani10d_50  ="grad_paviani10d_50.spd";
  string paviani10d_500 ="grad_paviani10d_500.spd";
  string paviani10d_2500="grad_paviani10d_2500.spd";
  string paviani10d_10K ="grad_paviani10d_10K.spd";


  nkm::MtxDbl yeval10(    10);
#ifndef __TIMING_BENCH__
  nkm::MtxDbl yeval100(  100);
  nkm::MtxDbl yeval500(  500);
#endif
  nkm::MtxDbl yeval10K(10000);

  int jout;

  std::map< std::string, std::string> km_params;
  km_params["constraint_type"] = "r";
  //km_params["order"] = "linear";
  km_params["order"] = "2";
  km_params["reduced_polynomial"]=nkm::toString<bool>(true);


#ifndef __TIMING_BENCH__  
  printf("*****************************************************************\n");
  printf("*** running rosenbrock 2D tests *********************************\n");
  printf("*****************************************************************\n");

  nkm::MtxDbl roserror(3,4); roserror.zero();
  
  jout=0; //the 0th output column is Rosenbrock  
  sd2d10.setJOut( jout);
  sd2d100.setJOut(jout);
  sd2d500.setJOut(jout);
  sd2d10K.setJOut(jout);

  km_params["lower_bounds"]="-2.0 -2.0";
  km_params["upper_bounds"]="2.0 2.0";
  km_params["optimization_method"]="local";
  //km_params["optimization_method"]="none";
  //km_params["nugget_formula"]="2";

  nkm::KrigingModel kmros10( sd2d10 , km_params); kmros10.create();
  nkm::KrigingModel kmros100(sd2d100, km_params); kmros100.create();
  nkm::KrigingModel kmros500(sd2d500, km_params); kmros500.create();

  //exit(0);


  //evaluate error the 10 pt rosenbrock kriging model at 10K points
  kmros10.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    roserror(0,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  roserror(0,3)=sqrt(roserror(0,2)/10000.0);

  //evaluate error the 100 pt rosenbrock kriging model at 10K points
  kmros100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    roserror(1,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  roserror(1,3)=sqrt(roserror(1,2)/10000.0);
  
  //evaluate error the 500 pt rosenbrock kriging model at 10K points
  kmros500.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    roserror(2,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  roserror(2,3)=sqrt(roserror(2,2)/10000.0);
  
  //sd2d10K.clear();

  //evaluate error the 500 pt rosenbrock kriging model at build points
  kmros500.evaluate(yeval500,sd2d500.xr);
  for(int i=0; i<500; ++i)
    roserror(2,0)+=pow(yeval500(i)-sd2d500.y(i,jout),2);
  roserror(2,1)=sqrt(roserror(2,0)/500.0);

  //evaluate error the 100 pt rosenbrock kriging model at build points
  kmros100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
    roserror(1,0)+=pow(yeval100(i)-sd2d100.y(i,jout),2);
  roserror(1,1)=sqrt(roserror(1,0)/100.0);

  //evaluate error the 10 pt rosenbrock kriging model at build points  
  kmros10.evaluate(yeval10,sd2d10.xr);
  for(int i=0; i<10; ++i)
    roserror(0,0)+=pow(yeval10(i)-sd2d10.y(i,jout),2);
  roserror(0,1)=sqrt(roserror(0,0)/10.0);
  
  printf("*****************************************************************\n");
  printf("*** running shubert 2D tests ************************************\n");
  printf("*****************************************************************\n");

  nkm::MtxDbl shuerror(3,4); shuerror.zero();
  jout=1;
  sd2d10.setJOut( jout);
  sd2d100.setJOut(jout);
  sd2d500.setJOut(jout);
  sd2d10K.setJOut(jout);

  nkm::KrigingModel kmshu10( sd2d10 , km_params); kmshu10.create();
  nkm::KrigingModel kmshu100(sd2d100, km_params); kmshu100.create();
  nkm::KrigingModel kmshu500(sd2d500, km_params); kmshu500.create();


  //evaluate error the 10 pt shubert kriging model at 10K points
  kmshu10.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    shuerror(0,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  shuerror(0,3)=sqrt(shuerror(0,2)/10000.0);

  //evaluate error the 100 pt shubert kriging model at 10K points
  kmshu100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    shuerror(1,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  shuerror(1,3)=sqrt(shuerror(1,2)/10000.0);
  
  //evaluate error the 500 pt shubert kriging model at 10K points
  kmshu500.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    shuerror(2,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  shuerror(2,3)=sqrt(shuerror(2,2)/10000.0);
  

  //evaluate error the 500 pt shubert kriging model at build points
  kmshu500.evaluate(yeval500,sd2d500.xr);
  for(int i=0; i<500; ++i)
    shuerror(2,0)+=pow(yeval500(i)-sd2d500.y(i,jout),2);
  shuerror(2,1)=sqrt(shuerror(2,0)/500.0);
 

  //evaluate error the 100 pt shubert kriging model at build points
  kmshu100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
    shuerror(1,0)+=pow(yeval100(i)-sd2d100.y(i,jout),2);
  shuerror(1,1)=sqrt(shuerror(1,0)/100.0);

  //evaluate error the 10 pt shubert kriging model at build points  
  kmshu10.evaluate(yeval10,sd2d10.xr);
  for(int i=0; i<10; ++i)
    shuerror(0,0)+=pow(yeval10(i)-sd2d10.y(i,jout),2);
  shuerror(0,1)=sqrt(shuerror(0,0)/10.0);

  printf("*****************************************************************\n");
  printf("*** running herbie 2D tests *************************************\n");
  printf("*****************************************************************\n");

  nkm::MtxDbl herberror(3,4); herberror.zero();

  jout=2;
  sd2d10.setJOut( jout);
  sd2d100.setJOut(jout);
  sd2d500.setJOut(jout);
  sd2d10K.setJOut(jout);

  nkm::KrigingModel kmherb10( sd2d10 , km_params); kmherb10.create();
  nkm::KrigingModel kmherb100(sd2d100, km_params); kmherb100.create();
  nkm::KrigingModel kmherb500(sd2d500, km_params); kmherb500.create();

  //evaluate error the 10 pt herbie kriging model at 10K points
  kmherb10.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    herberror(0,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  herberror(0,3)=sqrt(herberror(0,2)/10000.0);

  //evaluate error the 100 pt herbie kriging model at 10K points
  kmherb100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    herberror(1,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  herberror(1,3)=sqrt(herberror(1,2)/10000.0);
  
  //evaluate error the 500 pt herbie kriging model at 10K points
  kmherb500.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    herberror(2,2)+=pow(yeval10K(i)-sd2d10K.y(i,jout),2);
  herberror(2,3)=sqrt(herberror(2,2)/10000.0);
  

  //evaluate error the 500 pt herbie kriging model at build points
  kmherb500.evaluate(yeval500,sd2d500.xr);
  for(int i=0; i<500; ++i)
    herberror(2,0)+=pow(yeval500(i)-sd2d500.y(i,jout),2);
  herberror(2,1)=sqrt(herberror(2,0)/500.0);

  //evaluate error the 100 pt herbie kriging model at build points
  kmherb100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
    herberror(1,0)+=pow(yeval100(i)-sd2d100.y(i,jout),2);
  herberror(1,1)=sqrt(herberror(1,0)/100.0);

  //evaluate error the 10 pt herbie kriging model at build points  
  kmherb10.evaluate(yeval10,sd2d10.xr);
  for(int i=0; i<10; ++i)
    herberror(0,0)+=pow(yeval10(i)-sd2d10.y(i,jout),2);
  herberror(0,1)=sqrt(herberror(0,0)/10.0);
  
  sd2d10.clear();
  sd2d100.clear();
  sd2d500.clear();
  sd2d10K.clear();
  yeval10.clear();
  yeval100.clear();
#endif

  printf("*****************************************************************\n");
  printf("*** running paviani 10D tests ***********************************\n");
  printf("*****************************************************************\n");

  km_params["lower_bounds"]=" 2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0";
  km_params["upper_bounds"]="10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0";


  nkm::MtxDbl paverror(3,4); paverror.zero();
#ifndef __TIMING_BENCH__
  nkm::SurfData sdpav50( paviani10d_50 , 10, 0, 1, 0, 1, 0);
  nkm::SurfData sdpav500(paviani10d_500, 10, 0, 1, 0, 1, 0);
#endif
#ifndef __FAST_TEST__
  nkm::SurfData sdpav2500(paviani10d_2500, 10, 0, 1, 0, 1, 0);
#endif
  nkm::SurfData sdpav10K(paviani10d_10K, 10, 0, 1, 0, 1, 0);
#ifndef __TIMING_BENCH__
  nkm::KrigingModel kmpav50( sdpav50 , km_params); kmpav50.create();
  nkm::KrigingModel kmpav500(sdpav500, km_params); kmpav500.create();
#endif
#ifndef __FAST_TEST__
  nkm::KrigingModel kmpav2500(sdpav2500, km_params); kmpav2500.create();
  //cout << kmpav2500.model_summary_string();
#endif

#ifndef __TIMING_BENCH__
  nkm::MtxDbl yeval50(50);
#endif
#ifndef __FAST_TEST__
  nkm::MtxDbl yeval2500(2500);
#endif

#ifndef __TIMING_BENCH__
  //evaluate error the 10 pt paviani10d kriging model at 10K points
  kmpav50.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(0,2)+=pow(yeval10K(i)-sdpav10K.y(i),2);
  paverror(0,3)=sqrt(paverror(0,2)/10000.0);

  //evaluate error the 100 pt paviani10d kriging model at 10K points
  kmpav500.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(1,2)+=pow(yeval10K(i)-sdpav10K.y(i),2);
  paverror(1,3)=sqrt(paverror(1,2)/10000.0);
#endif

#ifndef __FAST_TEST__      
  //evaluate error the 2500 pt paviani10d kriging model at 10K points
  kmpav2500.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(2,2)+=pow(yeval10K(i)-sdpav10K.y(i),2);
  paverror(2,3)=sqrt(paverror(2,2)/10000.0);
#endif

#ifdef __TIMING_BENCH__  
  sdpav10K.y.copy(yeval10K);
  string pav10Kout="paviani10d_10K_nkm_out.spd";
  sdpav10K.write(pav10Kout);
#endif

  sdpav10K.clear();

#ifndef __TIMING_BENCH__  
#ifndef __FAST_TEST__
  //evaluate error the 2500 pt paviani10d kriging model at build points
  kmpav2500.evaluate(yeval2500,sdpav2500.xr);
  for(int i=0; i<2500; ++i)
    paverror(2,0)+=pow(yeval2500(i)-sdpav2500.y(i),2);
  paverror(2,1)=sqrt(paverror(2,0)/2500.0);
  
  sdpav2500.clear();
  yeval2500.clear();
#endif
  
  //evaluate error the 500 pt paviani10d kriging model at build points
  kmpav500.evaluate(yeval500,sdpav500.xr);
  for(int i=0; i<500; ++i)
    paverror(1,0)+=pow(yeval500(i)-sdpav500.y(i),2);
  paverror(1,1)=sqrt(paverror(1,0)/500.0);

  sdpav500.clear();

  //evaluate error the 50 pt paviani10d kriging model at build points  
  kmpav50.evaluate(yeval50,sdpav50.xr);
  for(int i=0; i<50; ++i)
    paverror(0,0)+=pow(yeval50(i)-sdpav50.y(i),2);
  paverror(0,1)=sqrt(paverror(0,0)/50.0);

  sdpav50.clear();
  yeval50.clear();

  yeval500.clear();
#endif
   
  printf("*****************************************************************\n");
  printf("*** writing output **********************************************\n");
  printf("*****************************************************************\n");

  FILE *fpout=fopen("new_Kriging.validate","w");

#ifndef __TIMING_BENCH__
  fprintf(fpout,"rosenbrock\n");
  fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",10,roserror(0,0),roserror(0,1),roserror(0,2),roserror(0,3));
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",100,roserror(1,0),roserror(1,1),roserror(1,2),roserror(1,3));
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",500,roserror(2,0),roserror(2,1),roserror(2,2),roserror(2,3));
  
  fprintf(fpout,"shubert\n");
  fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",10,shuerror(0,0),shuerror(0,1),shuerror(0,2),shuerror(0,3));
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",100,shuerror(1,0),shuerror(1,1),shuerror(1,2),shuerror(1,3));
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",500,shuerror(2,0),shuerror(2,1),shuerror(2,2),shuerror(2,3));

  fprintf(fpout,"herbie\n");
  fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",10,herberror(0,0),herberror(0,1),herberror(0,2),herberror(0,3));
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",100,herberror(1,0),herberror(1,1),herberror(1,2),herberror(1,3));
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",500,herberror(2,0),herberror(2,1),herberror(2,2),herberror(2,3));
#endif

  fprintf(fpout,"paviani\n");
  fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",50,paverror(0,0),paverror(0,1),paverror(0,2),paverror(0,3));
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",500,paverror(1,0),paverror(1,1),paverror(1,2),paverror(1,3));
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",2500,paverror(2,0),paverror(2,1),paverror(2,2),paverror(2,3));
  
  fclose(fpout);
  
  return;
}
