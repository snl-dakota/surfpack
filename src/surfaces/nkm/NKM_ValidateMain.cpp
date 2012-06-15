#include "NKM_SurfMat.hpp"
#include "NKM_SurfData.hpp"
#include "NKM_KrigingModel.hpp"
#include "NKM_GradKrigingModel.hpp"
#include <iostream>
#include <iomanip>
using namespace std;
using std::ostringstream;

#include <cstdlib>

//#define __PROFILING_TEST__ //not iplemented yet
//#define __TIMING_BENCH__
#define __FAST_TEST__
//#define __WITH_PAV_500__
//#define __FASTER_TEST__
//#define __EVEN_FASTER_TEST__
//#define __VALGRIND_TEST__
//#define __GKM_USE_KM_CORR_LEN_
#define __CORR_FUNC0__
#define __CORR_FUNC1__ "matern"
//#define __CORR_FUNC2__ "0.5"
#define __CORR_FUNC2__ "infinity"
#define __GPAIS_NDIM__ 6
#define __GPAIS_GLOBAL__
//#define __GPAIS_LOCAL__
//#define __GPAIS_CORRLEN__

using std::cout;
using std::endl;
using std::string;

void validate();
void validate_grad();
void validate_grad2();
void hack();
void check_matrix();
void compare_sample_designs();
void compare_sample_designs_pav(int nvarsr);
void nested_Krig_vs_GEK_herbie_smooth_herbie_2D_4D_8D();
void boost_save_loadtest();
void global_build_and_eval_mean_adjvar();
void local_corrlen_build_and_eval_mean_adjvar();
void corrlen_build_and_eval_mean_adjvar();
void build_dont_eval();
void local_corrlen_build_dont_eval();

void gen_sample_design_by_pivoted_cholesky() {
  int NumGuesses=100;
  int Npts=2048;
  int NptsGuess=4096;
  int Ndim=2;
  std::string filename="unit_hypercube_nested_design_2D_2048pts.txt";
  int imod=104395303; //a large prime number
  double dmod=static_cast<double>(imod);
  double L=0.155*std::pow(1.0/static_cast<double>(Npts),1.0/static_cast<double>(Ndim));
  double negtheta=-0.5/(L*L);
  double dtemp;

  nkm::MtxDbl x(NptsGuess,Ndim);
  nkm::MtxDbl xfinal(Npts,Ndim);
  nkm::MtxDbl R(NptsGuess,NptsGuess);
  for(int idim=0; idim<Ndim; ++idim) {
    x(0,idim)=0.5;
    for(int i=1; i<NptsGuess; ++i)
      x(i,idim)=static_cast<double>(std::rand()%imod)/dmod;
  }
  for(int j=0; j<NptsGuess-1; ++j)
    for(int i=j+1; i<NptsGuess; ++i) {
      dtemp=x(i,0)-x(j,0);
      R(i,j)=dtemp*dtemp;
    }
  for(int idim=1; idim<Ndim-1; ++idim)
    for(int j=0; j<NptsGuess-1; ++j)
      for(int i=j+1; i<NptsGuess; ++i) {
	dtemp=x(i,idim)-x(j,idim);
	R(i,j)+=dtemp*dtemp;
      }      
  for(int j=0; j<NptsGuess; ++j) {
    R(j,j)=1.0;
    for(int i=j+1; i<NptsGuess; ++i) {
      dtemp=x(i,Ndim)-x(j,Ndim);
      R(j,i)=R(i,j)=std::exp(negtheta*(R(i,j)+dtemp*dtemp));
    }
  }
  int info=0;
  char uplo='B';
  nkm::MtxInt ipiv(NptsGuess,1);  
  int ld_R=R.getNRowsAct();
  double min_allowed_rcond=std::pow(2.0,-40.0);
  int rank=-Npts;

  PIVOTCHOL_F77(&uplo, &NptsGuess, R.ptr(0,0), &ld_R,
    		ipiv.ptr(0,0), &rank, &min_allowed_rcond, &info); 
  printf("iLoop=0 rank=%d/%d\n",rank,NptsGuess);
  if(rank>Npts)
    rank=Npts;
  for(int idim=0; idim<Ndim; ++idim)
    for(int i=0; i<rank; ++i)
      xfinal(i,idim)=x(ipiv(i,0)-1,idim);

  for(int iLoop=1; iLoop<NumGuesses; ++iLoop) {

    for(int idim=0; idim<Ndim; ++idim) {
      //don't do i=0 because it doesn't change
      for(int i=1; i<rank; ++i)
	x(i,idim)=xfinal(i,idim);
      for(int i=rank; i<NptsGuess; ++i)
	x(i,idim)=static_cast<double>(std::rand()%imod)/dmod;
    }
    
    for(int j=0; j<NptsGuess-1; ++j)
      for(int i=j+1; i<NptsGuess; ++i) {
	dtemp=x(i,0)-x(j,0);
	R(i,j)=dtemp*dtemp;
      }
    for(int idim=1; idim<Ndim-1; ++idim)
      for(int j=0; j<NptsGuess-1; ++j)
	for(int i=j+1; i<NptsGuess; ++i) {
	  dtemp=x(i,idim)-x(j,idim);
	  R(i,j)+=dtemp*dtemp;
	}      
    for(int j=0; j<NptsGuess; ++j) {
      R(j,j)=1.0;
      for(int i=j+1; i<NptsGuess; ++i) {
	dtemp=x(i,Ndim)-x(j,Ndim);
	R(j,i)=R(i,j)=std::exp(negtheta*(R(i,j)+dtemp*dtemp));
      }
    }
    
    info=0; rank=-Npts;
    PIVOTCHOL_F77(&uplo, &NptsGuess, R.ptr(0,0), &ld_R,
		  ipiv.ptr(0,0), &rank, &min_allowed_rcond, &info); 
    int itemp=rank;
    if(rank>Npts)
      rank=Npts;
    for(int idim=0; idim<Ndim; ++idim)
      for(int i=0; i<rank; ++i)
	xfinal(i,idim)=x(ipiv(i,0)-1,idim);        
    printf("iLoop=%d done rank=%d/%d\n",iLoop,itemp,NptsGuess);
  }
  
  FILE* fpout=fopen(filename.c_str(),"w");
  for(int i=0; i<Npts; ++i) {
    fprintf(fpout,"%14.12f",xfinal(i,0));
    for(int idim=1; idim<Ndim; ++idim)
      fprintf(fpout," %14.12f",xfinal(i,idim));
    fprintf(fpout,"\n");
  }
  fclose(fpout);

  return;
}

void time_build_grad() {  
  std::map< std::string, std::string> gkm_params;
  gkm_params["constraint_type"] = "r";
  gkm_params["order"] = "2";
  gkm_params["reduced_polynomial"]=nkm::toString<bool>(true);
  
  
  gkm_params["lower_bounds"]="-2.0 -2.0";
  gkm_params["upper_bounds"]="2.0 2.0";
  string paviani10d_500 ="grad_validate2d_10.spd";
  nkm::SurfData sdpav500( paviani10d_500 , 2, 0, 3, 0, 1, 0);
  /*
  gkm_params["lower_bounds"]=" 2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0";
  gkm_params["upper_bounds"]="10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0 10.0";
  string paviani10d_500 ="grad_paviani10d_500.spd";
  nkm::SurfData sdpav500( paviani10d_500 , 10, 0, 1, 0, 1, 0);
  */
  nkm::GradKrigingModel gkmpav500(sdpav500, gkm_params); 
  gkmpav500.create();
  printf("pav10D 500pt GKM: time_spent_on_pivot_cholesky=%g time_spent_on_rcond_in_pivot_cholesky=%g\n",  
	 gkmpav500.time_spent_on_pivot_cholesky,
	 gkmpav500.time_spent_on_rcond_in_pivot_cholesky);
  /*
  printf("pav10D 500pt GKM: time_spent_on_pivot_cholesky={%g,%g,(%g),%g,%g} time_spent_on_rcond_in_pivot_cholesky=%g n_pivot_cholesky_calls=%d nrcond_calls_in_pivot_cholesky=%d\n",  
	 gkmpav500.time_spent_on_pivot_cholesky_block1,
	 gkmpav500.time_spent_on_pivot_cholesky_blocks1_2,
	 gkmpav500.time_spent_on_pivot_cholesky_block4,
	 gkmpav500.time_spent_on_pivot_cholesky_blocks1_2_3,
	 gkmpav500.time_spent_on_pivot_cholesky,
	 gkmpav500.time_spent_on_rcond_in_pivot_cholesky,
	 gkmpav500.n_pivot_cholesky_calls,
	 gkmpav500.n_rcond_calls_in_pivot_cholesky);
  */
}

int main(int argc, char* argv[])
{
  //boost_save_loadtest();
  //time_build_grad();
  //compare_sample_designs_pav(2);
  //compare_sample_designs_pav(4);
  //compare_sample_designs_pav(8);
  //compare_sample_designs();
  //hack();
  //validate();
  //validate_grad();
  //validate_grad2();
  //check_matrix();
  //gen_sample_design_by_pivoted_cholesky();
  //nested_Krig_vs_GEK_herbie_smooth_herbie_2D_4D_8D();
#ifdef __GPAIS_GLOBAL__
  global_build_and_eval_mean_adjvar();
#endif
  //build_dont_eval();
#ifdef __GPAIS_LOCAL__
  local_corrlen_build_and_eval_mean_adjvar();
#endif
  //local_corrlen_build_dont_eval();
#ifdef __GPAIS_CORRLEN__
  corrlen_build_and_eval_mean_adjvar();
#endif
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
  for(int j=0; j<ncols; ++j)
    for(int i=0; i<nrows; ++i)
      E1(i,j)=B(i,j);
  //for(int ij=0; ij<nrows*ncols; ++ij)
  //E1(ij)=B(ij);

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

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
void boost_save_loadtest(){

  printf("testing boost save and load\n");

  //filenames
  string validate2d_100="grad_validate2d_100.spd";
  string validate2d_10K="grad_validate2d_10K.spd";

  nkm::SurfData sd2d100(validate2d_100, 2, 0, 3, 0, 1, 0);
  nkm::SurfData sd2d10K(validate2d_10K, 2, 0, 3, 0, 1, 0);

  nkm::MtxDbl yevalOrig(    10000,1);
  nkm::MtxDbl yevalRestored(10000,1);

  int jout;
  std::map< std::string, std::string> km_params;
  km_params["constraint_type"] = "r";
  km_params["order"] = "2";
  km_params["reduced_polynomial"]=nkm::toString<bool>(true);

  jout=0; //the 0th output column is Rosenbrock  
  sd2d100.setJOut(jout);
  sd2d10K.setJOut(jout);

  km_params["lower_bounds"]="-2.0 -2.0";
  km_params["upper_bounds"]="2.0 2.0";
  //km_params["optimization_method"]="local";
  //km_params["optimization_method"]="none";
  //km_params["nugget_formula"]="2";
  
  nkm::KrigingModel kmOrig( sd2d100 , km_params); kmOrig.create();
  nkm::SurfPackModel *kmOrigPtr=&kmOrig;
  {
    std::ofstream nkm_km_ofstream("km.sav");
  
    boost::archive::text_oarchive output_archive(nkm_km_ofstream);
    output_archive << kmOrigPtr;
  }
  nkm::SurfPackModel *kmRestored;
  {
    std::ifstream nkm_km_ifstream("km.sav");
    boost::archive::text_iarchive input_archive(nkm_km_ifstream);
    input_archive >> kmRestored;
  }
  std::cout <<"Restored: " << kmRestored->model_summary_string() << std::endl;
    
  nkm::MtxDbl xr(sd2d10K.xr);
  kmOrig.evaluate(yevalOrig,sd2d10K.xr);
  kmRestored->evaluate(yevalRestored,xr);
  delete kmRestored;

  double rmse=0;
  double mean=0;
  double var=0;
  for(int i=0; i<10000; ++i) {
    mean+=yevalOrig(i,0);
    var+=yevalOrig(i,0)*yevalOrig(i,0);
    double temp=yevalOrig(i,0)-yevalRestored(i,0);
    rmse+=temp*temp;
  }
  mean/=10000.0;
  var=var/10000.0-mean*mean;
  rmse=std::sqrt(rmse/10000.0);
  printf("rmse=%22.16g stddev=%22.16g\n",rmse,std::sqrt(var));

  sd2d100.clear();
  sd2d10K.clear();
  yevalOrig.clear();
  yevalRestored.clear();

  return;
}
#endif

void nested_Krig_vs_GEK_herbie_smooth_herbie_2D_4D_8D(){
  //string build_herbie_2D="gradHerbie_NestedLHS_2D_1024pts.spd";
  //string build_smooth_2D="gradSmoothHerbie_NestedLHS_2D_1024pts.spd";
  string build_herbie_2D="gradHerbie_NestedLHS_2D_2048pts.spd";
  string build_smooth_2D="gradSmoothHerbie_NestedLHS_2D_2048pts.spd";
  //string build_herbie_2D="gradHerbie_PivotChol_2D_1024pts.spd";
  //string build_smooth_2D="gradSmoothHerbie_PivotChol_2D_1024pts.spd";
  string valid_herbie_2D="gradHerbie_NestedLHS_2D_16384pts.spd";
  string valid_smooth_2D="gradSmoothHerbie_NestedLHS_2D_16384pts.spd";

  //string build_herbie_4D="gradHerbie_NestedLHS_4D_1024pts.spd";
  //string build_smooth_4D="gradSmoothHerbie_NestedLHS_4D_1024pts.spd";
  string build_herbie_4D="gradHerbie_NestedLHS_4D_4096pts.spd";
  string build_smooth_4D="gradSmoothHerbie_NestedLHS_4D_4096pts.spd";
  //string build_herbie_4D="gradHerbie_PivotChol_4D_2048pts.spd";
  //string build_smooth_4D="gradSmoothHerbie_PivotChol_4D_2048pts.spd";
  string valid_herbie_4D="gradHerbie_NestedLHS_4D_16384pts.spd";
  string valid_smooth_4D="gradSmoothHerbie_NestedLHS_4D_16384pts.spd";

  //string build_herbie_8D="gradHerbie_NestedLHS_8D_512pts.spd";
  //string build_smooth_8D="gradSmoothHerbie_NestedLHS_8D_512pts.spd";
  string build_herbie_8D="gradHerbie_NestedLHS_8D_2048pts.spd";
  string build_smooth_8D="gradSmoothHerbie_NestedLHS_8D_2048pts.spd";
  //string build_herbie_8D="gradHerbie_PivotChol_8D_1024pts.spd";
  //string build_smooth_8D="gradSmoothHerbie_PivotChol_8D_1024pts.spd";
  string valid_herbie_8D="gradHerbie_NestedLHS_8D_16384pts.spd";
  string valid_smooth_8D="gradSmoothHerbie_NestedLHS_8D_16384pts.spd";

  string build_herbie_filename, valid_herbie_filename;
  string build_smooth_filename, valid_smooth_filename;

  FILE *fpout1=fopen("GradKrigingPaperHerbieEffectOfDimensionStudyTableNestedLHSPivotCholKrigR.txt","w");
  FILE *fpout2=fopen("GradKrigingPaperSmoothHerbieEffectOfDimensionStudyTableNestedLHSPivotCholKrigR.txt","w");

  std::map< std::string, std::string> herbie_krig_params;    
  std::map< std::string, std::string> herbie_GEK_params;    
  std::map< std::string, std::string> smooth_krig_params;    
  std::map< std::string, std::string> smooth_GEK_params;    
  herbie_krig_params["order"] = "2";
  herbie_krig_params["reduced_polynomial"]=nkm::toString<bool>(true);
  herbie_krig_params["optimization_method"]="none";

  for(int ndimpow=1; ndimpow<=3; ++ndimpow) { //loop over the number of
    //dimensions for the effect of dimension study
    
    int Ndim=static_cast<int> (std::pow(2.0,static_cast<double>(ndimpow)));
    switch(Ndim){
    case 2:
      herbie_krig_params["lower_bounds"]="-2.0 -2.0";
      herbie_krig_params["upper_bounds"]="2.0 2.0";
      herbie_GEK_params=herbie_krig_params;
      smooth_krig_params=herbie_krig_params;
      smooth_GEK_params =herbie_krig_params;

      build_herbie_filename=build_herbie_2D;
      build_smooth_filename=build_smooth_2D;
      valid_herbie_filename=valid_herbie_2D;
      valid_smooth_filename=valid_smooth_2D;
      
      break;
    case 4:
      herbie_krig_params["lower_bounds"]="-2.0 -2.0 -2.0 -2.0";
      herbie_krig_params["upper_bounds"]="2.0 2.0 2.0 2.0";
      herbie_GEK_params=herbie_krig_params;
      smooth_krig_params=herbie_krig_params;
      smooth_GEK_params =herbie_krig_params;

      build_herbie_filename=build_herbie_4D;
      build_smooth_filename=build_smooth_4D;
      valid_herbie_filename=valid_herbie_4D;
      valid_smooth_filename=valid_smooth_4D;

      break;
    case 8:
      herbie_krig_params["lower_bounds"]="-2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0";
      herbie_krig_params["upper_bounds"]="2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0";
      herbie_GEK_params=herbie_krig_params;
      smooth_krig_params=herbie_krig_params;
      smooth_GEK_params =herbie_krig_params;

      build_herbie_filename=build_herbie_8D;
      build_smooth_filename=build_smooth_8D;
      valid_herbie_filename=valid_herbie_8D;
      valid_smooth_filename=valid_smooth_8D;

      break;
    default:
      std::cerr << "Error: haven't coded for the " << Ndim << " dimensional case." << std::endl;
      assert(false);
    } 
    nkm::SurfData sd_build_herbie(build_herbie_filename, Ndim, 0, 1, 0, 1, 0);
    nkm::SurfData sd_build_smooth(build_smooth_filename, Ndim, 0, 1, 0, 1, 0);
    nkm::SurfData sd_valid_herbie(valid_herbie_filename, Ndim, 0, 1, 0, 1, 0);
    nkm::SurfData sd_valid_smooth(valid_smooth_filename, Ndim, 0, 1, 0, 1, 0);
    nkm::SurfData sd_build_temp;
    
    
    int NptsBuild=sd_build_herbie.getNPts();
    int NptsValid=sd_valid_herbie.getNPts();
    assert((NptsBuild==sd_build_smooth.getNPts())&&
	   (NptsValid==sd_valid_smooth.getNPts()));
    int Nref=static_cast<int>(std::log(static_cast<double>(NptsBuild)/static_cast<double>(2*Ndim))/std::log(2.0));
    //if(Nref>2) Nref=2; //fast test for debug
    nkm::MtxDbl yeval(NptsValid,1);
    nkm::MtxInt ipts(NptsBuild,1);
    nkm::MtxDbl error_metric(Nref+1,9);  error_metric.zero();
    for(int i=0; i<NptsBuild; ++i)
      ipts(i,0)=i;
    for(int iref=0; iref<=Nref; ++iref) { //nested sample design loop
      int NptsThis=2*Ndim*static_cast<int>(std::pow(2.0,static_cast<double>(iref)));
      if(true) {
	//make this run faster by feeding it the correlation lengths generated from the nested LHS design the first time I rank it (so can quickly calculate the mean absolute error, originally I only calculated the RMSE)
	switch(Ndim){
	case 2:
	  switch(NptsThis){
	  case 4:
	    herbie_GEK_params["correlation_lengths"] ="0.500044 0.500132";
	    smooth_GEK_params["correlation_lengths"] ="0.500044 0.500132";
	    break;
	  case 8:
	    herbie_krig_params["correlation_lengths"]="1.26045  0.353585";
	    herbie_GEK_params["correlation_lengths"] ="0.490945 0.399543";
	    smooth_krig_params["correlation_lengths"]="2.01378  0.353585";
	    smooth_GEK_params["correlation_lengths"] ="0.896517 0.918406";
	    break;
	  case 16:
	    herbie_krig_params["correlation_lengths"]="2.11395  0.481163";
	    herbie_GEK_params["correlation_lengths"] ="0.269972 0.327727";
	    smooth_krig_params["correlation_lengths"]="0.250022 0.535061";
	    smooth_GEK_params["correlation_lengths"] ="0.901693 0.916258";
	    break;
	  case 32:
	    herbie_krig_params["correlation_lengths"]="0.337192 0.405882";
	    herbie_GEK_params["correlation_lengths"] ="0.364354 0.349339";
	    smooth_krig_params["correlation_lengths"]="0.801452 0.818999";
	    smooth_GEK_params["correlation_lengths"] ="0.954075 0.964208";
	    break;
	  case 64:
	    herbie_krig_params["correlation_lengths"]="0.362288 0.386812";
	    herbie_GEK_params["correlation_lengths"] ="0.394655 0.396955";
	    smooth_krig_params["correlation_lengths"]="0.942504 0.946496";
	    smooth_GEK_params["correlation_lengths"] ="0.708602 0.843246";
	    break;
	  case 128:
	    herbie_krig_params["correlation_lengths"]="0.341277 0.345084";
	    herbie_GEK_params["correlation_lengths"] ="0.393593 0.378305";
	    smooth_krig_params["correlation_lengths"]="0.713824 0.793784";
	    smooth_GEK_params["correlation_lengths"] ="0.681042 0.641354";
	    break;
	  case 256:
	    herbie_krig_params["correlation_lengths"]="0.427906 0.428056";
	    herbie_GEK_params["correlation_lengths"] ="0.343611 0.345492";
	    smooth_krig_params["correlation_lengths"]="0.409262 0.631571";
	    smooth_GEK_params["correlation_lengths"] ="0.662672 0.702565";
	    break;
	  case 512:
	    herbie_krig_params["correlation_lengths"]="0.332171 0.335345";
	    herbie_GEK_params["correlation_lengths"] ="0.328333 0.323684";
	    smooth_krig_params["correlation_lengths"]="0.223043 0.532758";
	    smooth_GEK_params["correlation_lengths"] ="0.686387 0.673577";
	    break;
	  case 1024:
	    herbie_krig_params["correlation_lengths"]="0.182501 0.255595";
	    herbie_GEK_params["correlation_lengths"] ="0.295556 0.261237";
	    smooth_krig_params["correlation_lengths"]="0.107725 0.473198";
	    smooth_GEK_params["correlation_lengths"] ="0.636192 0.639449";
	    break;
	  case 2048:
	    herbie_krig_params["correlation_lengths"]="0.176487 0.118631";
	    smooth_krig_params["correlation_lengths"]="0.176487 0.118631";
	    break;
	  default:
	    std::cerr << "NptsThis is not a known size" << std::endl;
	    assert(false);
	  } //switch(NptsThis)
	  break;
	case 4:
	  switch(NptsThis){
	  case 8:
	    herbie_GEK_params["correlation_lengths"] ="0.986542 1.25723 0.752366  1.20457";
	    smooth_GEK_params["correlation_lengths"] ="1.20457  1.10578 0.763173  1.05946";
	    break;
	  case 16:
	    herbie_krig_params["correlation_lengths"]="0.632662 0.503578 15.6613  1.30938";
	    herbie_GEK_params["correlation_lengths"] ="0.817832 0.903696 0.606163 0.817832";
	    smooth_krig_params["correlation_lengths"]="0.510812 0.510812 1.57612  0.719313";
	    smooth_GEK_params["correlation_lengths"] ="0.903696 0.817832 0.750758 0.817832";
	    break;
	  case 32:
	    herbie_krig_params["correlation_lengths"]="0.524469 13.1696  13.1696  0.579533";
	    herbie_GEK_params["correlation_lengths"] ="0.687712 0.770831 0.579533 0.658908";
	    smooth_krig_params["correlation_lengths"]="0.50972  4.32949  13.3587  0.55526";
	    smooth_GEK_params["correlation_lengths"] ="0.888995 0.816084 0.816084 0.839698";
	    break;
	  case 64:
	    herbie_krig_params["correlation_lengths"]="0.530866 1.75907  0.473623 0.508631";
	    herbie_GEK_params["correlation_lengths"] ="0.473623 0.578295 0.603575 0.629961";
	    smooth_krig_params["correlation_lengths"]="0.603575 0.376989 11.2333  0.508631";
	    smooth_GEK_params["correlation_lengths"] ="0.887095 0.912763 0.887095 0.925875";
	    break;
	  case 128:
	    herbie_krig_params["correlation_lengths"]="0.745955 0.937167 0.317008 0.602285";
	    herbie_GEK_params["correlation_lengths"] ="0.446403 0.43385  0.446403 0.507544";
	    smooth_krig_params["correlation_lengths"]="0.602285 1.61135  0.360427 0.593756";
	    smooth_GEK_params["correlation_lengths"] ="0.964285 0.950629 0.964285 0.964285";
	    break;
	  case 256:
	    herbie_krig_params["correlation_lengths"]="0.776901 0.683312 0.683312 0.799381";
	    herbie_GEK_params["correlation_lengths"] ="0.408916 0.408916 0.445449 0.458339";
	    smooth_krig_params["correlation_lengths"]="0.788061 0.776901 0.776901 0.810863";
	    smooth_GEK_params["correlation_lengths"] ="0.962224 0.962224 0.962224 0.962224";
	    break;
	  case 512:
	    herbie_krig_params["correlation_lengths"]="0.62593  0.574595 0.653293 0.672196";
	    herbie_GEK_params["correlation_lengths"] ="0.408042 0.390952 0.408042 0.396567";
	    smooth_krig_params["correlation_lengths"]="0.919951 0.919951 0.919951 0.906924";
	    smooth_GEK_params["correlation_lengths"] ="1.00214  0.960168 0.98795  1.00214";
	    break;
	  case 1024:
	    herbie_krig_params["correlation_lengths"]="0.549352 0.462937 0.504297 0.565248";
	    herbie_GEK_params["correlation_lengths"] ="0.373776 0.373776 0.379145 0.373776";
	    smooth_krig_params["correlation_lengths"]="0.958116 0.944548 0.958116 0.958116";
	    smooth_GEK_params["correlation_lengths"] ="0.842697 0.842697 0.842697 0.89217";
	    break;
	  case 2048:
	    herbie_krig_params["correlation_lengths"]="0.482142 0.357356 0.412136 0.482142";
	    herbie_GEK_params["correlation_lengths"] ="0.4426   0.436332 0.4426   0.4426";
	    smooth_krig_params["correlation_lengths"]="0.956068 0.956068 0.956068 0.956068";
	    smooth_GEK_params["correlation_lengths"] ="0.840896 0.840896 0.840896 0.840896";
	    break;
	  case 4096:
	    herbie_krig_params["correlation_lengths"]="0.39969  0.38845  0.37218  0.405432";
	    smooth_krig_params["correlation_lengths"]="0.875781 0.875781 0.875781 0.851153";
	    break;
	  default:
	    std::cerr << "NptsThis is not a known size" << std::endl;
	    assert(false);
	  } //switch(NptsThis)
	  break;
	case 8:
	  switch(NptsThis){
	  case 16:
	    herbie_GEK_params["correlation_lengths"] ="1.25992  1.25992  1.62868  1.25992  1.25992  1.25992  1.85175  1.25992";
	    smooth_GEK_params["correlation_lengths"] ="0.753977 1.25992  1.25992  0.857244 1.25992  1.25992  1.25992  1.25992";
	    break;
	  case 32:
	    herbie_krig_params["correlation_lengths"]="0.786096 17.1154  0.6914   17.1154  17.1154  17.1154  1.15535  17.1154";
	    herbie_GEK_params["correlation_lengths"] ="0.786096 0.786096 0.786096 0.786096 0.6914   0.786096 0.786096 0.786096";
	    smooth_krig_params["correlation_lengths"]="1.69806  17.1154  0.6914   1.69806  17.1154  17.1154  0.786096 1.15535";
	    smooth_GEK_params["correlation_lengths"] ="0.786096 1.15535  1.15535  0.786096 1.15535  1.15535  0.893762 0.786096";
	    break;
	  case 64:
	    herbie_krig_params["correlation_lengths"]="0.634017 4.94358  4.94358  2.28857  15.6949  1.05946  0.720853 15.6949";
	    herbie_GEK_params["correlation_lengths"] ="0.720853 0.720853 0.720853 0.720853 0.720853 0.720853 0.720853 0.819584";
	    smooth_krig_params["correlation_lengths"]="0.634017 4.94358  4.94358  4.94358  15.6949  1.05946  0.720853 15.6949";
	    smooth_GEK_params["correlation_lengths"] ="0.720853 1.05946  1.05946  1.05946  1.05946  1.05946  0.931836 1.05946";
	    break;
	  case 128:
	    herbie_krig_params["correlation_lengths"]="0.661025 0.581396 1.42789  1.42789  0.661025 1.42789  0.661025 14.3923";
	    herbie_GEK_params["correlation_lengths"] ="0.661025 0.661025 0.661025 0.661025 0.661025 0.661025 0.661025 0.751561";
	    smooth_krig_params["correlation_lengths"]="0.661025 0.661025 0.661025 1.42789  0.661025 0.661025 0.661025 0.581396";
	    smooth_GEK_params["correlation_lengths"] ="0.971532 0.971532 0.971532 0.971532 0.971532 0.971532 0.854497 0.971532";
	    break;
	  case 256:
	    herbie_krig_params["correlation_lengths"]="0.606163 1.30938  1.30938  0.606163 0.606163 1.30938  0.606163 0.533142";
	    herbie_GEK_params["correlation_lengths"] ="0.890899 0.890899 0.890899 0.890899 0.890899 0.890899 0.783578 0.890899";
	    smooth_krig_params["correlation_lengths"]="1.30938  0.533142 1.30938  0.606163 0.606163 0.606163 13.1978  0.606163";
	    smooth_GEK_params["correlation_lengths"] ="0.890899 0.890899 0.890899 1.01292  0.890899 0.890899 0.890899 0.890899";
	    break;
	  case 512:
	    herbie_krig_params["correlation_lengths"]="0.488894 1.20071  1.20071  0.555854 1.20071  1.20071  0.555854 0.555854";
	    herbie_GEK_params["correlation_lengths"] ="0.816958 0.816958 0.816958 0.816958 0.816958 0.816958 0.718544 0.816958";
	    smooth_krig_params["correlation_lengths"]="0.555854 0.555854 1.20071  1.20071  1.20071  1.20071  0.555854 0.488894";
	    smooth_GEK_params["correlation_lengths"] ="0.816958 0.816958 0.816958 0.928851 0.816958 0.816958 0.816958 0.816958";
	    break;
	  case 1024:
	    herbie_krig_params["correlation_lengths"]="1.25186  0.50972  3.49564  0.50972  1.10106  0.50972  1.10106  0.50972";
	    herbie_GEK_params["correlation_lengths"] ="0.749154 0.749154 0.749154 0.749154 0.749154 0.749154 0.749154 0.749154";
	    smooth_krig_params["correlation_lengths"]="2.37841  0.50972  5.13766  0.50972  1.25186  0.749154 0.749154 0.50972";
	    smooth_GEK_params["correlation_lengths"] ="0.749154 0.749154 0.749154 0.85176  0.749154 0.749154 0.749154 0.749154";
	    break;
	  case 2048:
	    herbie_krig_params["correlation_lengths"]="0.686977 0.686977 0.531434 0.467416 4.71125  3.20551  0.686977 0.686977";
	    smooth_krig_params["correlation_lengths"]="2.18102  0.467416 0.531434 0.467416 3.20551  2.18102  0.686977 0.686977";
	    break;
	  default:
	    std::cerr << "NptsThis is not a known size" << std::endl;
	    assert(false);
	  } //switch(NptsThis) 
	  break;
	default:
	  std::cerr << "Ndim is not a known number of dimensions" << std::endl;
	  assert(false);
	} //switch{Ndim) 
      }


      error_metric(iref,0)=static_cast<double>(NptsThis);
      ipts.resize(NptsThis,1); //relies on actual and apparent sizes of the matrix 
      //class being different and that resize() doesn't copy or overwrite or shrink 
      //or enlarge unless it needs a bigger size than it actually has OR the user 
      //"forces" it to resize, neither of these 2 cases are true in the original 
      //implementation.of this function

      {//limit km and gkm for herbie to this scope

	sd_build_herbie.getPoints(sd_build_temp,ipts);
	printf("Herbie: Ndim=%d Npts=%4d/%-4d",Ndim,NptsThis,NptsBuild);
	fprintf(fpout1,"Herbie: Ndim=%d Npts=%4d/%-4d",Ndim,NptsThis,NptsBuild);

	if(iref>0) { //if have enough equations for a reduced quadratic trend
	  nkm::KrigingModel km(sd_build_temp,herbie_krig_params); km.create();
	  km.evaluate(yeval,sd_valid_herbie.xr);
	  for(int i=0; i<NptsValid; ++i) {
	    double tmpdbl=yeval(i,0)-sd_valid_herbie.y(i,0);
	    error_metric(iref,1)+=std::fabs(tmpdbl);
	    error_metric(iref,2)+=std::pow(tmpdbl,2.0);
	  }
	  error_metric(iref,1)/=static_cast<double>(NptsValid);
	  error_metric(iref,2)=std::sqrt(error_metric(iref,2)/static_cast<double>(NptsValid));
	  printf(" Krig_MAE=%12.6g Krig_RMSE=%12.6g",error_metric(iref,1),error_metric(iref,2));
	  fprintf(fpout1," Krig_MAE=%12.6g Krig_RMSE=%12.6g",error_metric(iref,1),error_metric(iref,2));
	}
	else{
	  printf(" Krig_MAE=NaN          Krig_RMSE=NaN         ");
	  fprintf(fpout1," Krig_MAE=NaN          Krig_RMSE=NaN         ");
	}
	fflush(fpout1);
	if(iref<Nref) { //Npts*(1+Ndim) equations makes for a BIG correlation matrix 
	  //(slow emulator construction) and I don't need the largest Npts for the
	  //Gradient Enhanced Kriging Paper
	  nkm::GradKrigingModel gkm(sd_build_temp,herbie_GEK_params); gkm.create();
	  gkm.evaluate(yeval,sd_valid_herbie.xr);
	  for(int i=0; i<NptsValid; ++i) {
	    double tmpdbl=yeval(i,0)-sd_valid_herbie.y(i,0);
	    error_metric(iref,3)+=std::fabs(tmpdbl);
	    error_metric(iref,4)+=std::pow(tmpdbl,2.0);
	  }
	  error_metric(iref,3)/=static_cast<double>(NptsValid);	
	  error_metric(iref,4)=std::sqrt(error_metric(iref,4)/static_cast<double>(NptsValid));	
	  printf(" GEK_MAE=%12.6g GEK_RMSE=%12.6g %5d/%-5d\n",
		 error_metric(iref,3),error_metric(iref,4),
		 gkm.getNumEqnKeep(),gkm.getNumEqnAvail());
	  fprintf(fpout1," GEK_MAE=%12.6g GEK_RMSE=%12.6g %5d/%-5d\n",
		  error_metric(iref,3),error_metric(iref,4),
		  gkm.getNumEqnKeep(),gkm.getNumEqnAvail());
	}
	else{
	  printf(" GEK_MAE=NaN          GEK_RMSE=NaN            NaN/NaN  \n");
	  fprintf(fpout1," GEK_MAE=NaN          GEK_RMSE=NaN            NaN/NaN  \n");
	}	
	fflush(fpout1);
      } //end herbie scope

      {//limit km and gkm for SMOOTH herbie to this scope
	sd_build_smooth.getPoints(sd_build_temp,ipts);
	printf("Smooth: Ndim=%d Npts=%4d/%-4d",Ndim,NptsThis,NptsBuild);
	fprintf(fpout2,"Smooth: Ndim=%d Npts=%4d/%-4d",Ndim,NptsThis,NptsBuild);
	
	if(iref>0) { //if have enough equations for a reduced quadratic trend
	  nkm::KrigingModel km(sd_build_temp,smooth_krig_params); km.create();
	  km.evaluate(yeval,sd_valid_smooth.xr);
	  for(int i=0; i<NptsValid; ++i) {
	    double tmpdbl=yeval(i,0)-sd_valid_smooth.y(i,0);
	    error_metric(iref,5)+=std::fabs(tmpdbl);
	    error_metric(iref,6)+=std::pow(tmpdbl,2.0);
	  }
	  error_metric(iref,5)/=static_cast<double>(NptsValid);
	  error_metric(iref,6)=std::sqrt(error_metric(iref,6)/static_cast<double>(NptsValid));
	  printf(" Krig_MAE=%12.6g Krig_RMSE=%12.6g",error_metric(iref,5),error_metric(iref,6));
	  fprintf(fpout2," Krig_MAE=%12.6g Krig_RMSE=%12.6g",error_metric(iref,5),error_metric(iref,6));
	}
	else{
	  printf(" Krig_MAE=NaN          Krig_RMSE=NaN         ");
	  fprintf(fpout2," Krig_MAE=NaN          Krig_RMSE=NaN         ");
	}
	fflush(fpout2);	
	if(iref<Nref) { //Npts*(1+Ndim) equations makes for a BIG correlation matrix 
	  //(slow emulator construction) and I don't need the largest Npts for the
	  //Gradient Enhanced Kriging Paper

	  nkm::GradKrigingModel gkm(sd_build_temp,smooth_GEK_params); gkm.create();
	  gkm.evaluate(yeval,sd_valid_smooth.xr);
	  for(int i=0; i<NptsValid; ++i) {
	    double tmpdbl=yeval(i,0)-sd_valid_smooth.y(i,0);
	    error_metric(iref,7)+=std::fabs(tmpdbl);
	    error_metric(iref,8)+=std::pow(tmpdbl,2.0);
	  }
	  error_metric(iref,7)/=static_cast<double>(NptsValid);	
	  error_metric(iref,8)=std::sqrt(error_metric(iref,8)/static_cast<double>(NptsValid));	
	  printf(" GEK_MAE=%12.6g GEK_RMSE=%12.6g %5d/%-5d\n",
		 error_metric(iref,7),error_metric(iref,8),
		 gkm.getNumEqnKeep(),gkm.getNumEqnAvail());
	  fprintf(fpout2," GEK_MAE=%12.6g GEK_RMSE=%12.6g %5d/%-5d\n",
		  error_metric(iref,7),error_metric(iref,8),
		  gkm.getNumEqnKeep(),gkm.getNumEqnAvail());
	}
	else{
	  printf(" GEK_MAE=NaN          GEK_RMSE=NaN            NaN/NaN  \n");
	  fprintf(fpout2," GEK_MAE=NaN          GEK_RMSE=NaN            NaN/NaN  \n");
	}	
	fflush(fpout2);
      } //end SMOOTH herbie scope
    } //for(int iref=0; iref<=Nref; ++iref)
    
  } //for(int ndimpow=1; ndimpow<3; ++ndimpow)

  fclose(fpout1);
  fclose(fpout2);
  return;
} //end of function

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

void compare_sample_designs() {
  
  string build_filename ="build_file.spd";
  nkm::SurfData sd2dbuild(build_filename , 2, 0, 3, 0, 1, 0);
  string valid_filename ="valid_file.spd";
  nkm::SurfData sd2dvalid(valid_filename , 2, 0, 3, 0, 1, 0);
  FILE* fpout=fopen("compare_out.txt","w");
  
  nkm::MtxDbl yeval(16384,1);
  int jout; //the 0th output column is Rosenbrock    
  double rmse;


  std::map< std::string, std::string> km_params;
  km_params["lower_bounds"]="-2.0 -2.0";
  km_params["upper_bounds"]="2.0 2.0";

  km_params["order"] = "2";
  km_params["reduced_polynomial"]=nkm::toString<bool>(true);
  

  jout=0;
  sd2dbuild.setJOut(jout);
  sd2dvalid.setJOut(jout);
  nkm::KrigingModel kmros( sd2dbuild, km_params); kmros.create();

  //evaluate error the rosenbrock kriging model at 2^14=16384 validation points
  kmros.evaluate(yeval,sd2dvalid.xr);
  rmse=0.0;
  for(int i=0; i<16384; ++i)
    rmse+=std::pow(yeval(i,0)-sd2dvalid.y(i,jout),2);
  rmse=std::sqrt(rmse/16384.0);
  fprintf(fpout,"%22.16g\n",rmse);


  jout=1;
  sd2dbuild.setJOut(jout);
  sd2dvalid.setJOut(jout);
  nkm::KrigingModel kmshu( sd2dbuild, km_params); kmshu.create();

  //evaluate error the shubert kriging model at 2^14=16384 validation points
  kmshu.evaluate(yeval,sd2dvalid.xr);
  rmse=0.0;
  for(int i=0; i<16384; ++i)
    rmse+=std::pow(yeval(i,0)-sd2dvalid.y(i,jout),2);
  rmse=std::sqrt(rmse/16384.0);
  fprintf(fpout,"%22.16g\n",rmse);


  jout=2;
  sd2dbuild.setJOut(jout);
  sd2dvalid.setJOut(jout);
  nkm::KrigingModel kmherb( sd2dbuild, km_params); kmherb.create();

  //evaluate error the herbie kriging model at 2^14=16384 validation points
  kmherb.evaluate(yeval,sd2dvalid.xr);
  rmse=0.0;
  for(int i=0; i<16384; ++i)
    rmse+=std::pow(yeval(i,0)-sd2dvalid.y(i,jout),2);
  rmse=std::sqrt(rmse/16384.0);
  fprintf(fpout,"%22.16g\n",rmse);

  fclose(fpout);


  return;
}

void compare_sample_designs_pav(int nvarsr) {
  
  string build_filename ="build_file.spd";
  nkm::SurfData sd_pav_build(build_filename , nvarsr, 0, 1, 0, 1, 0);
  string valid_filename ="valid_file.spd";
  nkm::SurfData sd_pav_valid(valid_filename , nvarsr, 0, 1, 0, 1, 0);
  FILE* fpout=fopen("compare_out.txt","w");

  nkm::MtxDbl yeval(16384,1);
  double rmse;

  std::map< std::string, std::string> km_params;
  {  
    ostringstream os;  
    os << "2.0";
    for(int i=1; i<nvarsr; ++i)
      os << " 2.0";
    km_params["lower_bounds"]=os.str();
  }
  {  
    ostringstream os;  
    os << "10.0";
    for(int i=1; i<nvarsr; ++i)
      os << " 10.0";
    km_params["upper_bounds"]=os.str();
  }

  km_params["order"] = "2";
  km_params["reduced_polynomial"]=nkm::toString<bool>(true);
  
  nkm::KrigingModel km( sd_pav_build, km_params); km.create();

  //evaluate error the rosenbrock kriging model at 2^14=16384 validation points
  km.evaluate(yeval,sd_pav_valid.xr);
  rmse=0.0;
  for(int i=0; i<16384; ++i)
    rmse+=std::pow(yeval(i,0)-sd_pav_valid.y(i,0),2);
  rmse=std::sqrt(rmse/16384.0);
  fprintf(fpout,"%22.16g\n",rmse);


  fclose(fpout);


  return;
}


void validate_grad2() {
  printf("validating Gradient Enhanced Kriging Model\n");
  //string buildfilename ="grad_validate2d_8a.spd";
  //string buildfilename ="grad_validate2d_16a.spd";
  string buildfilename ="grad_validate2d_32a.spd";
  //string buildfilename ="dakota_sbo_rosen_10_first11.spd";
  //nkm::SurfData sdbuild(buildfilename, 2, 0, 1, 0, 1, 0);
  nkm::SurfData sdbuild(buildfilename, 2, 0, 3, 0, 1, 0);
  string validfilename="grad_validate2d_10K.spd";
  nkm::SurfData sdvalid(validfilename, 2, 0, 3, 0, 1, 0);
  
  int NptsBuild=sdbuild.getNPts();
  nkm::MtxDbl yevalbuild(NptsBuild,1);
  nkm::MtxDbl yeval10K(10000,1);
  int jout=0; //the 0th output column is Rosenbrock  
  //int jout=2; //the 2nd output column is herbie  
  
  nkm::MtxDbl roserror_grad(1,4); roserror_grad.zero();
  nkm::MtxDbl roserror(1,4); roserror.zero();
  sdbuild.setJOut(jout);
  sdvalid.setJOut(jout);
  
  std::map< std::string, std::string> km_params;
  //km_params["lower_bounds"]="-1.4 0.8";
  //km_params["upper_bounds"]="-1.0 1.2";
  km_params["lower_bounds"]="-2 -2";
  km_params["upper_bounds"]= "2 2";
  km_params["order"] = "2";
  km_params["reduced_polynomial"]=nkm::toString<bool>(true);
  
  nkm::GradKrigingModel gkmros(sdbuild, km_params); gkmros.create();
  nkm::KrigingModel kmros(sdbuild, km_params); kmros.create();

  
  //evaluate error the NptsBuild pt rosenbrock grad kriging model at 10K points
  gkmros.evaluate(yeval10K,sdvalid.xr);
  for(int i=0; i<10000; ++i)
    roserror_grad(0,2)+=std::pow(yeval10K(i,0)-sdvalid.y(i,jout),2);
  roserror_grad(0,3)=std::sqrt(roserror_grad(0,2)/10000.0);

  //evaluate error the NptsBuild pt rosenbrock Grad kriging model at build points  
  gkmros.evaluate(yevalbuild,sdbuild.xr);
  for(int i=0; i<NptsBuild; ++i)
    roserror_grad(0,0)+=std::pow(yevalbuild(i,0)-sdbuild.y(i,0),2);
  roserror_grad(0,1)=std::sqrt(roserror_grad(0,0)/NptsBuild);

  //evaluate error the NptsBuild pt rosenbrock grad kriging model at 10K points
  kmros.evaluate(yeval10K,sdvalid.xr);
  for(int i=0; i<10000; ++i)
    roserror(0,2)+=std::pow(yeval10K(i,0)-sdvalid.y(i,jout),2);
  roserror(0,3)=std::sqrt(roserror(0,2)/10000.0);

  //evaluate error the NptsBuild pt rosenbrock Grad kriging model at build points  
  kmros.evaluate(yevalbuild,sdbuild.xr);
  for(int i=0; i<NptsBuild; ++i)
    roserror(0,0)+=std::pow(yevalbuild(i,0)-sdbuild.y(i,0),2);
  roserror(0,1)=std::sqrt(roserror(0,0)/NptsBuild);


  FILE *fpout=fopen("grad_Kriging.validate","w");
  
  fprintf(fpout,"rosenbrock\n");
  fprintf(fpout,"Grad Kriging:\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",NptsBuild,roserror_grad(0,0),roserror_grad(0,1),roserror_grad(0,2),roserror_grad(0,3));
  fprintf(fpout,"Kriging:\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",NptsBuild,roserror(0,0),roserror(0,1),roserror(0,2),roserror(0,3));

  
  fclose(fpout);

  return;
}


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



  nkm::MtxDbl yeval100(    10,1);
  nkm::MtxDbl yeval10K(10000,1);

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
    shuerror(0,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  shuerror(0,3)=std::sqrt(shuerror(0,2)/10000.0);

  //evaluate error the 10 pt shubert kriging model at build points  
  gkmshu100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
    shuerror(0,0)+=std::pow(yeval100(i,0)-sd2d100.y(i,jout),2);
  shuerror(0,1)=std::sqrt(shuerror(0,0)/100.0);

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


  nkm::MtxDbl yeval50(50,1);

  //evaluate error the 10 pt paviani10d kriging model at 10K points
  gkmpav50.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(0,2)+=std::pow(yeval10K(i,0)-sdpav10K.y(i,0),2);
  paverror(0,3)=std::sqrt(paverror(0,2)/10000.0);

  sdpav10K.clear();

  //evaluate error the 50 pt paviani10d kriging model at build points  
  gkmpav50.evaluate(yeval50,sdpav50.xr);
  for(int i=0; i<50; ++i)
    paverror(0,0)+=std::pow(yeval50(i,0)-sdpav50.y(i,0),2);
  paverror(0,1)=std::sqrt(paverror(0,0)/50.0);

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
  oss << md(0,0);
  for(int i=1; i<md.getNRows(); ++i)
    oss << " " << md(i,0);
  for(int j=1; j<md.getNCols(); ++j)
    for(int i=0; i<md.getNElems(); ++i)
      oss << " " << md(i,j);
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
  nkm::MtxDbl yeval10(    10,1);
#endif
#ifndef __TIMING_BENCH__
  nkm::MtxDbl yeval100(  100,1);
#ifndef __PROFILING_TEST__
  nkm::MtxDbl yeval500(  500,1);
#endif
#endif
  nkm::MtxDbl yeval10K(10000,1);

  int jout;

  std::map< std::string, std::string> km_params;
  km_params["constraint_type"] = "r";
  //km_params["order"] = "linear";
  //km_params["order"] = "2";
  //km_params["reduced_polynomial"]=nkm::toString<bool>(true);


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
    roserror(0,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  roserror(0,3)=std::sqrt(roserror(0,2)/10000.0);
#endif  //__PROFILING_TEST__

#ifndef __VALGRIND_TEST__
#ifndef __EVEN_FASTER_TEST__
  //evaluate error the 100 pt rosenbrock kriging model at 10K points
  gkmros100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    roserror(1,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  roserror(1,3)=std::sqrt(roserror(1,2)/10000.0);

#ifndef __FASTER_TEST__  
#ifndef __PROFILING_TEST__
  //evaluate error the 500 pt rosenbrock kriging model at 10K points
  gkmros500.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    roserror(2,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  roserror(2,3)=std::sqrt(roserror(2,2)/10000.0);
  
  //sd2d10K.clear();

  //evaluate error the 500 pt rosenbrock kriging model at build points
  gkmros500.evaluate(yeval500,sd2d500.xr);
  for(int i=0; i<500; ++i)
    roserror(2,0)+=std::pow(yeval500(i,0)-sd2d500.y(i,jout),2);
  roserror(2,1)=std::sqrt(roserror(2,0)/500.0);
#endif //__PROFILING_TEST__
#endif //__FASTER_TEST__


  //evaluate error the 100 pt rosenbrock kriging model at build points
  gkmros100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
    roserror(1,0)+=std::pow(yeval100(i,0)-sd2d100.y(i,jout),2);
  roserror(1,1)=std::sqrt(roserror(1,0)/100.0);
#endif //__EVEN_FASTER_TEST__
#endif //__VALGRIND_TEST__

#ifndef __PROFILING_TEST__
  //evaluate error the 10 pt rosenbrock kriging model at build points  
  gkmros10.evaluate(yeval10,sd2d10.xr);
  for(int i=0; i<10; ++i)
    roserror(0,0)+=std::pow(yeval10(i,0)-sd2d10.y(i,jout),2);
  roserror(0,1)=std::sqrt(roserror(0,0)/10.0);
  
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
    shuerror(0,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  shuerror(0,3)=std::sqrt(shuerror(0,2)/10000.0);

#ifndef __EVEN_FASTER_TEST__
  //evaluate error the 100 pt shubert kriging model at 10K points
  gkmshu100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    shuerror(1,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  shuerror(1,3)=std::sqrt(shuerror(1,2)/10000.0);

#ifndef __FASTER_TEST__   
  //evaluate error the 500 pt shubert kriging model at 10K points
  gkmshu500.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    shuerror(2,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  shuerror(2,3)=std::sqrt(shuerror(2,2)/10000.0);
  

  //evaluate error the 500 pt shubert kriging model at build points
  gkmshu500.evaluate(yeval500,sd2d500.xr);
  for(int i=0; i<500; ++i)
    shuerror(2,0)+=std::pow(yeval500(i,0)-sd2d500.y(i,jout),2);
  shuerror(2,1)=std::sqrt(shuerror(2,0)/500.0);
#endif //FASTER_TEST

  //evaluate error the 100 pt shubert kriging model at build points
  gkmshu100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
    shuerror(1,0)+=std::pow(yeval100(i,0)-sd2d100.y(i,jout),2);
  shuerror(1,1)=std::sqrt(shuerror(1,0)/100.0);
#endif //EVEN_FASTER_TEST

  //evaluate error the 10 pt shubert kriging model at build points  
  gkmshu10.evaluate(yeval10,sd2d10.xr);
  for(int i=0; i<10; ++i)
    shuerror(0,0)+=std::pow(yeval10(i,0)-sd2d10.y(i,jout),2);
  shuerror(0,1)=std::sqrt(shuerror(0,0)/10.0);

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
    herberror(0,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  herberror(0,3)=std::sqrt(herberror(0,2)/10000.0);

#ifndef __EVEN_FASTER_TEST__
  //evaluate error the 100 pt herbie kriging model at 10K points
  gkmherb100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    herberror(1,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  herberror(1,3)=std::sqrt(herberror(1,2)/10000.0);

#ifndef __FASTER_TEST__   
  //evaluate error the 500 pt herbie kriging model at 10K points
  gkmherb500.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    herberror(2,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  herberror(2,3)=std::sqrt(herberror(2,2)/10000.0);
  

  //evaluate error the 500 pt herbie kriging model at build points
  gkmherb500.evaluate(yeval500,sd2d500.xr);
  for(int i=0; i<500; ++i)
    herberror(2,0)+=std::pow(yeval500(i,0)-sd2d500.y(i,jout),2);
  herberror(2,1)=std::sqrt(herberror(2,0)/500.0);
#endif //FASTER_TEST

  //evaluate error the 100 pt herbie kriging model at build points
  gkmherb100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
    herberror(1,0)+=std::pow(yeval100(i,0)-sd2d100.y(i,jout),2);
  herberror(1,1)=std::sqrt(herberror(1,0)/100.0);
#endif //EVEN_FASTER_TEST

  //evaluate error the 10 pt herbie kriging model at build points  
  gkmherb10.evaluate(yeval10,sd2d10.xr);
  for(int i=0; i<10; ++i)
    herberror(0,0)+=std::pow(yeval10(i,0)-sd2d10.y(i,jout),2);
  herberror(0,1)=std::sqrt(herberror(0,0)/10.0);
  
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
  //printf("pav10D 50pt GKM: time_spent_on_pivot_cholesky=%g time_spent_on_rcond_in_pivot_cholesky=%g n_pivot_cholesky_calls=%d nrcond_calls_in_pivot_cholesky=%d\n",gkmpav50.time_spent_on_pivot_cholesky,gkmpav50.time_spent_on_rcond_in_pivot_cholesky,gkmpav50.n_pivot_cholesky_calls,gkmpav50.n_rcond_calls_in_pivot_cholesky);

#endif //TIMING_BENCH


#ifndef __FASTER_TEST__
#ifndef __TIMING_BENCH__
#ifdef __WITH_PAV_500__
#ifdef __GKM_USE_KM_CORR_LEN__
  nkm::KrigingModel kmpav500(sdpav500,km_params); kmpav500.create();
  gkm_params["correlation_lengths"]=mtxdbl_2_string(kmpav500.get_correlation_lengths(corr_lengths));
#endif
  nkm::GradKrigingModel gkmpav500(sdpav500, gkm_params); gkmpav500.create();
  printf("pav10D 500pt GKM: time_spent_on_pivot_cholesky={%g,%g,(%g),%g,%g} time_spent_on_rcond_in_pivot_cholesky=%g n_pivot_cholesky_calls=%d nrcond_calls_in_pivot_cholesky=%d\n",  
	 gkmpav500.time_spent_on_pivot_cholesky_block1,
	 gkmpav500.time_spent_on_pivot_cholesky_blocks1_2,
	 gkmpav500.time_spent_on_pivot_cholesky_block4,
	 gkmpav500.time_spent_on_pivot_cholesky_blocks1_2_3,
	 gkmpav500.time_spent_on_pivot_cholesky,
	 gkmpav500.time_spent_on_rcond_in_pivot_cholesky,
	 gkmpav500.n_pivot_cholesky_calls,
	 gkmpav500.n_rcond_calls_in_pivot_cholesky);

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
  nkm::MtxDbl yeval50(50,1);
#endif
#ifndef __FASTER_TEST__
#ifndef __FAST_TEST__
  nkm::MtxDbl yeval2500(2500,1);
#endif //FAST_TEST
#endif //FASTER_TEST

#ifndef __TIMING_BENCH__
  //evaluate error the 10 pt paviani10d kriging model at 10K points
  gkmpav50.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(0,2)+=std::pow(yeval10K(i,0)-sdpav10K.y(i,0),2);
  paverror(0,3)=std::sqrt(paverror(0,2)/10000.0);

#ifndef __FASTER_TEST__
#ifdef __WITH_PAV_500__
  //evaluate error the 100 pt paviani10d kriging model at 10K points
  gkmpav500.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(1,2)+=std::pow(yeval10K(i,0)-sdpav10K.y(i,0),2);
  paverror(1,3)=std::sqrt(paverror(1,2)/10000.0);
#endif //__WITH_PAV_500__
#endif //FASTER_TEST
#endif //TIMING_BENCH

#ifndef __FASTER_TEST__
#ifndef __FAST_TEST__      
  //evaluate error the 2500 pt paviani10d kriging model at 10K points
  gkmpav2500.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(2,2)+=std::pow(yeval10K(i,0)-sdpav10K.y(i,0),2);
  paverror(2,3)=std::sqrt(paverror(2,2)/10000.0);
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
    paverror(2,0)+=std::pow(yeval2500(i,0)-sdpav2500.y(i,0),2);
  paverror(2,1)=std::sqrt(paverror(2,0)/2500.0);
  
  sdpav2500.clear();
  yeval2500.clear();
#endif //__FAST_TEST__

#ifndef __TIMING_BENCH__    
#ifdef __WITH_PAV_500__
  //evaluate error the 500 pt paviani10d kriging model at build points
  gkmpav500.evaluate(yeval500,sdpav500.xr);
  for(int i=0; i<500; ++i)
    paverror(1,0)+=std::pow(yeval500(i,0)-sdpav500.y(i,0),2);
  paverror(1,1)=std::sqrt(paverror(1,0)/500.0);

  sdpav500.clear();
#endif //__WITH_PAV_500__
#endif //__TIMING_BENCH__
#endif //__FASTER_TEST__

#ifndef __TIMING_BENCH__  
  //evaluate error the 50 pt paviani10d kriging model at build points  
  gkmpav50.evaluate(yeval50,sdpav50.xr);
  for(int i=0; i<50; ++i)
    paverror(0,0)+=std::pow(yeval50(i,0)-sdpav50.y(i,0),2);
  paverror(0,1)=std::sqrt(paverror(0,0)/50.0);

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
#ifndef __VALGRIND_TEST__
  string validate2d_10 ="grad_validate2d_10.spd";
  string validate2d_500="grad_validate2d_500.spd";
  nkm::SurfData sd2d10( validate2d_10 , 2, 0, 3, 0, 1, 0);
  nkm::SurfData sd2d500(validate2d_500, 2, 0, 3, 0, 1, 0);
#endif
  string validate2d_100="grad_validate2d_100.spd";
  string validate2d_10K="grad_validate2d_10K.spd";
  nkm::SurfData sd2d100(validate2d_100, 2, 0, 3, 0, 1, 0);
  nkm::SurfData sd2d10K(validate2d_10K, 2, 0, 3, 0, 1, 0);
#endif

#ifndef __VALGRIND_TEST__
  string paviani10d_50  ="grad_paviani10d_50.spd";
  string paviani10d_500 ="grad_paviani10d_500.spd";
  string paviani10d_2500="grad_paviani10d_2500.spd";
  string paviani10d_10K ="grad_paviani10d_10K.spd";
#endif

#ifndef __VALGRIND_TEST__
  nkm::MtxDbl yeval10(    10,1);
#endif
#ifndef __TIMING_BENCH__
  nkm::MtxDbl yeval100(  100,1);
#ifndef __VALGRIND_TEST__
  nkm::MtxDbl yeval500(  500,1);
#endif
#endif
  nkm::MtxDbl yeval10K(10000,1);

  int jout;

  std::map< std::string, std::string> km_params;
  km_params["constraint_type"] = "r";
  //km_params["order"] = "linear";
  km_params["order"] = "2";
  km_params["reduced_polynomial"]=nkm::toString<bool>(true);
  //km_params["powered_exponential"]="1";
  //km_params["powered_exponential"]="1.5";
  //km_params["powered_exponential"]="2";
  //km_params["matern"] = "0.5";
  //km_params["matern"] = "1.5";
  //km_params["matern"] = "2.5";
  //km_params["matern"] = "infinity";

#ifndef __TIMING_BENCH__  
  printf("*****************************************************************\n");
  printf("*** running rosenbrock 2D tests *********************************\n");
  printf("*****************************************************************\n");

  nkm::MtxDbl roserror(3,4); roserror.zero();
  
  jout=0; //the 0th output column is Rosenbrock  
#ifndef __VALGRIND_TEST__
  sd2d10.setJOut( jout);
#endif
  sd2d100.setJOut(jout);
#ifndef __VALGRIND_TEST__
  sd2d500.setJOut(jout);
#endif

  sd2d10K.setJOut(jout);

  km_params["lower_bounds"]="-2.0 -2.0";
  km_params["upper_bounds"]="2.0 2.0";
  //km_params["optimization_method"]="local";
  //km_params["optimization_method"]="none";
  //km_params["nugget_formula"]="2";

#ifndef __VALGRIND_TEST__
  nkm::KrigingModel kmros10( sd2d10 , km_params); kmros10.create();
#endif
  nkm::KrigingModel kmros100(sd2d100, km_params); kmros100.create();
#ifndef __VALGRIND_TEST__
  nkm::KrigingModel kmros500(sd2d500, km_params); kmros500.create();
#endif


  //exit(0);


#ifndef __VALGRIND_TEST__
  //evaluate error the 10 pt rosenbrock kriging model at 10K points
  kmros10.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    roserror(0,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  roserror(0,3)=std::sqrt(roserror(0,2)/10000.0);
#endif

  //evaluate error the 100 pt rosenbrock kriging model at 10K points
  kmros100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    roserror(1,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  roserror(1,3)=std::sqrt(roserror(1,2)/10000.0);

#ifndef __VALGRIND_TEST__  
  //evaluate error the 500 pt rosenbrock kriging model at 10K points
  kmros500.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    roserror(2,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  roserror(2,3)=std::sqrt(roserror(2,2)/10000.0);
  
  //sd2d10K.clear();

  //evaluate error the 500 pt rosenbrock kriging model at build points
  kmros500.evaluate(yeval500,sd2d500.xr);
  for(int i=0; i<500; ++i)
    roserror(2,0)+=std::pow(yeval500(i,0)-sd2d500.y(i,jout),2);
  roserror(2,1)=std::sqrt(roserror(2,0)/500.0);
#endif

  //evaluate error the 100 pt rosenbrock kriging model at build points
  kmros100.evaluate(yeval100,sd2d100.xr);
  int N100=sd2d100.getNPts();
  for(int i=0; i<N100; ++i)
    roserror(1,0)+=std::pow(yeval100(i,0)-sd2d100.y(i,jout),2);
  roserror(1,1)=std::sqrt(roserror(1,0)/N100);

#ifndef __VALGRIND_TEST__
  //evaluate error the 10 pt rosenbrock kriging model at build points  
  kmros10.evaluate(yeval10,sd2d10.xr);
  for(int i=0; i<10; ++i)
    roserror(0,0)+=std::pow(yeval10(i,0)-sd2d10.y(i,jout),2);
  roserror(0,1)=std::sqrt(roserror(0,0)/10.0);
  
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
    shuerror(0,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  shuerror(0,3)=std::sqrt(shuerror(0,2)/10000.0);

  //evaluate error the 100 pt shubert kriging model at 10K points
  kmshu100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    shuerror(1,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  shuerror(1,3)=std::sqrt(shuerror(1,2)/10000.0);
  
  //evaluate error the 500 pt shubert kriging model at 10K points
  kmshu500.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    shuerror(2,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  shuerror(2,3)=std::sqrt(shuerror(2,2)/10000.0);
  

  //evaluate error the 500 pt shubert kriging model at build points
  kmshu500.evaluate(yeval500,sd2d500.xr);
  for(int i=0; i<500; ++i)
    shuerror(2,0)+=std::pow(yeval500(i,0)-sd2d500.y(i,jout),2);
  shuerror(2,1)=std::sqrt(shuerror(2,0)/500.0);
 

  //evaluate error the 100 pt shubert kriging model at build points
  kmshu100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
    shuerror(1,0)+=std::pow(yeval100(i,0)-sd2d100.y(i,jout),2);
  shuerror(1,1)=std::sqrt(shuerror(1,0)/100.0);

  //evaluate error the 10 pt shubert kriging model at build points  
  kmshu10.evaluate(yeval10,sd2d10.xr);
  for(int i=0; i<10; ++i)
    shuerror(0,0)+=std::pow(yeval10(i,0)-sd2d10.y(i,jout),2);
  shuerror(0,1)=std::sqrt(shuerror(0,0)/10.0);

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
    herberror(0,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  herberror(0,3)=std::sqrt(herberror(0,2)/10000.0);

  //evaluate error the 100 pt herbie kriging model at 10K points
  kmherb100.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    herberror(1,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  herberror(1,3)=std::sqrt(herberror(1,2)/10000.0);
  
  //evaluate error the 500 pt herbie kriging model at 10K points
  kmherb500.evaluate(yeval10K,sd2d10K.xr);
  for(int i=0; i<10000; ++i)
    herberror(2,2)+=std::pow(yeval10K(i,0)-sd2d10K.y(i,jout),2);
  herberror(2,3)=std::sqrt(herberror(2,2)/10000.0);
  

  //evaluate error the 500 pt herbie kriging model at build points
  kmherb500.evaluate(yeval500,sd2d500.xr);
  for(int i=0; i<500; ++i)
    herberror(2,0)+=std::pow(yeval500(i,0)-sd2d500.y(i,jout),2);
  herberror(2,1)=std::sqrt(herberror(2,0)/500.0);

  //evaluate error the 100 pt herbie kriging model at build points
  kmherb100.evaluate(yeval100,sd2d100.xr);
  for(int i=0; i<100; ++i)
    herberror(1,0)+=std::pow(yeval100(i,0)-sd2d100.y(i,jout),2);
  herberror(1,1)=std::sqrt(herberror(1,0)/100.0);

  //evaluate error the 10 pt herbie kriging model at build points  
  kmherb10.evaluate(yeval10,sd2d10.xr);
  for(int i=0; i<10; ++i)
    herberror(0,0)+=std::pow(yeval10(i,0)-sd2d10.y(i,jout),2);
  herberror(0,1)=std::sqrt(herberror(0,0)/10.0);
  
  sd2d10.clear();
  sd2d500.clear();
  yeval10.clear();
#endif
  sd2d100.clear();
  sd2d10K.clear();
  yeval100.clear();
#endif
#ifndef __VALGRIND_TEST__
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
  nkm::MtxDbl yeval50(50,1);
#endif
#ifndef __FAST_TEST__
  nkm::MtxDbl yeval2500(2500,1);
#endif

#ifndef __TIMING_BENCH__
  //evaluate error the 10 pt paviani10d kriging model at 10K points
  kmpav50.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(0,2)+=std::pow(yeval10K(i,0)-sdpav10K.y(i,0),2);
  paverror(0,3)=std::sqrt(paverror(0,2)/10000.0);

  //evaluate error the 100 pt paviani10d kriging model at 10K points
  kmpav500.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(1,2)+=std::pow(yeval10K(i,0)-sdpav10K.y(i,0),2);
  paverror(1,3)=std::sqrt(paverror(1,2)/10000.0);
#endif

#ifndef __FAST_TEST__      
  //evaluate error the 2500 pt paviani10d kriging model at 10K points
  kmpav2500.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(2,2)+=std::pow(yeval10K(i,0)-sdpav10K.y(i,0),2);
  paverror(2,3)=std::sqrt(paverror(2,2)/10000.0);
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
    paverror(2,0)+=std::pow(yeval2500(i,0)-sdpav2500.y(i,0),2);
  paverror(2,1)=std::sqrt(paverror(2,0)/2500.0);
  
  sdpav2500.clear();
  yeval2500.clear();
#endif
  
  //evaluate error the 500 pt paviani10d kriging model at build points
  kmpav500.evaluate(yeval500,sdpav500.xr);
  for(int i=0; i<500; ++i)
    paverror(1,0)+=std::pow(yeval500(i,0)-sdpav500.y(i,0),2);
  paverror(1,1)=std::sqrt(paverror(1,0)/500.0);

  sdpav500.clear();

  //evaluate error the 50 pt paviani10d kriging model at build points  
  kmpav50.evaluate(yeval50,sdpav50.xr);
  for(int i=0; i<50; ++i)
    paverror(0,0)+=std::pow(yeval50(i,0)-sdpav50.y(i,0),2);
  paverror(0,1)=std::sqrt(paverror(0,0)/50.0);

  sdpav50.clear();
  yeval50.clear();

  yeval500.clear();
#endif
#endif   
  printf("*****************************************************************\n");
  printf("*** writing output **********************************************\n");
  printf("*****************************************************************\n");

  FILE *fpout=fopen("new_Kriging.validate","w");

#ifndef __TIMING_BENCH__
  fprintf(fpout,"rosenbrock\n");
  fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
#ifndef __VALGRIND_TEST__
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",10,roserror(0,0),roserror(0,1),roserror(0,2),roserror(0,3));
#endif
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",100,roserror(1,0),roserror(1,1),roserror(1,2),roserror(1,3));
#ifndef __VALGRIND_TEST__
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
#endif
#ifndef __VALGRIND_TEST__
  fprintf(fpout,"paviani\n");
  fprintf(fpout,"# of samples, SSE at build points, RMSE at build points, SSE at 10K points, RMSE at 10K points\n");
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",50,paverror(0,0),paverror(0,1),paverror(0,2),paverror(0,3));
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",500,paverror(1,0),paverror(1,1),paverror(1,2),paverror(1,3));
  fprintf(fpout,"%12d, %19.6g, %20.6g, %17.6g, %18.6g\n",2500,paverror(2,0),paverror(2,1),paverror(2,2),paverror(2,3));
#endif  
  fclose(fpout);
  
  return;
}

//used for quick develop/test of matlab implementation of adaptive importance sampling
void global_build_and_eval_mean_adjvar() {
  int Nvarsr=__GPAIS_NDIM__;
  std::map< std::string, std::string> km_params;
  km_params["constraint_type"] = "r";
  km_params["order"] = "0";
  km_params["optimization_method"]="global";
#ifdef __CORR_FUNC0__
  km_params[__CORR_FUNC1__]=__CORR_FUNC2__;
#endif  

  //km_params["reduced_polynomial"]=nkm::toString<bool>(true);
  string buildfile="buildfile.spd";
  string evalfile_in ="evalfile_in.spd";
  //string evalfile_out  ="evalfile.out";
  nkm::SurfData sdbuild(buildfile, Nvarsr, 0, 1, 0, 0, 0);
  nkm::SurfData sdeval(evalfile_in, Nvarsr, 0, 1, 0, 0, 0);  
  int Neval=sdeval.getNPts();
  nkm::KrigingModel km(sdbuild,km_params); km.create();
  nkm::MtxDbl y(Neval,1);
  nkm::MtxDbl vary(Neval,1);
  km.evaluate(y,sdeval.xr);
  km.eval_variance(vary,sdeval.xr);
  FILE* fp=fopen("evalfile.out","w");
  for(int i=0; i<Neval; ++i) {
    for(int j=0; j<Nvarsr; ++j)
      fprintf(fp,"%22.16g ",sdeval.xr(i,j));
    fprintf(fp,"%22.16g %22.16g\n",y(i,0),vary(i,0));
  }
  fclose(fp);
  return;
}

//used for quick develop/test of matlab implementation of adaptive importance sampling
void local_corrlen_build_dont_eval() {
  //the point of this is to a local optimization (starting from a good initial guess)
  //and have the calling function grab the resulting correlation lengths out of the 
  //command line model summary
  int Nvarsr=__GPAIS_NDIM__;
  ifstream infile("corrlen.txt",ios::in);
  string corrlen_str;
  getline(infile,corrlen_str);
  std::map< std::string, std::string> km_params;
  km_params["constraint_type"] = "r";
  km_params["order"] = "0";
  //km_params["reduced_polynomial"]=nkm::toString<bool>(true);
  km_params["optimization_method"]="local";
  km_params["correlation_lengths"]=corrlen_str;
#ifdef __CORR_FUNC0__
  km_params[__CORR_FUNC1__]=__CORR_FUNC2__;
#endif  

  string buildfile="buildfile.spd";
  nkm::SurfData sdbuild(buildfile, Nvarsr, 0, 1, 0, 0, 0);
  nkm::KrigingModel km(sdbuild,km_params); km.create();
  return;
}

//used for quick develop/test of matlab implementation of adaptive importance sampling
void local_corrlen_build_and_eval_mean_adjvar() {
  int Nvarsr=__GPAIS_NDIM__;
  ifstream infile("corrlen.txt",ios::in);
  string corrlen_str;
  getline(infile,corrlen_str);
  std::map< std::string, std::string> km_params;
  km_params["constraint_type"] = "r";
  km_params["order"] = "0";
  //km_params["reduced_polynomial"]=nkm::toString<bool>(true);
  km_params["optimization_method"]="local";
  km_params["correlation_lengths"]=corrlen_str;
#ifdef __CORR_FUNC0__
  km_params[__CORR_FUNC1__]=__CORR_FUNC2__;
#endif  

  string buildfile="buildfile.spd";
  string evalfile_in ="evalfile_in.spd";
  //string evalfile_out  ="evalfile.out";
  nkm::SurfData sdbuild(buildfile, Nvarsr, 0, 1, 0, 0, 0);
  nkm::SurfData sdeval(evalfile_in, Nvarsr, 0, 1, 0, 0, 0);  
  int Neval=sdeval.getNPts();
  nkm::KrigingModel km(sdbuild,km_params); km.create();
  nkm::MtxDbl y(Neval,1);
  nkm::MtxDbl vary(Neval,1);
  km.evaluate(y,sdeval.xr);
  km.eval_variance(vary,sdeval.xr);
  FILE* fp=fopen("evalfile.out","w");
  for(int i=0; i<Neval; ++i) {
    for(int j=0; j<Nvarsr; ++j)
      fprintf(fp,"%22.16g ",sdeval.xr(i,j));
    fprintf(fp,"%22.16g %22.16g\n",y(i,0),vary(i,0));
  }
  fclose(fp);
  return;
}

//used for quick develop/test of matlab implementation of adaptive importance sampling
void corrlen_build_and_eval_mean_adjvar() {
  int Nvarsr=__GPAIS_NDIM__;
  ifstream infile("corrlen.txt",ios::in);
  string corrlen_str;
  getline(infile,corrlen_str);
  std::map< std::string, std::string> km_params;
  km_params["constraint_type"] = "r";
  km_params["order"] = "0";
  //km_params["reduced_polynomial"]=nkm::toString<bool>(true);
  km_params["optimization_method"]="none";
  km_params["correlation_lengths"]=corrlen_str;
#ifdef __CORR_FUNC0__
  km_params[__CORR_FUNC1__]=__CORR_FUNC2__;
#endif  

  string buildfile="buildfile.spd";
  string evalfile_in ="evalfile_in.spd";
  //string evalfile_out  ="evalfile.out";
  nkm::SurfData sdbuild(buildfile, Nvarsr, 0, 1, 0, 0, 0);
  nkm::SurfData sdeval(evalfile_in, Nvarsr, 0, 1, 0, 0, 0);  
  int Neval=sdeval.getNPts();
  nkm::KrigingModel km(sdbuild,km_params); km.create();
  nkm::MtxDbl y(Neval,1);
  nkm::MtxDbl vary(Neval,1);
  km.evaluate(y,sdeval.xr);
  km.eval_variance(vary,sdeval.xr);
  FILE* fp=fopen("evalfile.out","w");
  for(int i=0; i<Neval; ++i) {
    for(int j=0; j<Nvarsr; ++j)
      fprintf(fp,"%22.16g ",sdeval.xr(i,j));
    fprintf(fp,"%22.16g %22.16g\n",y(i,0),vary(i,0));
  }
  fclose(fp);
  return;
}
