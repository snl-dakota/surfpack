#include "NKM_SurfData.hpp"
#include "NKM_KrigingModel.hpp"
#include <iostream>

//#define __TIMING_BENCH__
using std::cout;
using std::endl;
using std::string;

using namespace nkm;

void validate();

int main(int argc, char* argv[])
{
  validate();
  return 0;
}


void validate()
{
  printf("validating Kriging Model\n");

  //filenames

#ifndef __TIMING_BENCH__  
  string validate2d_10 ="validate2d_10.spd";
  string validate2d_100="validate2d_100.spd";
  string validate2d_500="validate2d_500.spd";
  string validate2d_10K="validate2d_10K.spd";

  SurfData sd2d10( validate2d_10 , 2, 0, 3, 0, 0);
  SurfData sd2d100(validate2d_100, 2, 0, 3, 0, 0);
  SurfData sd2d500(validate2d_500, 2, 0, 3, 0, 0);
  SurfData sd2d10K(validate2d_10K, 2, 0, 3, 0, 0);
#endif

  string paviani10d_50  ="paviani10d_50.spd";
  string paviani10d_500 ="paviani10d_500.spd";
  string paviani10d_2500="paviani10d_2500.spd";
  string paviani10d_10K ="paviani10d_10K.spd";


  MtxDbl yeval10(    10);
#ifndef __TIMING_BENCH__
  MtxDbl yeval100(  100);
  MtxDbl yeval500(  500);
#endif
  MtxDbl yeval10K(10000);

  int jout;

  std::map< std::string, std::string> km_params;
  km_params["constraint_type"] = "r";
  //km_params["order"] = "linear";
  km_params["order"] = "2";
  km_params["reduced_polynomial"]=toString<bool>(true);


#ifndef __TIMING_BENCH__  
  printf("*****************************************************************\n");
  printf("*** running rosenbrock 2D tests *********************************\n");
  printf("*****************************************************************\n");

  MtxDbl roserror(3,4); roserror.zero();
  
  jout=0; //the 0th output column is Rosenbrock  
  sd2d10.setJOut( jout);
  sd2d100.setJOut(jout);
  sd2d500.setJOut(jout);
  sd2d10K.setJOut(jout);

  km_params["lower_bounds"]="-2.0 -2.0";
  km_params["upper_bounds"]="2.0 2.0";

  KrigingModel kmros10( sd2d10 , km_params); kmros10.create();
  KrigingModel kmros100(sd2d100, km_params); kmros100.create();
  KrigingModel kmros500(sd2d500, km_params); kmros500.create();

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

  MtxDbl shuerror(3,4); shuerror.zero();
  jout=1;
  sd2d10.setJOut( jout);
  sd2d100.setJOut(jout);
  sd2d500.setJOut(jout);
  sd2d10K.setJOut(jout);

  KrigingModel kmshu10( sd2d10 , km_params); kmshu10.create();
  KrigingModel kmshu100(sd2d100, km_params); kmshu100.create();
  KrigingModel kmshu500(sd2d500, km_params); kmshu500.create();


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

  MtxDbl herberror(3,4); herberror.zero();

  jout=2;
  sd2d10.setJOut( jout);
  sd2d100.setJOut(jout);
  sd2d500.setJOut(jout);
  sd2d10K.setJOut(jout);

  KrigingModel kmherb10( sd2d10 , km_params); kmherb10.create();
  KrigingModel kmherb100(sd2d100, km_params); kmherb100.create();
  KrigingModel kmherb500(sd2d500, km_params); kmherb500.create();

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


  MtxDbl paverror(3,4); paverror.zero();
#ifndef __TIMING_BENCH__
  SurfData sdpav50( paviani10d_50 , 10, 0, 1, 0, 0);
  SurfData sdpav500(paviani10d_500, 10, 0, 1, 0, 0);
#endif
  SurfData sdpav2500(paviani10d_2500, 10, 0, 1, 0, 0);
  SurfData sdpav10K(paviani10d_10K, 10, 0, 1, 0, 0);
#ifndef __TIMING_BENCH__
  KrigingModel kmpav50( sdpav50 , km_params); kmpav50.create();
  KrigingModel kmpav500(sdpav500, km_params); kmpav500.create();
#endif
  KrigingModel kmpav2500(sdpav2500, km_params); kmpav2500.create();
  //cout << kmpav2500.model_summary_string();

#ifndef __TIMING_BENCH__
  MtxDbl yeval50(50);
#endif
  MtxDbl yeval2500(2500);

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
      
  //evaluate error the 2500 pt paviani10d kriging model at 10K points
  kmpav2500.evaluate(yeval10K,sdpav10K.xr);
  for(int i=0; i<10000; ++i)
    paverror(2,2)+=pow(yeval10K(i)-sdpav10K.y(i),2);
  paverror(2,3)=sqrt(paverror(2,2)/10000.0);

#ifdef __TIMING_BENCH__  
  sdpav10K.y.copy(yeval10K);
  string pav10Kout="paviani10d_10K_nkm_out.spd";
  sdpav10K.write(pav10Kout);
#endif

  sdpav10K.clear();

#ifndef __TIMING_BENCH__  
  //evaluate error the 2500 pt paviani10d kriging model at build points
  kmpav2500.evaluate(yeval2500,sdpav2500.xr);
  for(int i=0; i<2500; ++i)
    paverror(2,0)+=pow(yeval2500(i)-sdpav2500.y(i),2);
  paverror(2,1)=sqrt(paverror(2,0)/2500.0);
  
  sdpav2500.clear();
  yeval2500.clear();
  
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
