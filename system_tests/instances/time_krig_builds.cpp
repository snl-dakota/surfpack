#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#include "surfpack_system_headers.h"

#include "surfpack.h"
#include "SurfData.h"
#include "SurfpackInterpreter.h"
#include "AxesBounds.h"
#include "KrigingSurface.h"
#include "PolynomialSurface.h"

#if !defined(HAVE_GETTIMEOFDAY) && (defined(_MSC_VER) || defined(__MINGW32__))
#include <windows.h>
#endif

using namespace std;
using namespace SurfpackInterface;

// Modified from http://mywebpage.netscape.com/yongweiwu/timeval.h.txt
#if !defined(HAVE_GETTIMEOFDAY) && (defined(_MSC_VER) || defined(__MINGW32__))
int gettimeofday (struct timeval *tv, void* tz)
{
  union {
    __int64 ns100; /*time since 1 Jan 1601 in 100ns units */
    FILETIME ft;
  } now;

  GetSystemTimeAsFileTime (&now.ft);
  tv->tv_usec = (long) ((now.ns100 / 10LL) % 1000000LL);
  tv->tv_sec = (long) ((now.ns100 - 116444736000000000LL) / 10000000LL);
  return (0);
}
#endif

double time_difference(struct timeval& starttime, struct timeval& endtime)
{
  return (((double)endtime.tv_sec+(1.0e-06)*endtime.tv_usec) -
	((double)starttime.tv_sec+(1.0e-06)*starttime.tv_usec)) ;

}

int main(int argc, char** argv)
{
  AxesBounds* ab = 0;
  CreateAxes(ab,string("../system_tests/axes_bounds/generic_2d.axb"),
    AxesBounds::file);
  //AxesBounds ab(string("../system_tests/axes_bounds/generic_2d.axb"));
  vector<string> functions;
  vector<double> thetas(2,1.0);
  vector<unsigned> setsizes;
  vector<double> onepoint(1);
  const int num_trials = 1;
  vector<double> times_one_size(num_trials);
  vector<double> responses(2);
  double timeneeded;
  SurfData krigtimes;
  for (unsigned i = 3050; i <= 5000; i += 50) setsizes.push_back(i);
  setsizes.push_back(10000);
  //for (unsigned i = 500; i < 1000; i += 25) setsizes.push_back(i);
  //for (unsigned i = 1000; i <= 2500; i += 100) setsizes.push_back(i);
  functions.push_back(string("rosenbrock"));
  struct timeval t1;
  struct timeval t2;
  for (unsigned setsize = 0; setsize < setsizes.size(); setsize++) {
    cout << setw(8) << setsizes[setsize] ;
    for (unsigned trial = 0; trial < num_trials; trial++) {
      SurfData* sd = 0;
      CreateSample(sd,ab,setsizes[setsize],functions);
      //SurfData* sd = ab.sampleMonteCarlo(setsizes[setsize],functions);
      Surface* ks = 0;
      CreateSurface(ks,sd,string("kriging"));
      //KrigingSurface ks(sd);
      dynamic_cast<KrigingSurface*>(ks)->usePreComputedCorrelationVector(thetas);
      //ks.usePreComputedCorrelationVector(thetas);
      gettimeofday(&t1,NULL);
      ks->createModel();
      gettimeofday(&t2,NULL);
      timeneeded = time_difference(t1,t2); 
      times_one_size[trial] = timeneeded;
      cout <<  setw(15) << timeneeded;
      delete sd; sd = 0;
      delete ks; ks = 0;
    }
    sort(times_one_size.begin(),times_one_size.end());
    double avgtime = times_one_size[times_one_size.size()/2];
    double meantime = accumulate(times_one_size.begin(),times_one_size.end(),0.0)/times_one_size.size();
    onepoint[0] = setsizes[setsize];
    responses[0] = avgtime;
    responses[1] = meantime;
    krigtimes.addPoint(SurfPoint(onepoint,responses));
    cout << " avg: " << setw(15) 
         << avgtime << endl;
  }
  Save(&krigtimes,string("one_trial_3050to5000.spd"));
  //krigtimes.write(string("one_trial_3050to5000.txt"));
  Surface* ps = 0;
  CreateSurface(ps,&krigtimes,string("polynomial"));
  ///\todo configure to be order 3
  //PolynomialSurface ps(&krigtimes,3);

  //SurfData trainingtimes("trainingtimes.txt");
  //PolynomialSurface ps(&trainingtimes,3);
  ps->createModel();
  Save(ps,string("poly3_krigtimes.sps"));
  //ps.write("poly3_krigtimes.txt");
  // cleanup
  delete ab; ab = 0;
  return 0;
}
