#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <surfpack.h>
using namespace std;
int main(int argc, char** argv)
{
  ifstream infile(argv[1], ios::in);
  string line;
  getline(infile, line);
  istringstream is(line);
  unsigned numvars;
  is >> numvars;
  vector<double> vars(numvars);
  for (unsigned i = 0; i < vars.size(); i++) {
    getline(infile, line);
    istringstream isv(line);
    isv >> vars[i];
  }
  infile.close();
  double rval = surfpack::rastrigin(vars);   
  //for (unsigned j = 0; j < vars.size(); j++) {
  //  cout << vars[j] << endl;
  //}

  ofstream outfile(argv[2], ios::out);
  outfile << rval << endl;
  outfile.close();
    
  return 0;
}
