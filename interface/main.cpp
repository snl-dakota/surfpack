#include "surfpack_config.h"
#include "surfpack_system_headers.h"
#include "surfpack.h"
#include "SurfpackInterpreter.h"

using namespace std;

int main(int argc, char** argv)
{
  // If a command line argument is given, us it as the input file
  // The default is to just use standard input
  SurfpackInterpreter si;
  if (argc == 2) {
    std::string infile(argv[1]);
    si.execute(&infile);
    //infile.close();
  } else {
    si.execute();
  }
  return 0;
}
