// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

#include <iostream>
#include <fstream>
#include "SurfpackInterpreter.h"

using namespace std;

int main(int argc, char** argv)
{
  // If a command line argument is given, us it as the input file
  // The default is to just use standard input
  SurfpackInterpreter si;
  if (argc == 2) {
    ifstream infile(argv[1], ios::in);
    si.execute(infile);
    infile.close();
  } else {
    si.execute();
  }
  return 0;
}
