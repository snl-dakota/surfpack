// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines

class Surface;
class SurfData;


#include <iostream>
#include "SurfpackParser.h"

class SurfpackInterpreter
{
public:
  SurfpackInterpreter();
  ~SurfpackInterpereter();
  void execute(std::istream& is = cin, std::ostream& os = cout, 
    std::ostream& es = cerr);
  void commandLoop(std::ostream& os = cout, std::ostream& es = cerr);

  // individual surfpack commands
  executeLoadSurface(SurfpackParser::ParsedCommand& command);
  executeCreateSurface(SurfpackParser::ParsedCommand& command);
  
  
private:
  struct SymbolTable
  {
    std::map<std::string, SurfData*> dataVars;
    std::map<std::string, Surface*> surfaceVars;
  };

  struct SymbolTable symbol_table;
  SurfpackParser& parser;
};
    
