// These lines were automatically prepended
#ifndef SURFPACK_CONFIG_H
#define SURFPACK_CONFIG_H
#include "config.h"
#endif
// End of prepended lines
  
#include "SurfpackParser.h"

SurfpackInterpreter::SurfpackInterpreter() : parser(SurfpackParser::instance())
{

}

SurfpackInterpreter::~SurfpackInterpereter()
{

}

void SurfpackInterpreter::execute(std::istream& is, std::ostream& os, 
    std::ostream& es)
{
  if (parser.yyparse(&is, &os)) {
    commandLoop(os, es);
  } else {
    es << "Parse error.  Command(s) not executed." << endl;
  }
}

void SurfpackInterpreter::commandLoop(std::ostream& os, std::ostream& es)
{
  const vector<SurfpackParser::ParsedCommand>& commands = parser.commandList();
  for (unsigned i = 0; i < commands.size(); i++) {
    if (commands[i].name == "LoadSurface") {
      executeLoadSurface(commands[i]);
    else if (commands[i].name == "CreateSurface") {
      executeCreateSurface(commands[i]);
    else {
      es << "Unrecognized command: " << commands[i].name << endl;
    }
  }
}
  
void SurfpackInterpreter::executeLoadSurface(SurfpackParser::ParsedCommand& command)
{

}

void SurfpackInterpreter::executeCreateSurface(SurfpackParser::ParsedCommand& command)
{

}
