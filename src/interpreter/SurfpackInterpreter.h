#include "surfpack_config.h"

#include <iostream>
#include <map>
#include <set>
#include "SurfpackParser.h"
#include "SurfData.h"
#include "Surface.h"
#include "AxesBounds.h"


class SurfpackInterpreter
{
public:
  SurfpackInterpreter();
  ~SurfpackInterpreter();
  void execute(std::istream& is = std::cin, std::ostream& os = std::cout, 
    std::ostream& es = std::cerr);
  void commandLoop(std::ostream& os = std::cout, std::ostream& es = std::cerr);

  // individual surfpack commands
  void executeLoadData(const SurfpackParser::ParsedCommand& command);
  void executeLoadSurface(const SurfpackParser::ParsedCommand& command);
  void executeSaveData(const SurfpackParser::ParsedCommand& command);
  void executeSaveSurface(const SurfpackParser::ParsedCommand& command);
  void executeCreateSurface(const SurfpackParser::ParsedCommand& command);
  void executeConvertData(const SurfpackParser::ParsedCommand& command);
  void executeConvertSurface(const SurfpackParser::ParsedCommand& command);
  void executeEvaluate(const SurfpackParser::ParsedCommand& command);
  void executeFitness(const SurfpackParser::ParsedCommand& command);
  void executeAxesBounds(const SurfpackParser::ParsedCommand& command);
  void executeGridPoints(const SurfpackParser::ParsedCommand& command);
  void executeMonteCarloSample(const SurfpackParser::ParsedCommand& command);
  void executeShellCommand(const SurfpackParser::ParsedCommand& command);
  
  
protected:
  class command_error 
  {
  public:
    command_error(const std::string& msg_ = "", 
      const std::string& cmdstring_ = "") 
      : msg(msg_), cmdstring(cmdstring_) {}
    void print() { std::cerr << "Error in " << cmdstring << ":  " 
			     << msg << std::endl; }
  private:
    std::string msg;
    std::string cmdstring;
  };

public:
  typedef std::pair<std::string, SurfData*> SurfDataSymbol;
  typedef std::map<std::string, SurfData*> SurfDataMap;
  typedef std::pair<std::string, Surface*> SurfaceSymbol;
  typedef std::map<std::string, Surface*> SurfaceMap;
  typedef std::pair<std::string, AxesBounds*> AxesBoundsSymbol;
  typedef std::map<std::string, AxesBounds*> AxesBoundsMap;
  
  
private:
  struct SymbolTable
  {
    SurfDataMap dataVars;
    SurfaceMap surfaceVars;
    AxesBoundsMap pointDefinitionVars;
    ~SymbolTable() { 
      for (SurfDataMap::iterator iter = dataVars.begin();
	    iter != dataVars.end();
            ++iter) {
        delete iter->second; 
      }
      for (SurfaceMap::iterator siter = surfaceVars.begin();
	    siter != surfaceVars.end();
            ++siter) {
        delete siter->second; 
      }
      for (AxesBoundsMap::iterator pditer = pointDefinitionVars.begin();
	    pditer != pointDefinitionVars.end();
            ++pditer) {
        delete pditer->second; 
      }
    }
  };

  struct SymbolTable symbol_table;
  SurfpackParser& parser;
};
    
