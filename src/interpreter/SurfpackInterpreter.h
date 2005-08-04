#include "surfpack_config.h"
#include "AxesBounds.h"
class SurfData;
class Surface;
class ParsedCommand;
class SurfpackParser;


class SurfpackInterpreter
{
public:
  SurfpackInterpreter();
  ~SurfpackInterpreter();
  void execute(std::istream& is = std::cin, std::ostream& os = std::cout, 
    std::ostream& es = std::cerr);
  void commandLoop(std::ostream& os = std::cout, std::ostream& es = std::cerr);

  // individual surfpack commands
  void executeLoadData(const ParsedCommand& command);
  void executeLoadSurface(const ParsedCommand& command);
  void executeSaveData(const ParsedCommand& command);
  void executeSaveSurface(const ParsedCommand& command);
  void executeCreateSurface(const ParsedCommand& command);
  void executeConvertData(const ParsedCommand& command);
  void executeConvertSurface(const ParsedCommand& command);
  void executeEvaluate(const ParsedCommand& command);
  void executeFitness(const ParsedCommand& command);
  void executeAxesBounds(const ParsedCommand& command);
  void executeGridPoints(const ParsedCommand& command);
  void executeMonteCarloSample(const ParsedCommand& command);
  void executeShellCommand(const ParsedCommand& command);
  
  
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
    ~SymbolTable();  
  };

  struct SymbolTable symbol_table;
  SurfpackParser& parser;
};
    
