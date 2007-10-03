/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
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
  void execute(const std::string* input_string = 0, 
    const std::string* output_string = 0);
  void commandLoop(std::ostream& os = std::cout, std::ostream& es = std::cerr);

  // individual surfpack commands
  void executeConvertData(const ParsedCommand& command);
  void executeConvertSurface(const ParsedCommand& command);
  void executeCreateAxes(const ParsedCommand& command);
  void executeCreateSample(const ParsedCommand& command);
  void executeCreateSurface(const ParsedCommand& command);
  void executeEvaluate(const ParsedCommand& command);
  void executeFitness(const ParsedCommand& command);
  void executeLoad(const ParsedCommand& command);
  void executeLoadData(const ParsedCommand& command);
  void executeLoadSurface(const ParsedCommand& command);
  void executeSave(const ParsedCommand& command);
  void executeSaveData(const ParsedCommand& command);
  void executeSaveSurface(const ParsedCommand& command);
  void executeShellCommand(const ParsedCommand& command);
protected:
  class command_error 
  {
  public:
    command_error(const std::string& msg_in = "", 
      const std::string& cmdstring_in = "") 
      : msg(msg_in), cmdstring(cmdstring_in) {}
    void print() { std::cerr << "Error in " << cmdstring << ":  " 
			     << msg << std::endl; }
  protected:
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
    AxesBoundsMap axesVars;
    ~SymbolTable();  
    Surface* lookupSurface(const std::string);
    SurfData* lookupData(const std::string);
    AxesBounds* lookupAxes(const std::string);
  };

  struct SymbolTable symbolTable;
  SurfpackParser& parser;
};
    
