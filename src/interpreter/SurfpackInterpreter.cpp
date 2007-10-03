/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#include "surfpack.h"
#include "AxesBounds.h"
#include "SurfpackParserArgs.h"
#include "SurfpackParser.h"
#include "SurfpackInterface.h"
#include "SurfpackInterpreter.h"
#include "SurfData.h"
#include "Surface.h"
#include "SurfaceFactory.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::runtime_error;
using std::string;
using std::vector;
using std::ostringstream;
using SurfpackInterface::CreateAxes;
using SurfpackInterface::CreateSurface;
using SurfpackInterface::Fitness;
using SurfpackInterface::Load;
using SurfpackInterface::Save;
///////////////////////////////////////////////////////////////////////////////
/////		SurfpackInterpreter namespace functions			  /////
///////////////////////////////////////////////////////////////////////////////

SurfpackInterpreter::SurfpackInterpreter() : parser(SurfpackParser::instance())
{

}

SurfpackInterpreter::~SurfpackInterpreter()
{

}

void SurfpackInterpreter::execute(const string* input_string, 
  const string* output_string)
{
  if (parser.yyparse(input_string, output_string) == 0) {
    commandLoop(cout , cerr);
  } else {
    cerr << "Command(s) not executed." << endl;
    cerr << SurfpackParser::cmdstream.str() << endl; 
  }
}

void SurfpackInterpreter::commandLoop(ostream& os, ostream& es)
{
  const vector<ParsedCommand>& commands = parser.commandList();
  for (unsigned i = 0; i < commands.size(); i++) {
    try {
      if (commands[i].isShellCommand()) {
	executeShellCommand(commands[i]);
      //} else if (commands[i].name == "ConvertData") {
      //  executeConvertData(commands[i]);
      //} else if (commands[i].name == "ConvertSurface") {
      //  executeConvertSurface(commands[i]);
      } else if (commands[i].name == "CreateAxes") {
        executeCreateAxes(commands[i]);
      } else if (commands[i].name == "CreateSample") {
        executeCreateSample(commands[i]);
      } else if (commands[i].name == "CreateSurface") {
        executeCreateSurface(commands[i]);
      } else if (commands[i].name == "Evaluate") {
        executeEvaluate(commands[i]);
      } else if (commands[i].name == "Fitness") {
        executeFitness(commands[i]);
      } else if (commands[i].name == "Load") {
        executeLoad(commands[i]);
      //} else if (commands[i].name == "LoadData") {
      //  executeLoadData(commands[i]);
      //} else if (commands[i].name == "LoadSurface") {
      //  executeLoadSurface(commands[i]);
      } else if (commands[i].name == "Save") {
        executeSave(commands[i]);
      //} else if (commands[i].name == "SaveData") {
      //  executeSaveData(commands[i]);
      //} else if (commands[i].name == "SaveSurface") {
      //  executeSaveSurface(commands[i]);
      } else {
        es << "Unrecognized command: " << commands[i].name << endl;
      }
    } catch (command_error& ce) {
      ce.print();
    } catch (runtime_error& rte) {
      es << "FailedInstruction: " << commands[i].cmdstring << endl;
      es << rte.what() << endl;
    } catch (string& msg) {
      es << "FailedInstruction: " << commands[i].cmdstring << endl;
      es << msg << endl;
    } catch (...) {
      es << "FailedInstruction: " << commands[i].cmdstring << endl;
    }
  }
}
  
void SurfpackInterpreter::executeLoad(const ParsedCommand& c)
{
  string filename = SurfpackParser::parseStringLiteral("file",c.arglist);
  if (surfpack::hasExtension(filename,".sps")) {
    executeLoadSurface(c);
  } else if (surfpack::hasExtension(filename,".spd") ||
	     surfpack::hasExtension(filename,".dat")) {
    executeLoadData(c);
  } else {
    throw string("Expected file extension: .sps (surface) or .spd (data)");
  }
}

void SurfpackInterpreter::executeLoadData(const ParsedCommand& c)
{
  SurfData* data = 0;
  string name = SurfpackParser::parseIdentifier("name",c.arglist);
  string filename = SurfpackParser::parseStringLiteral("file",c.arglist);
  bool valid = false;
  int n_vars = 
    SurfpackParser::parseInteger("n_predictors", c.arglist,valid,false);
  if (valid) {
    // n_responses is required if n_vars is present
    int n_responses = 
      SurfpackParser::parseInteger("n_responses", c.arglist,valid,true);
    int n_cols_to_skip = 
      SurfpackParser::parseInteger("n_cols_to_skip", c.arglist,valid,false);
    if (!valid) n_cols_to_skip = 0;
    Load(data,filename,n_vars,n_responses,n_cols_to_skip);
  } else {
    Load(data,filename);
  }
  assert(data);
  symbolTable.dataVars.insert(SurfDataSymbol(name,data));
}

void SurfpackInterpreter::executeLoadSurface(const ParsedCommand& c)
{
  string name = SurfpackParser::parseIdentifier("name",c.arglist);
  string filename = SurfpackParser::parseStringLiteral("file",c.arglist);
  Surface* surf = SurfaceFactory::createSurface(filename);
  symbolTable.surfaceVars.insert(SurfaceSymbol(name,surf));
}

void SurfpackInterpreter::executeSaveData(const ParsedCommand& c)
{
  string data_name = SurfpackParser::parseIdentifier("data", c.arglist);
  string filename = SurfpackParser::parseStringLiteral("file", c.arglist);
  SurfData* sd = symbolTable.lookupData(data_name);
  // Call SaveData 
  Save(sd,filename);
}

void SurfpackInterpreter::executeConvertData(
  const ParsedCommand& c)
{
}

void SurfpackInterpreter::executeConvertSurface(
  const ParsedCommand& c)
{
}

void SurfpackInterpreter::executeSave(const ParsedCommand& c)
{
  string filename = SurfpackParser::parseStringLiteral("file", c.arglist);
  // Don't automatically fail if either of these isn't defined
  string data_name = SurfpackParser::parseIdentifier("data", c.arglist,false);
  string surf_name = 
    SurfpackParser::parseIdentifier("surface", c.arglist, false);
  if (data_name == "") {
    if (surf_name == "") {
      // Do fail if both are missing
      throw command_error(
	string("Save command requires either 'surface' or 'data' argument"), 
        c.cmdstring);
    } else {
      Surface* surface = symbolTable.lookupSurface(surf_name);
      surface->write(filename);
    }
  } else {
    if (surf_name != "") {
      // Fail if both are specified 
      throw command_error(
	string("Save command may not have both 'surface' and 'data' arguments"),
        c.cmdstring);
    } else {
      SurfData* sd = symbolTable.lookupData(data_name);
      sd->write(filename);
    }
  }
}

void SurfpackInterpreter::executeSaveSurface(
  const ParsedCommand& c)
{
  string surf_name = SurfpackParser::parseIdentifier("surface", c.arglist);
  string filename = SurfpackParser::parseStringLiteral("file", c.arglist);
  Surface* surface = symbolTable.lookupSurface(surf_name);
  Save(surface,filename);
  // Call save surface
}

int getResponseIndex(const ArgList& arglist, const SurfData& sd)
{
  bool valid;
  bool is_response;
  unsigned int response_index;
  string response_name = SurfpackParser::parseIdentifier("response",
	arglist,false);
  if (response_name == "") {
      response_index = 
      SurfpackParser::parseInteger("response_index",arglist,valid,false);
      if (!valid) {
        // Neither response nor response_index specified
        // Default index is 0
        return 0;
      } else {
        return response_index;
      }
  } else {
    ostringstream os;
    valid = sd.varIndex(response_name, response_index, is_response);
    if (!valid) {
      os << "No response named '" << response_name << "' found." << endl;
      throw(os.str());
    } else if (!is_response) {
      os << "'" << response_name << "' is a predictor variable, but a"
	   << " response variable was requested" << endl;
      throw(os.str());
    } else {
      return response_index;
    }
  }
}

void SurfpackInterpreter::executeCreateSurface(const ParsedCommand& c)
{
  // Extract the variable name for this SurfData object
  string name = SurfpackParser::parseIdentifier("name", c.arglist);
  string data = SurfpackParser::parseIdentifier("data", c.arglist);
  string type = SurfpackParser::parseIdentifier("type", c.arglist);
  bool valid = false;
  SurfData* sd = symbolTable.lookupData(data);
  int response_index = getResponseIndex(c.arglist,*sd);
  // Call CreateSurface
  Surface* surface = 0; 
  CreateSurface(surface,sd,type,response_index);
  surface->configList(c.arglist);
  surface->createModel();
  symbolTable.surfaceVars.insert(SurfaceSymbol(name,surface));
}

void SurfpackInterpreter::executeCreateSample(const ParsedCommand& c)
{
  // Extract the variable name for this SurfData object
  string name = SurfpackParser::parseIdentifier("name", c.arglist);
  string axes = SurfpackParser::parseIdentifier("axes", c.arglist);
  vector<unsigned> grid_points = 
    SurfpackParser::parseUnsignedTuple("grid_points", c.arglist, false);
  vector<string> test_functions = 
    SurfpackParser::parseStringTuple("test_functions", c.arglist, false);
  vector<string> x_labels = 
    SurfpackParser::parseStringTuple("labels", c.arglist, false);
  bool valid = false;
  int n_points = SurfpackParser::parseInteger("size",c.arglist,valid,false);
  AxesBounds* ab = symbolTable.lookupAxes(axes);
  SurfData* sd = 0;
  if (valid) {
    if (!grid_points.empty()) { // both size and grid_points specified
      cerr << "Cannot specify both size and grid_points variables" << endl;
      exit(1);
    } else { // only size specified
	// MonteCarlo sample
      SurfpackInterface::CreateSample(sd,*ab,n_points,test_functions); 
    }
  } else {
    if (!grid_points.empty()) { // only grid_points specified
	// grid sample
      SurfpackInterface::CreateSample(sd,*ab,grid_points,test_functions); 
    } else { // neither specified
      cerr << "Must specify either size or grid_points" << endl;
      exit(1);
    }
  }
  if (!x_labels.empty()) sd->setXLabels(x_labels);
  symbolTable.dataVars.insert(SurfDataSymbol(name,sd));
}

void SurfpackInterpreter::executeEvaluate(const ParsedCommand& c)
{
  // Extract the variable name for this SurfData object
  Surface* surf = 0;
  string surf_name = SurfpackParser::parseIdentifier("surface", c.arglist);
  string data = SurfpackParser::parseIdentifier("data", c.arglist);
  Surface* surface = symbolTable.lookupSurface(surf_name);
  SurfData* sd = symbolTable.lookupData(data);
  // Call Evaluate
  unsigned new_index = surface->getValue(*sd);  
  string response_name = 
    SurfpackParser::parseIdentifier("label", c.arglist, false);
  if (response_name != "") {
    sd->setFLabel(new_index,response_name);
  }
}

void SurfpackInterpreter::executeFitness(const ParsedCommand& c)
{
  string surf_name = SurfpackParser::parseIdentifier("surface",c.arglist);
  string data = SurfpackParser::parseIdentifier("data", c.arglist, false);
  Surface* surface = symbolTable.lookupSurface(surf_name);
  SurfData* sd = 0;
  if (data != "") {
    sd = symbolTable.lookupData(data);
  }
  string metric = SurfpackParser::parseIdentifier("metric", c.arglist);
  // Extract the response index for the input data set
  bool valid = false;
  int response_index = 
    SurfpackParser::parseInteger("response_index", c.arglist,valid,false);
  if (!valid) { // No response_index was specified, use 0
    response_index = 0;
  }
  double fitness;
  if (metric == "cv") {
    unsigned n = 
      static_cast<unsigned>(SurfpackParser::parseInteger("n",c.arglist,valid));
    fitness = Fitness(surface,n,sd,response_index);
  } else {
    fitness = Fitness(surface,metric,sd,response_index); 
  }
  cout << metric << " for " << surf_name;
  if (data != "") cout << " on " << data;
  cout << ": " << fitness << endl;
}
  
void SurfpackInterpreter::executeCreateAxes(const ParsedCommand& c)
{
  string name = SurfpackParser::parseIdentifier("name", c.arglist);
  string bounds = 
    SurfpackParser::parseStringLiteral("bounds", c.arglist, true);
  CreateAxes(sd,bounds);
  }
  symbolTable.axesVars.insert(AxesBoundsSymbol(name,sd));
}

void SurfpackInterpreter::
  executeShellCommand(const ParsedCommand& c)
{
  system(c.cmdstring.c_str());
}

//********************************************************************
SurfpackInterpreter::SymbolTable::~SymbolTable() 
{ 
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
  for (AxesBoundsMap::iterator pditer = axesVars.begin();
        pditer != axesVars.end();
        ++pditer) {
    delete pditer->second; 
  }
}

Surface* SurfpackInterpreter::SymbolTable::lookupSurface(string name)
{
  SurfaceMap::iterator iter = surfaceVars.find(name);
  if (iter == surfaceVars.end()) {
    string msg = "Surface variable " + name + " not found in symbol table."; 
    throw msg;
  }
  assert(iter->second);
  return iter->second;
}

SurfData* SurfpackInterpreter::SymbolTable::lookupData(string name)
{
  SurfDataMap::iterator iter = dataVars.find(name);
  if (iter == dataVars.end()) {
    string msg = "Data variable " + name + " not found in symbol table."; 
    throw msg;
  }
  assert(iter->second);
  return iter->second;
}

AxesBounds* SurfpackInterpreter::SymbolTable::lookupAxes(string name)
{
  AxesBoundsMap::iterator iter = axesVars.find(name);
  if (iter == axesVars.end()) {
    string msg = "Axes variable " + name + " not found in symbol table."; 
    throw msg;
  }
  assert(iter->second);
  return iter->second;
}
