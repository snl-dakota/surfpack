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
#include "SurfpackModel.h"

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
    //commandLoopOld(cout , cerr);
  } else {
    cerr << "Command(s) not executed." << endl;
    cerr << SurfpackParser::cmdstream.str() << endl; 
  }
}


void SurfpackInterpreter::commandLoopOld(ostream& os, ostream& es)
{
  cout << "Old Surfpack" << endl;
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

void SurfpackInterpreter::commandLoop(ostream& os, ostream& es)
{
  const vector<ParsedCommand>& fullCommands = parser.commandList();
  vector<Command>& commands = parser.comms;
  for (unsigned i = 0; i < commands.size(); i++) {
    try {
      //if (commands[i].isShellCommand()) {
      //  executeShellCommand(commands[i]);
      if (commands[i].first == "CreateSample") {
        execCreateSample(commands[i].second);
      } else if (commands[i].first== "CreateSurface") {
        execCreateSurface(commands[i].second);
      } else if (commands[i].first == "Evaluate") {
        execEvaluate(commands[i].second);
      } else if (commands[i].first == "Fitness") {
        execFitness(commands[i].second);
      } else if (commands[i].first == "Load") {
        execLoad(commands[i].second);
      } else if (commands[i].first == "Save") {
        execSave(commands[i].second);
      } else if (commands[i].first == "CreateAxes") {
        execCreateAxes(commands[i].second);
      } else {
        es << "Unrecognized command: " << commands[i].first << endl;
      }
    } catch (string& msg) {
      es << "FailedInstruction: " << fullCommands[i].cmdstring << endl;
      es << msg << endl;
    } catch (...) {
      es << "Unknown Exception.  FailedInstruction: " 
	 << fullCommands[i].cmdstring << endl;
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
  AxesBounds* ab;
  SurfpackInterface::CreateAxes(ab,bounds);
  symbolTable.axesVars.insert(AxesBoundsSymbol(name,ab));
}

void SurfpackInterpreter::
  executeShellCommand(const ParsedCommand& c)
{
  system(c.cmdstring.c_str());
}

///////////////////////////////////////////////////////////////////////////////
/////		SurfpackInterpreter namespace functions			  /////
///////////////////////////////////////////////////////////////////////////////

void SurfpackInterpreter::execLoad(ParamMap& args)
{
  bool valid;
  string filename = asStr(args["file"]);
  if (surfpack::hasExtension(filename,".sps")) {
    throw string("Read surface feature not currently supported.");
  //  executeLoadSurface(c);
  } else if (surfpack::hasExtension(filename,".spd") ||
	     surfpack::hasExtension(filename,".dat")) {
    execLoadData(args);
  } else {
    throw string("Expected file extension: .sps (surface) or .spd (data)");
  }
}

void SurfpackInterpreter::execLoadData(ParamMap& args)
{
  SurfData* data = 0;
  string name = asStr(args["name"]); 
  string filename = asStr(args["file"]); 
  bool valid = false;
  int n_vars = asInt(args["n_predictors"],valid);
  if (valid) {
    // n_responses is required if n_vars is present
    int n_responses = asInt(args["n_responses"]); // fail on error
    int n_cols_to_skip = asInt(args["n_cols_to_skip"],valid);
    if (!valid) n_cols_to_skip = 0;
    data = NewInterface::LoadData(filename,n_vars,n_responses,n_cols_to_skip);
  } else {
    data = NewInterface::LoadData(filename);
  }
  assert(data);
  symbolTable.dataVars.insert(SurfDataSymbol(name,data));
}

///\todo Add support for LoadSurface in interpreter
//void SurfpackInterpreter::execLoadSurface(ParamMap& args)
//{
//  string name = SurfpackParser::parseIdentifier("name",c.arglist);
//  string filename = SurfpackParser::parseStringLiteral("file",c.arglist);
//  Surface* surf = SurfaceFactory::createSurface(filename);
//  symbolTable.surfaceVars.insert(SurfaceSymbol(name,surf));
//}
//
void SurfpackInterpreter::execSaveData(ParamMap& args)
{
  string data_name = asStr(args["data"]);
  string filename = asStr(args["file"]);
  SurfData* sd = symbolTable.lookupData(data_name);
  // Call SaveData 
  NewInterface::Save(sd,filename);
}

void SurfpackInterpreter::execSave(ParamMap& args)
{
  try {
    string filename = asStr(args["file"]); 
    // Don't automatically fail if either of these isn't defined
    bool valid_data;
    string data_name = asStr(args["data"],valid_data); 
    bool valid_surface;
    string surf_name = asStr(args["surface"],valid_surface);
    if (!valid_data) {
      if (!valid_surface) {
        // Do fail if both are missing
        throw string("Save command requires either 'surface' or 'data' argument");
      } else {
        std::cout << "Model written" << std::endl;
        //Surface* surface = symbolTable.lookupSurface(surf_name);
        //surface->write(filename);
      }
    } else {
      if (valid_surface) {
        // Fail if both are specified 
        throw string("Save command may not have both 'surface' and 'data' arguments");
      } else {
        SurfData* sd = symbolTable.lookupData(data_name);
        sd->write(filename);
      }
    }
  } catch (string e) {
    cout << e << endl;
  } catch (...) {
    cout << "Caught some other error" << endl;
  }
}
//
//void SurfpackInterpreter::execSaveSurface(
//  ParamMap& args)
//{
//  string surf_name = SurfpackParser::parseIdentifier("surface", c.arglist);
//  string filename = SurfpackParser::parseStringLiteral("file", c.arglist);
//  Surface* surface = symbolTable.lookupSurface(surf_name);
//  Save(surface,filename);
//  // Call save surface
//}
//
//int getResponseIndex(const ArgList& arglist, const SurfData& sd)
//{
//  bool valid;
//  bool is_response;
//  unsigned int response_index;
//  string response_name = SurfpackParser::parseIdentifier("response",
//	arglist,false);
//  if (response_name == "") {
//      response_index = 
//      SurfpackParser::parseInteger("response_index",arglist,valid,false);
//      if (!valid) {
//        // Neither response nor response_index specified
//        // Default index is 0
//        return 0;
//      } else {
//        return response_index;
//      }
//  } else {
//    ostringstream os;
//    valid = sd.varIndex(response_name, response_index, is_response);
//    if (!valid) {
//      os << "No response named '" << response_name << "' found." << endl;
//      throw(os.str());
//    } else if (!is_response) {
//      os << "'" << response_name << "' is a predictor variable, but a"
//	   << " response variable was requested" << endl;
//      throw(os.str());
//    } else {
//      return response_index;
//    }
//  }
//}
//
void SurfpackInterpreter::execCreateSurface(ParamMap& args)
{
  // Extract the variable name for this SurfData object
  string name = asStr(args["name"]);
  string data = asStr(args["data"]);
  bool valid = false;
  SurfData* sd = symbolTable.lookupData(data);
  // Call CreateSurface
  SurfpackModelFactory* smf = SurfaceFactory::createModelFactory(args);
  SurfpackModel* model = smf->Build(*sd);
  assert(model);
  symbolTable.modelVars.insert(SurfpackModelSymbol(name,model));
}

void SurfpackInterpreter::execCreateSample(ParamMap& args)
{
  // Extract the variable name for this SurfData object
  string name = asStr(args["name"]);
  string axes = asStr(args["axes"]);
  bool valid_grid_points;
  VecUns grid_points = asVecUns(args["grid_points"], valid_grid_points);
  bool valid_size = false;
  int n_points = asInt(args["size"],valid_size); 
  AxesBounds* ab = symbolTable.lookupAxes(axes);
  SurfData* sd = 0;
  if (valid_size) {
    if (valid_grid_points) { // both size and grid_points specified
      throw string("Cannot specify both size and grid_points");
    } else { // only size specified
	// MonteCarlo sample
      sd = NewInterface::CreateSample(ab,n_points); 
    }
  } else {
    if (!grid_points.empty()) { // only grid_points specified
	// grid sample
      sd = NewInterface::CreateSample(ab,grid_points); 
    } else { // neither specified
      throw string("Must specify either size or grid_points");
    }
  }
  bool valid_functions;
  VecStr test_functions = asVecStr(args["test_functions"], valid_functions);
  if (valid_functions) {
    NewInterface::Evaluate(sd,test_functions);
  }
  symbolTable.dataVars.insert(SurfDataSymbol(name,sd));
}

void SurfpackInterpreter::execEvaluate(ParamMap& args)
{
  // Extract the variable name for this SurfData object
  string surf_name = asStr(args["surface"]);
  string data = asStr(args["data"]);
  SurfpackModel* model = symbolTable.lookupModel(surf_name);
  SurfData* sd = symbolTable.lookupData(data);
  // Call Evaluate
  VecDbl results = (*model)(*sd);
  sd->addResponse(results,args["label"]); 
}

void SurfpackInterpreter::execFitness(ParamMap& args)
{
  string surface = asStr(args["surface"]);
  bool valid_data;
  string data = asStr(args["data"],valid_data);
  SurfpackModel* model = symbolTable.lookupModel(surface);
  SurfData* sd = 0;
  if (valid_data) {
    sd = symbolTable.lookupData(data);
  }
  string metric = asStr(args["metric"]);
  // Extract the response index for the input data set
  bool valid_response_index;
  int response_index = asInt(args["response_index"],valid_response_index);
  bool valid_n;
  int n = asInt(args["n"],valid_n);
  if (!valid_response_index) { // No response_index was specified, use 0
    response_index = 0;
  }
  double fitness;
  if (valid_data) {
    fitness = NewInterface::Fitness(model,sd,metric,response_index,n);
  } else {
    fitness = NewInterface::Fitness(model,metric,response_index,n); 
  }
  cout << metric << " for " << surface;
  if (data != "") cout << " on " << data;
  cout << ": " << fitness << endl;
}
  
void SurfpackInterpreter::execCreateAxes(ParamMap& args)
{
  string name = asStr(args["name"]);
  string bounds = asStr(args["bounds"]);
  AxesBounds* ab = NewInterface::CreateAxes(bounds);
  symbolTable.axesVars.insert(AxesBoundsSymbol(name,ab));
}
//
//void SurfpackInterpreter::
//  executeShellCommand(ParamMap& args)
//{
//  system(c.cmdstring.c_str());
//}

int SurfpackInterpreter::asInt(const string& arg)
{
  if (arg == "") throw string("Expected integer value");
  return atoi(arg.c_str());
}

std::string SurfpackInterpreter::asStr(const string& arg)
{
  if (arg == "") throw string("Expected string value");
  return arg;
}

double SurfpackInterpreter::asDbl(const string& arg)
{
  if (arg == "") throw string("Expected double value");
  return atof(arg.c_str());
}

VecUns SurfpackInterpreter::asVecUns(const string& arg)
{
  if (arg == "") throw string("Expected vector unsigned");
  return surfpack::toVec<unsigned>(arg); 
}

VecStr SurfpackInterpreter::asVecStr(const string& arg)
{
  if (arg == "") throw string("Expected vector string");
  return surfpack::toVec<std::string>(arg); 
}

int SurfpackInterpreter::asInt(const string& arg, bool& valid)
{
  if (arg == "") {
    valid = false;
    return 0;
  }
  valid = true;
  return atoi(arg.c_str());
}

std::string SurfpackInterpreter::asStr(const string& arg, bool& valid)
{
  if (arg == "") {
    valid = false;
    return "";
  }
  valid = true;
  return arg; 
}

double SurfpackInterpreter::asDbl(const string& arg, bool& valid)
{
  if (arg == "") {
    valid = false;
    return 0;
  }
  valid = true;
  return atof(arg.c_str()); 
}

VecUns SurfpackInterpreter::asVecUns(const string& arg, bool& valid)
{
  if (arg == "") {
    valid = false;
    return VecUns();
  }
  valid = true;
  return surfpack::toVec<unsigned>(arg); 
}

VecStr SurfpackInterpreter::asVecStr(const string& arg, bool& valid)
{
  if (arg == "") {
    valid = false;
    return VecStr();
  }
  valid = true;
  return surfpack::toVec<std::string>(arg); 
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

SurfpackModel* 
SurfpackInterpreter::SymbolTable::lookupModel(const std::string name)
{
  SurfpackModel* result = modelVars[name];
  if (result == 0) {
    cout << "Bad lookup; table size:  " << modelVars.size() << endl;
    for (SurfpackModelMap::iterator itr = modelVars.begin();
        itr != modelVars.end(); ++itr) {
      cout << itr->first << " " << itr->second << endl;
    }
    
    string msg = "Model variable " + name + " not found in symbol table."; 
    throw msg;
  }
  return result;
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
