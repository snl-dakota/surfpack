#include "surfpack_config.h"
#include "surfpack.h"
#include "AxesBounds.h"
#include "SurfpackParserArgs.h"
#include "SurfpackParser.h"
#include "SurfpackInterpreter.h"
#include "SurfData.h"
#include "Surface.h"
#include "SurfaceFactory.h"

using namespace std;

SurfpackInterpreter::SurfpackInterpreter() : parser(SurfpackParser::instance())
{

}

SurfpackInterpreter::~SurfpackInterpreter()
{

}

void SurfpackInterpreter::execute(std::istream& is, std::ostream& os, 
    std::ostream& es)
{
  if (parser.yyparse(is, os) == 0) {
    commandLoop(os, es);
  } else {
    es << "Parse error.  Command(s) not executed." << endl;
    es << SurfpackParser::cmdstream.str() << endl; 
  }
}

void SurfpackInterpreter::commandLoop(std::ostream& os, std::ostream& es)
{
  const vector<ParsedCommand>& commands = parser.commandList();
  for (unsigned i = 0; i < commands.size(); i++) {
    try {
      if (commands[i].isShellCommand()) {
	executeShellCommand(commands[i]);
      } else if (commands[i].name == "LoadSurface") {
        os << commands[i].cmdstring << endl;
        executeLoadSurface(commands[i]);
      } else if (commands[i].name == "LoadData") {
        executeLoadData(commands[i]);
      } else if (commands[i].name == "Save") {
        executeSave(commands[i]);
      } else if (commands[i].name == "SaveSurface") {
        executeSaveSurface(commands[i]);
      } else if (commands[i].name == "SaveData") {
        executeSaveData(commands[i]);
      } else if (commands[i].name == "ConvertData") {
        executeConvertData(commands[i]);
      } else if (commands[i].name == "ConvertSurface") {
        executeConvertSurface(commands[i]);
      } else if (commands[i].name == "CreateSurface") {
        executeCreateSurface(commands[i]);
      } else if (commands[i].name == "Evaluate") {
        executeEvaluate(commands[i]);
      } else if (commands[i].name == "Fitness") {
        executeFitness(commands[i]);
      } else if (commands[i].name == "CreateAxes") {
        executeCreateAxes(commands[i]);
      } else if (commands[i].name == "GridPoints") {
        executeGridPoints(commands[i]);
      } else if (commands[i].name == "MonteCarloSample") {
        executeMonteCarloSample(commands[i]);
      } else {
        es << "Unrecognized command: " << commands[i].name << endl;
      }
    } catch (command_error& ce) {
      ce.print();
    } catch (std::runtime_error& rte) {
      es << "FailedInstruction: " << commands[i].cmdstring << endl;
      es << rte.what() << endl;
    } catch (std::string& msg) {
      es << "FailedInstruction: " << commands[i].cmdstring << endl;
      es << msg << endl;
    } catch (...) {
      es << "FailedInstruction: " << commands[i].cmdstring << endl;
    }
  }
}
  
void SurfpackInterpreter::executeLoadData(
  const ParsedCommand& command)
{
  // Extract the variable name for this SurfData object
  string name = SurfpackParser::parseOutIdentifier(string("name"), command.arglist);
  if (name == "") {
    throw command_error(
      string("No name argument specified."), command.cmdstring);
  }
  //cout << "Variable name: " << name << endl;
  // Extract the name of the file to be read
  string filename = SurfpackParser::parseOutStringLiteral(string("file"), command.arglist);
  //cout << "Filename: " << filename << endl;
  if (filename == "") {
    throw command_error(
      string("No filename specified."), command.cmdstring);
  }
  // Read the file into a SurfData object and add it to the symbol table
  SurfData* sd = new SurfData(filename);
  symbol_table.dataVars.insert(SurfDataSymbol(name,sd));
  //cout << "Executed LoadData" << endl;
}

void SurfpackInterpreter::executeLoadSurface(
  const ParsedCommand& command)
{
  // Extract the variable name for this Surface object
  string name = SurfpackParser::parseOutIdentifier(string("name"), command.arglist);
  if (name == "") {
    throw command_error(
      string("No name argument specified."), command.cmdstring);
  }
  //cout << "Variable name: " << name << endl;
  // Extract the name of the file to be read
  string filename = SurfpackParser::parseOutStringLiteral(string("file"), command.arglist);
  //cout << "Filename: " << filename << endl;
  if (filename == "") {
    throw command_error(
      string("No filename specified."), command.cmdstring);
  }
  // Read the file into a Surface object and add it to the symbol table
  Surface* surf = SurfaceFactory::createSurface(filename);
  symbol_table.surfaceVars.insert(SurfaceSymbol(name,surf));
  
  //cout << "Executed LoadSurface" << endl;
}

void SurfpackInterpreter::executeSaveData(
  const ParsedCommand& command)
{
  // Extract the variable name for this SurfData object
  string data = SurfpackParser::parseOutIdentifier(string("data"), command.arglist);
  if (data == "") {
    throw command_error(
      string("No data argument specified."), command.cmdstring);
  }
  
  // Extract the name of the file to be written to 
  string filename = SurfpackParser::parseOutStringLiteral(string("file"), command.arglist);
  if (filename == "") {
    throw command_error(
      string("No filename specified."), command.cmdstring);
  }

  // Look up the data object in the symbol table
  SurfDataMap::iterator iter = symbol_table.dataVars.find(data);
  if (iter == symbol_table.dataVars.end()) {
    throw command_error(
      string("Symbol not found"), command.cmdstring);
  } else {
    // Data object found in symbol table, write the object to file
    SurfData* sd = iter->second;
    sd->write(filename);
  }
}

void SurfpackInterpreter::executeConvertData(
  const ParsedCommand& command)
{
  // Extract the variable name for this SurfData object
  string inputfile = SurfpackParser::parseOutStringLiteral(string("input"), command.arglist);
  if (inputfile == "") {
    throw command_error(
      string("No input filename specified."), command.cmdstring);
  }
  //cout << "Input name: " << inputfile << endl;
  // Extract the name of the file to be read
  string outputfile = SurfpackParser::parseOutStringLiteral(string("output"), command.arglist);
  //cout << "Output name: " << outputfile << endl;
  if (outputfile == "") {
    throw command_error(
      string("No output filename specified."), command.cmdstring);
  }
  SurfData* sd = new SurfData(inputfile);
  sd->write(outputfile);
  //cout << "Converted " << inputfile << " to " << outputfile << endl;
}
void SurfpackInterpreter::executeConvertSurface(
  const ParsedCommand& command)
{
  // Extract the variable name for this SurfData object
  string inputfile = SurfpackParser::parseOutStringLiteral(string("input"), command.arglist);
  if (inputfile == "") {
    throw command_error(
      string("No input filename specified."), command.cmdstring);
  }
  //cout << "Input name: " << inputfile << endl;
  // Extract the name of the file to be read
  string outputfile = SurfpackParser::parseOutStringLiteral(string("output"), command.arglist);
  //cout << "Output name: " << outputfile << endl;
  if (outputfile == "") {
    throw command_error(
      string("No output filename specified."), command.cmdstring);
  }
  Surface* surf = SurfaceFactory::createSurface(inputfile);
  surf->write(outputfile);
  //cout << "Converted " << inputfile << " to " << outputfile << endl;
}

void SurfpackInterpreter::executeSave(
  const ParsedCommand& command)
{
  // Extract the name of the file to be written to 
  string filename = SurfpackParser::parseOutStringLiteral(string("file"), command.arglist);
  //cout << "Filename: " << filename << endl;
  if (filename == "") {
    throw command_error(
      string("No filename specified."), command.cmdstring);
  }
  // Must have either surface or data arg
  // Extract the variable name for this Surface object
  string surface = SurfpackParser::parseOutIdentifier(string("surface"), command.arglist);
  string data = SurfpackParser::parseOutIdentifier(string("data"), command.arglist);
  if (surface == "") {
    if (data == "") {
      throw command_error(
        string("No surface or data argument specified."), command.cmdstring);
    } else {
      // save the data
      // Look up the data object in the symbol table
      SurfDataMap::iterator iter = symbol_table.dataVars.find(data);
      if (iter == symbol_table.dataVars.end()) {
        throw command_error(
          string("Data variable not found in symbol table"), command.cmdstring);
      } else {
        // Data object found in symbol table, write the object to file
        SurfData* sd = iter->second;
        sd->write(filename);
      }
    }
  } else if (data != "") {
    throw command_error(
      string("Cannot specify both data and surface."), command.cmdstring);
  } else {
    // save the surface
    SurfaceMap::iterator iter = symbol_table.surfaceVars.find(surface);
    if (iter == symbol_table.surfaceVars.end()) {
      ostringstream os;
      os << "Surface not found in symbol table: " << surface;
      throw command_error(os.str() , command.cmdstring);
    } else {
      Surface* surf = iter->second;
      surf->write(filename);
    }
  }
}
void SurfpackInterpreter::executeSaveSurface(
  const ParsedCommand& command)
{
  // Extract the variable name for this Surface object
  string surface = SurfpackParser::parseOutIdentifier(string("surface"), command.arglist);
  if (surface == "") {
    throw command_error(
      string("No surface argument specified."), command.cmdstring);
  }
  //cout << "Variable surface: " << surface << endl;
  // Extract the name of the file to be written to 
  string filename = SurfpackParser::parseOutStringLiteral(string("file"), command.arglist);
  //cout << "Filename: " << filename << endl;
  if (filename == "") {
    throw command_error(
      string("No filename specified."), command.cmdstring);
  }
  // Look up the surface in the symbol table 
  SurfaceMap::iterator iter = symbol_table.surfaceVars.find(surface);
  if (iter == symbol_table.surfaceVars.end()) {
    ostringstream os;
    os << "Surface not found in symbol table: " << surface;
    throw command_error(os.str() , command.cmdstring);
  } else {
    Surface* surf = iter->second;
    surf->write(filename);
  }
  //cout << "Executed Load Surface" << endl;
}

void SurfpackInterpreter::executeCreateSurface(const ParsedCommand& command)
{
  // Extract the variable name for this SurfData object
  string name = SurfpackParser::parseOutIdentifier(string("name"), command.arglist);
  if (name == "") {
    throw command_error(
      string("No name argument specified."), command.cmdstring);
  }
  //cout << "Surface Variable name: " << name << endl;

  // Extract the type of surface 
  string type = SurfpackParser::parseOutIdentifier(string("type"), command.arglist);
  if (type == "") {
    throw command_error(
      string("No surface type specified."), command.cmdstring);
  }
  //cout << "Surface type: " << type << endl;

  // Extract the name of the SurfData object (it should already be in the symbol table) 
  string dataName = SurfpackParser::parseOutIdentifier(string("data"), command.arglist);
  //cout << "Data name: " << dataName << endl;
  if (dataName == "") {
    throw command_error(
      string("No data object specified."), command.cmdstring);
  } else {
    SurfDataMap::iterator iter = symbol_table.dataVars.find(dataName);
    if (iter == symbol_table.dataVars.end()) {
      throw command_error(
        string("Data object not found"), command.cmdstring);
    }
    SurfData* sd = iter->second;
    Surface* surf = SurfaceFactory::createSurface(type, sd);
    surf->config(Arg(string("xsize"),
      new RvalInteger(static_cast<int>(sd->xSize()))));
    surf->configList(command.arglist);
    surf->createModel();
    symbol_table.surfaceVars.insert(SurfaceSymbol(name,surf));
  }
}

void SurfpackInterpreter::executeEvaluate(const ParsedCommand& command)
{
  // Extract the variable name for this SurfData object
  Surface* surf = 0;
  string surfaceName = SurfpackParser::parseOutIdentifier(string("surface"), command.arglist);
  if (surfaceName == "") {
    throw command_error(
      string("No existing surface specified."), command.cmdstring);
  } else {
    SurfaceMap::iterator iter = symbol_table.surfaceVars.find(surfaceName);
    if (iter == symbol_table.surfaceVars.end()) {
      throw command_error(
        string("Surface name not found."), command.cmdstring);
    } else {
      surf = iter->second;
    }
  }
  
  //cout << "Surface Variable name: " << surfaceName << endl;

  // Extract the name of the input data set (must be in the symbol table already) 
  string data = SurfpackParser::parseOutIdentifier(string("data"), command.arglist);
  SurfData* isd = 0;
  if (data == "") {
    throw command_error(
      string("No data specified."), command.cmdstring);
  } else {
    SurfDataMap::iterator iter = symbol_table.dataVars.find(data);
    if (iter == symbol_table.dataVars.end()) {
      throw command_error(
        string("Data object not found"), command.cmdstring);
    }
    isd = iter->second;
  }
  //cout << "data: " << data << endl;

  // Extract the name of the output data set 
  SurfData* osd = 0;
  string name = SurfpackParser::parseOutIdentifier(string("name"), command.arglist);
  if (name == "") {
    osd = isd;
  } else {
    //SurfDataMap::iterator iterForOutput = symbol_table.dataVars.find(dataName);
    //if (iterForOutput != symbol_table.dataVars.end()) {
    //  osd = iterForOutput->second;
    //} else {
      osd = new SurfData(*isd);
      symbol_table.dataVars.insert(SurfDataSymbol(name, osd));
    //}
    //cout << "name: " << name << endl;
  }

  // Now actually do the evaluation
  surf->getValue(*osd);
  
}

void SurfpackInterpreter::executeFitness(const ParsedCommand& command)
{
  // Extract the variable name for this Surface object
  Surface* surf = 0;
  string surfaceName = SurfpackParser::parseOutIdentifier(string("surface"), command.arglist);
  if (surfaceName == "") {
    throw command_error(
      string("No existing surface specified."), command.cmdstring);
  } else {
    SurfaceMap::iterator iter = symbol_table.surfaceVars.find(surfaceName);
    if (iter == symbol_table.surfaceVars.end()) {
      throw command_error(
        string("Surface name not found."), command.cmdstring);
    } else {
      surf = iter->second;
    }
  }
  
  //cout << "Surface Variable name: " << surfaceName << endl;

  // Extract the name of the SurfData object to be used (if any)
  string data = SurfpackParser::parseOutIdentifier(string("data"), command.arglist);
  SurfData* isd = 0;
  if (data != "") {
    SurfDataMap::iterator iter = symbol_table.dataVars.find(data);
    if (iter == symbol_table.dataVars.end()) {
      throw command_error(
        string("Data object not found"), command.cmdstring);
    }
    isd = iter->second;
  }
  //cout << "data: " << data << endl;

  // Extract the name of the output data set 
  string metric = SurfpackParser::parseOutIdentifier(string("metric"), command.arglist);
  if (metric == "") {
      throw command_error(
        string("No fitness metric specified."), command.cmdstring);
  }

  // Extract the response index for the input data set
  bool has_response_index = false;
  int response_index = SurfpackParser::parseOutInteger(string("response_index"), 
    command.arglist, has_response_index);
  if (has_response_index) {
    isd->setDefaultIndex(response_index);
  }

  double fitness = surf->goodnessOfFit(metric, isd);
  cout << metric << " for " << surfaceName;
  if (data != "") cout << " on " << data;
  cout << ": " << fitness << endl;
  
}
  
void SurfpackInterpreter::executeCreateAxes(const ParsedCommand& command)
{
  AxesBounds* sd = 0;
  // Extract the variable name for this AxesBounds object
  string name = 
    SurfpackParser::parseOutIdentifier(string("name"), command.arglist);
  if (name == "") {
    throw command_error(
      string("No name argument specified."), command.cmdstring);
  }

  //cout << "Variable name: " << name << endl;
  string bounds = 
    SurfpackParser::parseOutStringLiteral(string("bounds"), command.arglist);
  if (bounds == "") {
    // Extract the name of the file to be read
    string filename = 
      SurfpackParser::parseOutStringLiteral(string("file"), command.arglist);
    //cout << "Filename: " << filename << endl;
    if (filename == "") {
      throw command_error(
        string("No filename specified."), command.cmdstring);
    }
    sd = new AxesBounds(filename);
  } else {
    sd = new AxesBounds(bounds,AxesBounds::data);
  }
  // Read the file into a SurfData object and add it to the symbol table
  symbol_table.pointDefinitionVars.insert(AxesBoundsSymbol(name,sd));
  //cout << "Executed AxesBounds" << endl;
}

void SurfpackInterpreter::executeGridPoints(const ParsedCommand& command)
{
  // Extract the variable name for this SurfData object
  string axes = SurfpackParser::parseOutIdentifier(string("axes"), command.arglist);
  if (axes == "") {
    throw command_error(
      string("No axes argument specified."), command.cmdstring);
  }
  AxesBoundsMap::iterator iter = symbol_table.pointDefinitionVars.find(axes);
  if (iter == symbol_table.pointDefinitionVars.end()) {
    throw command_error(
      string("Definition not found in symbol table."), command.cmdstring);
  }
  //cout << "Variable axes: " << axes << endl;
  AxesBounds* pd = iter->second;

  // Extract the name of the SurfData object (it may already be in the symbol table) 
  string dataName = SurfpackParser::parseOutIdentifier(string("name"), command.arglist);
  //cout << "Data name: " << dataName << endl;
  if (dataName == "") {
    throw command_error(
      string("No data object specified."), command.cmdstring);
  } 
  SurfDataMap::iterator sditer = symbol_table.dataVars.find(dataName);
  if (sditer != symbol_table.dataVars.end()) {
  // Delete the one that's there because we're going to replace it
    delete sditer->second;
    symbol_table.dataVars.erase(sditer);
  } 
  vector<string> testFunctions;
  const ArgList& args = command.arglist;
  for (unsigned i = 0; i < args.size(); i++) {
    if (args[i].name == "test_function") {
      testFunctions.push_back(args[i].getRVal()->getIdentifier());
    }
  }
  SurfData* gridData = pd->sampleGrid(testFunctions);
  symbol_table.dataVars.insert(SurfDataSymbol(dataName, gridData));
  //cout << "Executed GridPoints" << endl;
}

void SurfpackInterpreter::executeMonteCarloSample(const ParsedCommand& command)
{
  // Extract the variable name for this SurfData object
  string axes = SurfpackParser::parseOutIdentifier(string("axes"), command.arglist);
  if (axes == "") {
    throw command_error(
      string("No axes argument specified."), command.cmdstring);
  }
  AxesBoundsMap::iterator iter = symbol_table.pointDefinitionVars.find(axes);
  if (iter == symbol_table.pointDefinitionVars.end()) {
    throw command_error(
      string("Definition not found in symbol table."), command.cmdstring);
  }
  //cout << "Variable axes: " << axes << endl;
  AxesBounds* pd = iter->second;

  // Extract the name of the SurfData object (it may already be in the symbol table) 
  string name = SurfpackParser::parseOutIdentifier(string("name"), command.arglist);
  //cout << "Data name: " << name << endl;
  if (name == "") {
    throw command_error(
      string("No name for resulting data object specified."), command.cmdstring);
  } 
  SurfDataMap::iterator sditer = symbol_table.dataVars.find(name);
  if (sditer != symbol_table.dataVars.end()) {
  // Delete the one that's there because we're going to replace it
    delete sditer->second;
    symbol_table.dataVars.erase(sditer);
  } 
  vector<string> testFunctions;
  unsigned numSamples = 100; // default value
  const ArgList& args = command.arglist;
  for (unsigned i = 0; i < args.size(); i++) {
    if (args[i].name == "test_function") {
      testFunctions.push_back(args[i].getRVal()->getIdentifier());
    } else if (args[i].name == "size") {
      numSamples = args[i].getRVal()->getInteger();
    }
  }
  SurfData* gridData = pd->sampleMonteCarlo(numSamples, testFunctions);
  symbol_table.dataVars.insert(SurfDataSymbol(name, gridData));
  //cout << "Executed MonteCarloSample" << endl;

}

void SurfpackInterpreter::
  executeShellCommand(const ParsedCommand& command)
{
  system(command.cmdstring.c_str());
}

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
  for (AxesBoundsMap::iterator pditer = pointDefinitionVars.begin();
        pditer != pointDefinitionVars.end();
        ++pditer) {
    delete pditer->second; 
  }
}
