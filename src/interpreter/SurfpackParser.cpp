#include "config.h"

#include <iostream>
#include <cstdlib>
#include "SurfpackParser.h"
//#include <FlexLexer.h>
#include <sstream>
#include <string>

using namespace std;
extern ostringstream cmdstream;

SurfpackParser& SurfpackParser::instance()
{
  static SurfpackParser sp;
  return sp;
}

yyFlexLexer& SurfpackParser::globalLexer()
{
  return global_lexer;
}

int SurfpackParser::yyparse(istream* is, ostream* os)
{
  global_lexer.switch_streams(is,os);
  return ::yyparse();
}

SurfpackParser::SurfpackParser()
{
  init();
}

void SurfpackParser::init()
{
  commands.clear();
  currentArgList = 0;
  currentArgIndex = -1;
  currentTupleIndex = -1;
  cmdstream.str("");
}

void SurfpackParser::storeCommandString()
{
  if (commands.size() > 0) {
    int loc;
    string newcommand = cmdstream.str();
    if (newcommand.find("/*") == 0) {
      newcommand.erase(0,2);
    }
    if ((loc = newcommand.find("*/")) != string::npos) {
      newcommand.erase(loc,2);
    }
    if (newcommand.find("!") == 0) {
       newcommand.erase(0,1);
    }
    commands[commands.size()-1].cmdstring = newcommand;
    cmdstream.str("");
  }
}

std::vector<SurfpackParser::ParsedCommand>& SurfpackParser::commandList()
{
  return commands;
}

void SurfpackParser::print()
{
  for (unsigned i = 0; i < commands.size(); i++) {
    cout << commands[i].name << endl;
    for (unsigned j = 0; j < commands[i].arglist.size(); j++) {
      cout << "   " << commands[i].arglist[j].name << " ";
      if (commands[i].arglist[j].rval.tuple.size() != 0) {
        for (unsigned k = 0; k < commands[i].arglist[j].rval.tuple.size(); k++) {
          cout << commands[i].arglist[j].rval.tuple[k] << " ";
        }
      } else if (commands[i].arglist[j].rval.triplet.numPts != 0) {
        cout << "{" 
             <<	commands[i].arglist[j].rval.triplet.min
             << ","
             <<	commands[i].arglist[j].rval.triplet.max
             << ","
             <<	commands[i].arglist[j].rval.triplet.numPts
             << "}" ;
      } 
      cout << endl;
    }
  }
}

void SurfpackParser::addCommandName()
{
  //cout << "Add command name" << endl;
  ParsedCommand pc;
  commands.push_back(pc);
  currentArgList = &(commands[commands.size()-1].arglist);
  currentArgIndex = -1;
  commands[commands.size()-1].name = string(global_lexer.YYText());
}

void SurfpackParser::addArgName()
{
  //cout << "Add arg name" << endl;
  if (currentArgList == 0) {
    cerr << "currentArgList is NULL; cannot assign name" << endl;
  } else {
    Arg newArg;
    currentArgList->push_back(newArg);
    currentArgIndex++;
    (*currentArgList)[currentArgIndex].name = string(global_lexer.YYText());
  }
}

void SurfpackParser::addArgValIdent()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot assign Identifier" << endl;
  } else {
    (*currentArgList)[currentArgIndex].rval.identifier = string(global_lexer.YYText());
  }
}

void SurfpackParser::addArgValInt()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot assign Integer" << endl;
  } else {
    (*currentArgList)[currentArgIndex].rval.integer = atoi(global_lexer.YYText());
  }
}

void SurfpackParser::addArgValString()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot assign String" << endl;
  } else {
    string currentToken = string(global_lexer.YYText());
    // The token contains leading and trailing apostrophes; remove them.
    int pos;
    while ( (pos = currentToken.find('\'')) != string::npos) {
      currentToken.erase(pos,pos+1);
    }
    (*currentArgList)[currentArgIndex].rval.literal= currentToken;
    //cout << "Stripped string: " << currentToken << endl;
  }
}

void SurfpackParser::addArgValReal()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot assign Real" << endl;
  } else {
    (*currentArgList)[currentArgIndex].rval.real = atof(global_lexer.YYText());
  }
}

void SurfpackParser::addArgValTuple()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot assign Tuple" << endl;
  } else {
    currentTupleIndex = -1;
  }
}

void SurfpackParser::addArgValArgList()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot assign ArgList" << endl;
  } else {
    currentArgList = &((*currentArgList)[currentArgIndex].rval.arglist);
  }
}

void SurfpackParser::addNumberAsTriplet()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot addNumberAsTriplet" << endl;
  } else {
    double val = atof(global_lexer.YYText());
    (*currentArgList)[currentArgIndex].rval.triplet.min = val;
    (*currentArgList)[currentArgIndex].rval.triplet.max = val;
    (*currentArgList)[currentArgIndex].rval.triplet.numPts= 1;
  }
}

void SurfpackParser::addTriplet()
{

}

void SurfpackParser::addTripletMin()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot addTripletMin" << endl;
  } else {
    double val = atof(global_lexer.YYText());
    (*currentArgList)[currentArgIndex].rval.triplet.min = val;
  }
}

void SurfpackParser::addTripletMax()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot addTripletMax" << endl;
  } else {
    double val = atof(global_lexer.YYText());
    (*currentArgList)[currentArgIndex].rval.triplet.max = val;
  }
}

void SurfpackParser::addTripletNumPts()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot addTripletMax" << endl;
  } else {
    int val = atoi(global_lexer.YYText());
    (*currentArgList)[currentArgIndex].rval.triplet.numPts = val;
  }
}

void SurfpackParser::addTupleVal()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot addTupleVal" << endl;
  } else {
    double val = atof(global_lexer.YYText());
    (*currentArgList)[currentArgIndex].rval.tuple.push_back(val);
  }
}

std::string SurfpackParser::parseOutIdentifier(const string& argname,
  const ArgList& arglist)
{
  for (unsigned i = 0; i < arglist.size(); i++) {
    if (arglist[i].name == argname) {
      return arglist[i].rval.identifier;
    }
  }
  return string("");
}

std::string SurfpackParser::parseOutStringLiteral(const string& argname,
  const ArgList& arglist)
{
  for (unsigned i = 0; i < arglist.size(); i++) {
    if (arglist[i].name == argname) {
      return arglist[i].rval.literal;
    }
  }
  return string("");
}

int SurfpackParser::parseOutInteger(const string& argname,
  const ArgList& arglist, bool& valid)
{
  valid = false;
  for (unsigned i = 0; i < arglist.size(); i++) {
    if (arglist[i].name == argname) {
      valid = true;
      return arglist[i].rval.integer;
    }
  }
  return -1;
}

void SurfpackParser::shellCommand()
{
  ParsedCommand pc(true); // ctor arg specifies that it's a shell command
  commands.push_back(pc);
  storeCommandString();
}
