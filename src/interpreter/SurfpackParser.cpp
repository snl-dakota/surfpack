#include <iostream>
#include <cstdlib>
#include "SurfpackParser.h"
#include <FlexLexer.h>

using namespace std;

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
}

void SurfpackParser::print()
{
  for (unsigned i = 0; i < commands.size(); i++) {
    cout << commands[i].name << endl;
    for (unsigned j = 0; j < commands[i].arglist.size(); j++) {
      cout << "   " << commands[i].arglist[j].name << " ";
      if (commands[i].arglist[j].lval.tuple.size() != 0) {
        for (unsigned k = 0; k < commands[i].arglist[j].lval.tuple.size(); k++) {
          cout << commands[i].arglist[j].lval.tuple[k] << " ";
        }
      } else if (commands[i].arglist[j].lval.triplet.numPts != 0) {
        cout << "{" 
             <<	commands[i].arglist[j].lval.triplet.min
             << ","
             <<	commands[i].arglist[j].lval.triplet.max
             << ","
             <<	commands[i].arglist[j].lval.triplet.numPts
             << "}" ;
      } 
      cout << endl;
    }
  }
}

void SurfpackParser::addCommandName()
{
  cout << "Add command name" << endl;
  ParsedCommand pc;
  commands.push_back(pc);
  currentArgList = &(commands[commands.size()-1].arglist);
  currentArgIndex = -1;
  commands[commands.size()-1].name = string(global_lexer.YYText());
}

void SurfpackParser::addArgName()
{
  cout << "Add arg name" << endl;
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
    (*currentArgList)[currentArgIndex].lval.identifier = string(global_lexer.YYText());
  }
}

void SurfpackParser::addArgValInt()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot assign Integer" << endl;
  } else {
    (*currentArgList)[currentArgIndex].lval.integer = atoi(global_lexer.YYText());
  }
}

void SurfpackParser::addArgValString()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot assign String" << endl;
  } else {
    (*currentArgList)[currentArgIndex].lval.literal= string(global_lexer.YYText());
  }
}

void SurfpackParser::addArgValReal()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot assign Real" << endl;
  } else {
    (*currentArgList)[currentArgIndex].lval.real = atof(global_lexer.YYText());
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
    currentArgList = &((*currentArgList)[currentArgIndex].lval.arglist);
  }
}

void SurfpackParser::addNumberAsTriplet()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot addNumberAsTriplet" << endl;
  } else {
    double val = atof(global_lexer.YYText());
    (*currentArgList)[currentArgIndex].lval.triplet.min = val;
    (*currentArgList)[currentArgIndex].lval.triplet.max = val;
    (*currentArgList)[currentArgIndex].lval.triplet.numPts= 1;
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
    (*currentArgList)[currentArgIndex].lval.triplet.min = val;
  }
}

void SurfpackParser::addTripletMax()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot addTripletMax" << endl;
  } else {
    double val = atof(global_lexer.YYText());
    (*currentArgList)[currentArgIndex].lval.triplet.max = val;
  }
}

void SurfpackParser::addTripletNumPts()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot addTripletMax" << endl;
  } else {
    int val = atoi(global_lexer.YYText());
    (*currentArgList)[currentArgIndex].lval.triplet.numPts = val;
  }
}

void SurfpackParser::addTupleVal()
{
  if (currentArgIndex == -1) {
    cerr << "currentArgIndex = -1; cannot addTupleVal" << endl;
  } else {
    double val = atof(global_lexer.YYText());
    (*currentArgList)[currentArgIndex].lval.tuple.push_back(val);
  }
}
