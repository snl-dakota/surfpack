#include "config.h"

#ifndef MY_PARSING_STRUCTURES_H
#define MY_PARSING_STRUCTURES_H

#include <string>
#include <vector>
#include "FlexLexer.h"
#include "surfparse.tab.h"
extern int yyparse();

class SurfpackParser 
{
public:
  struct Arg;
  typedef std::vector<Arg> ArgList;
  typedef std::vector<double> Tuple;
  struct Triplet {
    Triplet() : min(0), max(0), numPts(0) {}
    double min;
    double max;
    unsigned numPts;
  };
  struct Lval
  {
    int integer;
    double real;
    Tuple tuple;
    Triplet triplet;
    std::string identifier;
    std::string literal;
    ArgList arglist;
  };
  struct Arg
  {
    std::string name;
    Lval lval;
  }; 
  struct ParsedCommand
  { 
    std::string name;
    ArgList arglist;
    std::string cmdstring;
  };

  // commands for use by the Bison parser
  void addCommandName();
  void addArgName();
  void addArgValIdent();
  void addArgValInt();
  void addArgValString();
  void addArgValReal();
  void addArgValTuple();
  void addArgValArgList();
  void addNumberAsTriplet();
  void addTriplet();
  void addTripletMin();
  void addTripletMax();
  void addTripletNumPts();
  void addTupleVal();
  void init();
  void print();
  void storeCommandString();


  static SurfpackParser& instance();
  yyFlexLexer& globalLexer();
  int yyparse(std::istream* is = &std::cin, std::ostream* os = &std::cout);
  std::vector<ParsedCommand>& commandList();

  // commands for use by interpreter to parse out certain argument types
  static std::string parseOutIdentifier(const std::string& argname,
    const ArgList& arglist);
  static std::string parseOutStringLiteral(const std::string& argname,
    const ArgList& arglist);
protected:
// Default constructor, copy constructor and assignment operator declared 
// protected and not implemented, in order to make sure that only one
// instance is created
  SurfpackParser();
  SurfpackParser(const SurfpackParser&);
  const SurfpackParser& operator=(const SurfpackParser&);

private:
  std::vector<ParsedCommand> commands;
  ArgList* currentArgList;
  int currentArgIndex;
  int currentTupleIndex;
  yyFlexLexer global_lexer;
};

#endif




