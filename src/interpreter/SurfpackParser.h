#include "config.h"

#ifndef MY_PARSING_STRUCTURES_H
#define MY_PARSING_STRUCTURES_H

#include <string>
#include <vector>
#include <FlexLexer.h>
#include "surfparse.h"
extern int yyparse();

class SurfpackParser 
{
public:
  class Arg;
  typedef std::vector<Arg> ArgList;
  typedef std::vector<double> Tuple;
  class Triplet 
  {
  public:
    Triplet() : min(0), max(0), numPts(0) {}
    double min;
    double max;
    unsigned numPts;
  };
  class Rval
  {
  public:
    int integer;
    double real;
    Tuple tuple;
    Triplet triplet;
    std::string identifier;
    std::string literal;
    ArgList arglist;
    Rval() {}
    Rval(int integer_in) : integer(integer_in) {}
    Rval(double real_in) : real(real_in) {}
    Rval(const Tuple& tuple_in) : tuple(tuple_in) {}
    Rval(const std::string& string_in, const std::string& type) 
    {
      if (type == "identifier") {
        this->identifier = string_in;
      } else if (type == "literal") {
        this->literal = string_in;
      } else {
        throw std::string("Bad 2nd parameter to Rval ctor");
      }
    }
  };
  class Arg
  {
  public:
    std::string name;
    Rval rval;
    Arg(const std::string& name_in, const Rval& rval_in) 
      : name(name_in), rval(rval_in) {}
    Arg() {}
  }; 
  class ParsedCommand
  { 
  public:
    ParsedCommand() : shellCommand(false) {}
    ParsedCommand(bool shell_command) : shellCommand(shell_command) {}
    bool isShellCommand() const { return shellCommand; }
    bool shellCommand;
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
  void shellCommand();


  static SurfpackParser& instance();
  yyFlexLexer& globalLexer();
  int yyparse(std::istream* is = &std::cin, std::ostream* os = &std::cout);
  std::vector<ParsedCommand>& commandList();

  // commands for use by interpreter to parse out certain argument types
  static std::string parseOutIdentifier(const std::string& argname,
    const ArgList& arglist);
  static std::string parseOutStringLiteral(const std::string& argname,
    const ArgList& arglist);
  static int parseOutInteger(const std::string& argname,
    const ArgList& arglist, bool& valid);
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




