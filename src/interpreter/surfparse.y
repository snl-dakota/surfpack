/* parser for surfpack */

%{
#define YYINITDEPTH 400
#define YYSTYPE double
int yylex();
int yyerror(const char* s);

#include <iostream>
#include "SurfpackParser.h"
#define parser SurfpackParser::instance()
using namespace std;

%}

%token STANDARD_COMMAND
%token IDENTIFIER
%token STRING
%token INTEGER
%token REAL

%%

input:		/* empty */
	|	statement input
;

statement: 	STANDARD_COMMAND {parser.addCommandName();} arglist	
;

arglist:	'[' arg arglisttail ']'				
;

arglisttail: 	/* empty */
	| 	',' arg arglisttail
;

arg:		IDENTIFIER {parser.addArgName();} '=' argval
;

argval:		IDENTIFIER 					{parser.addArgValIdent()}
	|	INTEGER						{parser.addArgValInt()}
	|	STRING						{parser.addArgValString()}
	|	REAL						{parser.addArgValReal()}
	|	tuple						
	|	triplet						
	| 	arglist						{parser.addArgValArgList()}
;

tuple:		'(' number {parser.addTupleVal();} dimlisttail ')'
;

number:		INTEGER
	|	REAL
;

dimlisttail:	/* empty */	
	|	',' number {parser.addTupleVal();} dimlisttail
;

triplet:	'{' number {parser.addTripletMin();} ',' 
		number {parser.addTripletMax();} ',' 
		INTEGER {parser.addTripletNumPts();} '}'
;

	
%%

int yyerror (const char* s) /* Called by yyparse on error */
{
  cout << s << endl;
}

int yylex()
{
  return parser.globalLexer().yylex();
}
