/* lexical analyzer for Surfpack */
%{
	#include <iostream>
	#include "surfparse.tab.h"
	#include <sstream>
	using namespace std;
        ostringstream cmdstream;
%}
%%
LoadSurface |
LoadData |
CreateSurface |
ConvertData |
ConvertSurface |
MonteCarloSample |
GridPoints |
Evaluate |
SaveData |
SaveSurface |
PointDefinition |
Fitness               			{ /*cmdstream << yytext; cout << "Keyword: " << yytext << endl;*/ return STANDARD_COMMAND; }
[[:alpha:]][[:alnum:]_]*		{ /*cmdstream << yytext; cout << "Identifier: " << yytext << endl;*/ return IDENTIFIER;}
'[^']*'					{ /*cmdstream << yytext; cout << "StringLiteral: " << yytext << endl;*/ return STRING; }
#.*\n					{ /*cmdstream << yytext; cout << "Comment: " << yytext << endl;*/ }
!.*\n					{ /*cmdstream << yytext; cout << "ShellCommand: " << yytext << endl;*/ } 
[+-]?[0-9]+					{ /*cmdstream << yytext; cout << "Integer: " << yytext << endl;*/ return INTEGER;}
[+-]?([0-9])*"."([0-9])*([eEdD][+-]?[0-9]{1,3})?	{ /*cmdstream << yytext; cout << "Real: " << yytext << endl;*/ return REAL;}
"("               			{ /*cmdstream << yytext; cout << "Punctuation: " << yytext << endl;*/ return '('; }
")"               			{ /*cmdstream << yytext; cout << "Punctuation: " << yytext << endl;*/ return ')'; }
"{"               			{ /*cmdstream << yytext; cout << "Punctuation: " << yytext << endl;*/ return '{'; }
"}"               			{ /*cmdstream << yytext; cout << "Punctuation: " << yytext << endl;*/ return '}'; }
","               			{ /*cmdstream << yytext; cout << "Punctuation: " << yytext << endl;*/ return ','; }
"="               			{ /*cmdstream << yytext; cout << "Punctuation: " << yytext << endl;*/ return '='; }
"["               			{ /*cmdstream << yytext; cout << "LeftBracket: " << yytext << endl;*/ return '['; }
"]"               			{ /*cmdstream << yytext; cout << "RightBracket: " << yytext << endl;*/ return ']'; }
[ \t\r]				        { cmdstream << yytext;}
"\n"
.               			{ /*cmdstream << yytext; cout << "Unrecognized: " << yytext << endl;*/ }
