/* lexical analyzer for Surfpack */
%{
	#include <iostream>
	#include "surfparse.tab.h"
	using namespace std;
%}
%%
LoadSurface |
LoadData |
CreateSurface |
Convert |
Sample |
GridPoints |
Evaluate |
SaveData |
SaveSurface |
PointDefinition |
Fitness               			{ cout << "Keyword: " << yytext << endl; return STANDARD_COMMAND; }
[[:alpha:]][[:alnum:]_]*		{ cout << "Identifier: " << yytext << endl; return IDENTIFIER;}
'[^']*'					{ cout << "StringLiteral: " << yytext << endl; return STRING; }
#.*\n					{ cout << "Comment: " << yytext << endl; }
!.*\n					{ cout << "ShellCommand: " << yytext << endl; } 
[+-]?[0-9]+					{ cout << "Integer: " << yytext << endl; return INTEGER;}
[+-]?([0-9])*"."([0-9])*([eE][+-]?[0-9]{0,3})?	{ cout << "Real: " << yytext << endl; return REAL;}
"("               			{ cout << "Punctuation: " << yytext << endl; return '('; }
")"               			{ cout << "Punctuation: " << yytext << endl; return ')'; }
"{"               			{ cout << "Punctuation: " << yytext << endl; return '{'; }
"}"               			{ cout << "Punctuation: " << yytext << endl; return '}'; }
","               			{ cout << "Punctuation: " << yytext << endl; return ','; }
"="               			{ cout << "Punctuation: " << yytext << endl; return '='; }
"["               			{ cout << "LeftBracket: " << yytext << endl; return '['; }
"]"               			{ cout << "RightBracket: " << yytext << endl; return ']'; }
[ \t\r\n]
.               			{ cout << "Unrecognized: " << yytext << endl; }

