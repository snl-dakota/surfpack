/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

/* lexical analyzer for Surfpack */
%{
	#include "surfpack_config.h"
	#include "surfparse.h"
	extern void appendToken(const char*);
%}
%%
CreateAxes |
CreateSurface |
CreateSample |
Evaluate |
Fitness |
Load |
Save               			{ appendToken(yytext); return STANDARD_COMMAND; }
[[:alpha:]][[:alnum:]_]*		{ appendToken(yytext); return IDENTIFIER; }
'[^']*'					{ appendToken(yytext); return STRING; }
#.*\n					{  }
!.*\n					{ appendToken(yytext); return SHELL_COMMAND; } 
\/\*(.*\n)*.*\*\/			{ appendToken(yytext); return SHELL_COMMAND; } 
[+-]?[0-9]+				{ appendToken(yytext); return INTEGER; }
[+-]?([0-9])*"."([0-9])*([eEdD][+-]?[0-9]{1,3})?	{ appendToken(yytext); return REAL; }
"("               			{ appendToken(yytext); return '('; }
")"               			{ appendToken(yytext); return ')'; }
"{"               			{ appendToken(yytext); return '{'; }
"}"               			{ appendToken(yytext); return '}'; }
","               			{ appendToken(yytext); return ','; }
"="               			{ appendToken(yytext); return '='; }
"["               			{ appendToken(yytext); return '['; }
"]"               			{ appendToken(yytext); return ']'; }
[ \t\r]				        { appendToken(yytext); }
"\n"
.               			{ appendToken(yytext); }
