#include "FlexWrapper.h"
#include <FlexLexer.h>


FlexWrapper::FlexWrapper()
{
  flexPtr = new yyFlexLexer;
}

FlexWrapper::~FlexWrapper()
{
  delete flexPtr;
}

void FlexWrapper::setParseStreams(std::istream& in, std::ostream& out)
{
  flexPtr->switch_streams(&in,&out);
}

int FlexWrapper::nextToken()
{
  flexPtr->yylex();
}

const char* FlexWrapper::currentToken()
{
  return flexPtr->YYText();
}
