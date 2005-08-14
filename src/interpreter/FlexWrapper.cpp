#include <FlexLexer.h>
#include <iostream>
#include <iostream.h>

#include "FlexWrapper.h"
//#define FLEX_USES_STD

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
#ifdef FLEX_USES_STD
  flexPtr->switch_streams(&in,&out);
#else
  flexPtr->switch_streams(
    reinterpret_cast< ::istream* >(&in),reinterpret_cast< ::ostream* >(&out)
  );
#endif
}

int FlexWrapper::nextToken()
{
  return flexPtr->yylex();
}

const char* FlexWrapper::currentToken()
{
  return flexPtr->YYText();
}
