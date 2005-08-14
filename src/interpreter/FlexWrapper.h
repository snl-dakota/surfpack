#ifndef FLEX_WRAPPER_H
#define FLEX_WRAPPER_H

class yyFlexLexer;

class FlexWrapper
{
public:
  FlexWrapper();
  ~FlexWrapper();
  void setParseStreams(std::istream& in, std::ostream& out);
  int nextToken();
  const char* currentToken();
private:
  yyFlexLexer* flexPtr;
};

#endif
