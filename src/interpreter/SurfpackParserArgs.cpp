#include "surfpack_config.h"
#include "SurfpackParserArgs.h"

Triplet::Triplet() : min(0), max(0), numPts(0) {}

int Rval::getInteger() const
{
  noSuchValue();
}

double Rval::getReal() const
{
  noSuchValue();
}

const Tuple& Rval::getTuple() const
{
  noSuchValue();
}

const Triplet& Rval::getTriplet() const
{
  noSuchValue();
}

const std::string& Rval::getIdentifier() const
{
  noSuchValue();
}

const std::string& Rval::getStringLiteral() const
{
  noSuchValue();
}

void Rval::noSuchValue() const
{
  throw std::string("This Rval class does not have such a value");
}

RvalInteger::RvalInteger(int value_in) : value(value_in) 
{

}

int RvalInteger::getInteger() const
{
  return value;
}

Rval* RvalInteger::clone() const
{
  return new RvalInteger(value);
}

RvalReal::RvalReal(double value_in) : value(value_in) 
{

}

double RvalReal::getReal() const
{
  return value;
}

Rval* RvalReal::clone() const
{
  return new RvalReal(value);
}

RvalTuple::RvalTuple(const Tuple& value_in) : value(value_in) 
{

}

const Tuple& RvalTuple::getTuple() const
{
  return value;
}

Rval* RvalTuple::clone() const
{
  return new RvalTuple(value);
}

RvalTriplet::RvalTriplet(const Triplet& value_in) : value(value_in) 
{

}

const Triplet& RvalTriplet::getTriplet() const
{
  return value;
}

Rval* RvalTriplet::clone() const
{
  return new RvalTriplet(value);
}

RvalIdentifier::RvalIdentifier(const std::string& value_in) : value(value_in) 
{

}

const std::string& RvalIdentifier::getIdentifier() const
{
  return value;
}

Rval* RvalIdentifier::clone() const
{
  return new RvalIdentifier(value);
}

RvalStringLiteral::RvalStringLiteral(const std::string& value_in) 
  : value(value_in) 
{

}

const std::string& RvalStringLiteral::getStringLiteral() const
{
  return value;
}

Rval* RvalStringLiteral::clone() const
{
  return new RvalStringLiteral(value);
}

Arg::Arg(const std::string& name_in, Rval* rval_in)
  : name(name_in), rval(rval_in)
{

}

Arg::Arg()
 : rval(0)
{

}

void Arg::setName(const std::string& name_in)
{
  name = name_in;
}

const Rval* Arg::getRVal() const
{
  return rval;
}

void Arg::setRVal(Rval* rval_in)
{
  delete rval;
  rval = rval_in;
}

Arg::~Arg()
{
  delete rval;
  rval = 0;
}

