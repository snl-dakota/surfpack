#ifndef SURFPACK_PARSER_ARGS_H
#define SURFPACK_PARSER_ARGS_H
#include "surfpack_config.h"
#include "surfpack_system_headers.h"

class Arg;
typedef std::vector<Arg> ArgList;
typedef std::vector<double> Tuple;
class Triplet 
{
public:
  Triplet();
  double min;
  double max;
  unsigned numPts;
};

class Rval
{
public:
  virtual int getInteger() const;
  virtual double getReal() const;
  virtual const Tuple& getTuple() const;
  virtual const Triplet& getTriplet() const;
  virtual const std::string& getIdentifier() const;
  virtual const std::string& getStringLiteral() const;
  virtual const ArgList& getArgList() const;
  virtual Rval* clone() const = 0;
  void noSuchValue() const;
  virtual ~Rval();
};

class RvalInteger: public Rval
{
public:
  RvalInteger(int value_in); 
  virtual Rval* clone() const;
  virtual int getInteger() const;
private:
  int value;
};

class RvalReal : public Rval
{
public:
  RvalReal(double value_in); 
  virtual Rval* clone() const;
  virtual double getReal() const;
private:
  double value;
};

class RvalTuple : public Rval
{
public:
  RvalTuple(const Tuple& value_in); 
  virtual Rval* clone() const;
  virtual const Tuple& getTuple() const;
private:
  Tuple value;
};

class RvalTriplet : public Rval
{
public:
  RvalTriplet(const Triplet& value_in); 
  virtual Rval* clone() const;
  virtual const Triplet& getTriplet() const;
private:
  Triplet value;
};

class RvalIdentifier : public Rval
{
public:
  RvalIdentifier(const std::string& value_in); 
  virtual Rval* clone() const;
  virtual const std::string& getIdentifier() const;
private:
  std::string value;
};

class RvalStringLiteral : public Rval
{
public:
  RvalStringLiteral(const std::string& value_in); 
  virtual Rval* clone() const;
  virtual const std::string& getStringLiteral() const;
private:
  std::string value;
};

class RvalArgList : public Rval
{
public:
  RvalArgList(const ArgList& value_in);
  virtual Rval* clone() const;
  virtual const ArgList& getArgList() const;
private:
  ArgList value;
};

class Arg
{
public:
  std::string name;
  const Rval* getRVal() const;
  Arg(const Arg& other);
  const Arg& operator=(const Arg& other);
  Arg(const std::string& name_in, Rval* rval_in);
  ~Arg();
  void setName(const std::string& name_in);
  void setRVal(Rval* rval_in);
  Arg(); 
private:
  Rval* rval;
}; 

#endif
