/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __DIRECT_ANN_MODEL_H__
#define __DIRECT_ANN_MODEL_H__

#include "surfpack_system_headers.h"
#include "SurfpackModel.h"

class DirectANNBasisSet
{
public:
  MtxDbl weights;
  /// default constructor used when reading from archive file
  DirectANNBasisSet() { /* empty ctor */ }
  DirectANNBasisSet(const MtxDbl& weights_in);
  double eval(unsigned index, const VecDbl& x) const;
  double deriv(unsigned index, const VecDbl& x, const VecUns& vars) const;
  double nodeSum(unsigned index, const VecDbl& x) const;
  std::string asString() const;

private:


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private serialize function
  friend class boost::serialization::access;
  /// serializer for derived class data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

};

class DirectANNModel : public SurfpackModel
{

public:

  DirectANNModel(const DirectANNBasisSet& bs_in, const VecDbl& coeffs_in);
  virtual VecDbl gradient(const VecDbl& x) const;
  virtual std::string asString() const;

protected:

  /// default constructor used when reading from archive file
  DirectANNModel() { /* empty ctor */ }

  virtual double evaluate(const VecDbl& x) const;
  DirectANNBasisSet bs;
  VecDbl coeffs;

private:

  /// disallow copy construction as not implemented
  DirectANNModel(const DirectANNModel& other);

  /// disallow assignment as not implemented
  DirectANNModel& operator=(const DirectANNModel& other);

friend class DirectANNModelTest;


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
    // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

};

///////////////////////////////////////////////////////////
///   Direct ANN Model Factory	
///////////////////////////////////////////////////////////

class DirectANNModelFactory : public SurfpackModelFactory 
{

public:
  DirectANNModelFactory();
  DirectANNModelFactory(const ParamMap& args);
  MtxDbl randomMatrix(unsigned nrows, unsigned ncols);

protected:

  /// Model-specific portion of creation process
  virtual SurfpackModel* Create(const SurfData& sd);

  /// set member data prior to build; appeals to SurfpackModel::config()
  virtual void config();

  unsigned nodes;
  double range;
  unsigned samples;

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class Model data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif

};


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION

template< class Archive >
void DirectANNBasisSet::serialize(Archive & archive, 
				  const unsigned int version)
{
  archive & weights;
}

template< class Archive >
void DirectANNModel::serialize(Archive & archive,
			       const unsigned int version)
{
  archive & boost::serialization::base_object<SurfpackModel>(*this);
  archive & bs;
  archive & coeffs;
}

#endif

#endif
