// ----------------------------------------------------------
// Project: SURFPACK++
//
// File:        SurfData.cpp
// Author:      Eric Cyr
// Created:     June 6, 2002
// Modified:    June 17, 2002
//              June 24 - 26, 2002
//              July 11, 2002
//              May 14, 2003 Mark Richards
//              August 13, 2003 Mark Richards
//              August 14, 2003 Mark Richards
//
// Description: 
// + SurfData class - this is the container for all the data 
//   from which a surrogate/model is created.
// + Left shift (<<) operator for SurfData. Essentially so I
//   I can print out informaion
// ----------------------------------------------------------

#include "SurfData.h"
#include "Surface.h"

using namespace std;

#ifdef __TESTING_MODE__
   int SurfData::constructCount = 0;
   int SurfData::destructCount = 0;
#endif

// Constructs a SurfData object
// 
SurfData::SurfData() : dimension(0), responses(0) 
{
#ifdef __TESTING_MODE__
   constructCount++;
#endif
}

// copy constructor
//
SurfData::SurfData(const SurfData & data) : dimension(data.dimension), responses(data.responses)
{

#ifdef __TESTING_MODE__
   constructCount++;
#endif

    points = data.points;
}

// cleans up after a SurfData object
// 
SurfData::~SurfData() 
{
   // clean up the related surfaces
   listeners.clear();

   // clean up all the created points
   deletePoints();

#ifdef __TESTING_MODE__
   destructCount++;
#endif
}
/// retrieve a point from the data set
SurfPoint* SurfData::getPoint(int x) 
{
   if (x < 0 || static_cast<unsigned>(x) >= points.size()) {
       cerr << "Out of range in SurfData::getPoint" << endl;
       return 0;
   }
   return points[x];
}


int SurfData::addResponse()
{
    int newindex;
    if (points.size() == 0) {
	cerr << "Cannot add another response because there are no points." << endl;
	return -1;
    } else {
	newindex = points[0]->addResponse(0);
	responses = newindex + 1;
        for (unsigned i = 1; i < points.size(); i++) {
	    if (points[i]->addResponse() != newindex) {
		cerr << "Size mismatch among data points in this SurfData object.\nThis will likely cause a fatal error." << endl;
	    }
	}
    }
    return newindex;
}

// add a point to the data set
//
// point: a refernce to the point to be added
//        This point will be copied 
//
void SurfData::addPoint(const SurfPoint & pt) 
{
   // create a new object and add it 
   // it will get deleted when 'this' is deleted
   shallowAddPoint(new SurfPoint(pt));
}

// add a point to the data set
//
// point: a vector to the point to be added
// fval: the function value at this point
//
//void SurfData::addPoint(const vector<double> & pt,double fval) 
//{
//   // create a new object and add it (remember to
//   // delete it!!!)
//   shallowAddPoint(new SurfPoint(pt,fval));
//   notifyListeners();
//}

//void SurfData::addPoint(const vector<double> & pt,double fval,const vector<double> & grad) 
//{
//   // create a new object and add it (remember to
//   // delete it!!!)
//   shallowAddPoint(new SurfPoint(pt,fval,grad));
//   notifyListeners();
//}

// add a point to the data set
//
// point: a vector to the point to be added
// fvals: the response values at this point
//
//void SurfData::addPoint(const vector<double> & point,const vector<double> & fvals)
//{
//   shallowAddPoint(new SurfPoint(point,fvals));
//   notifyListeners();
//}

// add a point to the data set
//
// point: a vector to the point to be added
// fvals: the response values at this point
//
//void SurfData::addPoint(const std::vector<double> & point,
//                      const std::vector<double> & fvals,
//                      const std::vector<std::vector<double> > & grads)
//{
//   shallowAddPoint(new SurfPoint(point,fvals,grads));
//   notifyListeners();
//}

// add a point to the data set
//
// point: an array of "dim" values representing the domain
//        of the point to be added
// fvals: an array of "respCnt" values representing the response
//        values at this point
// grads: an array of "dim * respCnt" values representing the
//        gradient values at this point
// dim: dimension of point (the dimension of the domain)
// respCnt: the number of responses in this point
//
//void SurfData::addPoint(const double * point,const double * fvals,const double * grads,
//                        int dim,int respCnt)
//{
//   shallowAddPoint(new SurfPoint(point,fvals,grads,dim,respCnt));
//   notifyListeners();
//}

// add a set of points to the data
// 
// points: the vector of points to be added 
// 
//void SurfData::addPoints(const vector<SurfPoint*> & pts) 
//{
//   // create a copy of all the points in the other SurfData object
//   // remember to delete them!!!
//   for(unsigned int i=0;i<pts.size();i++)
//      shallowAddPoint(new SurfPoint(*pts[i])); 
//   notifyListeners();
//}

// the "FORTRAN" style interface to the SurfData object
//
// domain: the x values
// response: the various response values
//
//void SurfData::addPoints(const std::vector<std::vector<double> > & domain,
//                       const std::vector<std::vector<double> > & response)
//   throw(SurfException)
//{
//   bool criteria = domain.size() == response.size() && domain.size();
//   if(!criteria) 
//      throw SurfException("number of points to be created is not consistent"); 
//
//   for(unsigned int i=0;i<domain.size();i++) // for each point
//   {
//      SurfPoint * point = new SurfPoint(domain[i],response[i]);
//
//      shallowAddPoint(point);
//   }
//
//   notifyListeners();
//}

// add a set of points in a "FORTRAN" like format
// domain: the X (Vector) values for this surface
// response: the Y (Vector) values each one representing
//           a different surface
// grad: the gradient (Matrix) values each one repesenting
//       a different surface
//
//void SurfData::addPoints(const vector<vector<double> > & domain,
//                         const vector<vector<double> > & response, 
//                         const vector<vector<vector<double> > > & grad)
//   throw(SurfException)
//{
//   bool criteria = (domain.size() == response.size()) && (domain.size() == grad.size());
// 
//   if(!criteria) 
//      throw SurfException("number of points to be created is not consistent"); 
//
//   for(unsigned int i=0;i<domain.size();i++) // for each point
//   {
//      SurfPoint * point = new SurfPoint(domain[i],response[i],grad[i]);
//
//      shallowAddPoint(point);
//   }
//
//   notifyListeners();
//}

// add a point to the data set
//
// point: an array of "pointCnt * dim" values representing the domain
//        of the point to be added
// fvals: an array of "pointCnt * respCnt" values representing the response
//        values at this point
// grads: an array of "pointCnt * dim * respCnt" values representing the
//        gradient values at this point
// dim: dimension of point (the dimension of the domain)
// respCnt: the number of responses in this point
// pointCnt: the number of points to be added
//
//
//void SurfData::addPoints(const double * domain,const double * fvals,const double * grads,
//                        int pointCnt,int dim,int respCnt)
//{
//   for(int i=0;i<pointCnt;i++) 
//   {
//      SurfPoint * point = new SurfPoint(domain + i * dim, 
//                                        fvals  + i * respCnt, 
//                                        grads!=0?grads  + i * dim * respCnt:0, 
//                                        dim,respCnt);
//
//      shallowAddPoint(point);
//   }
//
//   notifyListeners();
//}
   
// add a set of data to this object
//
// data: SurfData object to have its data
//       unioned with this object 
//
//void SurfData::appendData(const SurfData & data) 
//{
//   // create a copy of all the points in the other SurfData object
//   // remember to delete them!!!
//   for(unsigned int i=0;i<data.points.size();i++)
//      shallowAddPoint(new SurfPoint(*data.points[i])); 
//   notifyListeners();
//}
//
// get the number of data points
//
//int SurfData::getPointCount() const { return points.size(); }

// get the number of data points
int SurfData::size() const { return points.size(); }

// get the smallest dimension of the points
int SurfData::getDimension() const { return dimension; }

// get the smallest # of responses of the points
int SurfData::getResponseCount() const { return responses; }

/// return a reference to the vector of points
vector<SurfPoint*> & SurfData::getPoints() { return points; }

//void SurfData::getIterator(AbstractDataIterator & iterator) const 
//{ getIterator(&iterator); }
//
//void SurfData::getIterator(AbstractDataIterator * iterator) const 
//{ 
//   if(iterator)
//   {
//      DataIterator local;
//      local.assign(points);
// 
//      iterator->assign(local);
//   }
//}

// so that you don't have to make the <<
// operator a friend 
//
// os: stream to be written to
//
//ostream & SurfData::print(ostream & os) const 
//{
//   os << "SurfData:" 
//      << " Point count = "  
//      << getPointCount()
//      << ", Dimension = " 
//      << getDimension();
//   return os;
//}

// read a set of data points from a file or other input stream
istream& SurfData::read(istream & is, bool binary)
{
  // the first lines in the files should specify how many xvals and fvals there are per line
  int xvals, fvals, gvals, npts;
  if (!binary) {
      is >> npts;
      is >> xvals;
      is >> fvals;
      is >> gvals;
  } else {
      is.read((char*)&npts,sizeof(int));
      is.read((char*)&xvals,sizeof(int));
      is.read((char*)&fvals,sizeof(int));
      is.read((char*)&gvals,sizeof(int));
  }

  for (int i = 0; i < npts; i++) {
	 shallowAddPoint(new SurfPoint(xvals, fvals, gvals, is, binary));
  } 
  return is;
}

ostream& SurfData::write(ostream & os, bool binary) const
{
   
   if (!binary) {
       os << points.size() << endl
          << dimension << endl 
          << responses << endl;
       os << 0 << endl; // grad size
   } else {
       unsigned s = points.size();
       int zero = 0;
       os.write((char*)&s,sizeof(int));
       os.write((char*)&dimension,sizeof(int));
       os.write((char*)&responses,sizeof(int));
       os.write((char*)&zero,sizeof(int));


   }
   vector<SurfPoint*>::const_iterator itr;
   itr = points.begin();
   while (itr != points.end()) {
	   (*itr)->write(os,binary);
	   if (!binary) {
	       os << endl; // Surfpoint->write(os) does not write newline after each point
  	   } 
	   ++itr;
   }
   return os;

}
// inform the SurfData object that the surface wants to be notified
// when something changes
//
void SurfData::addListener(Surface * surface)
{
   // only add the listener if its not already there
   list<Surface*>::iterator itr = find(listeners.begin(),listeners.end(),surface);
   if(itr==listeners.end())
      listeners.push_back(surface);
}

// remove surface from the list of surfaces to be notified
// when the data changes
//
void SurfData::removeListener(Surface * surface)
{
   // make sure its OK to erase the object, then do so
   list<Surface*>::iterator itr = find(listeners.begin(),listeners.end(),surface);
   if(itr!=listeners.end())
      listeners.erase(itr);
}


// overloaded operators
////////////////////////////////////
   
// assignment operator
//
SurfData & SurfData::operator=(const SurfData & sd)
{
   // clean up whats already here
   deletePoints();

   this->dimension = sd.dimension;
   this->responses = sd.responses;
   this->points = sd.points;
   // now all the surfaces need to be update
   // or at least told of a change
   notifyListeners();

   return (*this);
}

// so a SurfData object can be printed
ostream & operator<<(ostream & os,const SurfData & data) 
{ 
    return data.write(os); 
}

// so a SurfException object can be printed
//
ostream & operator<<(ostream & os,const SurfException & e) 
{ 
    return e.print(os); 
}

// make a "shallow" copy of the point. That
// is only copy the pointer to the point
// and not the entire data structure
// pt: pointer to the SurfPoint to be added
//
void SurfData::shallowAddPoint(SurfPoint * pt)
{
    if (dimension == 0) {
        dimension = pt->getDimension();
	responses = pt->getResponseCount();
    } else if (pt->getDimension() != dimension) {
	cerr << "Point cannot be added: dimension does not match dimension for this object's data points" << endl;
	return;
    } else if (pt->getResponseCount() != responses) {
	cerr << "Point cannot be added: number of responses does not match number of responses in this object's data points" << endl;
	return;
    }
    // make a "shallow" copy of this pointer
    points.push_back(pt);
    notifyListeners();
}

// private methods
////////////////////////////////////

// when ever a point is added the surfaces need 
// to be notified of a change, so they can rebuild 
//
void SurfData::notifyListeners()
{
    list<Surface*>::iterator itr;
    for(itr=listeners.begin();itr!=listeners.end();++itr)
       (*itr)->notify();
  
}

// clean up the points stored in this object
//
void SurfData::deletePoints()
{
   vector<SurfPoint*>::iterator itr;
   for(itr=points.begin();itr!=points.end();++itr)
      delete *itr;
   points.clear();
}

//void SurfData::removePoints()
//{
//    points.clear();
//}
