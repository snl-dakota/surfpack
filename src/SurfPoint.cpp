// ----------------------------------------------------------------------
// SurfPoint (old SurrogateBasePoint) is the base for points used by a Surface.
// The point itself is stored as 'x'. We store the function value 'fvals'
// at the point
// ----------------------------------------------------------------------
//
// Project: SURFPACK++
//
// File:       SurfPoint.cpp
// Author:     ???
// Created:    ???
// Modified:   June 17, 2002
//             June 27, 2002
//             July 11, 2002
// 
// Description
// + SurfPoint class - a container class for a point, response value
//   and gradient information
// + Left shift (<<) operator for SurfPoint. Easy printing.
// ----------------------------------------------------------------------


#include "SurfPoint.h"

using namespace std;


#ifdef __TESTING_MODE__
   int SurfPoint::constructCount = 0;
   int SurfPoint::copyCount = 0;
   int SurfPoint::destructCount = 0;
#endif

//SurfPoint::SurfPoint() : x(0), fvals(1), grad(1), responseCount(0)
//{
//  // Default constructor.
//
//#ifdef __TESTING_MODE__
//   constructCount++;
//#endif
//
//  fvals[0] = 0;
//}
SurfPoint::SurfPoint(int nxvals, int nfvals, int ngvals, istream & is, bool binary) 
{

#ifdef __TESTING_MODE__
   constructCount++;
#endif

   x.resize(nxvals);
   fvals.resize(nfvals);
//   grad.resize(ngvals);
	   
   if (!is.eof()) {
	if (!binary) {
		for (int i = 0; i < nxvals; i++) {
		   is >> x[i];
		}
		for (int j = 0; j < nfvals; j++) {
		   is >> fvals[j];
		}
	//	for (int k = 0; k < ngvals; k++) {
	////	   is >> grad[k];
	//	}
	} else { // read the point in binary format
		for (int i = 0; i < nxvals; i++) {
		   is.read((char*)&x[i],sizeof(double));
		}
		for (int j = 0; j < nfvals; j++) {
		   is.read((char*)&fvals[j],sizeof(double));
		}
	//	for (int k = 0; k < ngvals; k++) {
	////	   is >> grad[k];
	//	}
	}
   }   
   responseCount = fvals.size();
}
SurfPoint::SurfPoint(const SurfPoint& pt) : responseCount(pt.responseCount)
{
  // Copy constructor. Copies 'x' and 'f' from 'pt'.

#ifdef __TESTING_MODE__
   constructCount++;
   copyCount++;
#endif

  x = pt.x;
  fvals = pt.fvals;
//  grad = pt.grad;
}

SurfPoint::SurfPoint(const vector<double>& x_in, double f_in) : fvals(1), responseCount(1)
{
  // Constructor with 'x' and 'f' specified.

#ifdef __TESTING_MODE__
   constructCount++;
#endif

  x = x_in;
  fvals[0] = f_in;
}

SurfPoint::SurfPoint(const vector<double>& x_in) : responseCount(0)
{
    x = x_in;
}
//SurfPoint::SurfPoint(const vector<double>& x_in, double f_in,const vector<double> & grad_in) 
//   : fvals(1), grad(1), responseCount(1)
//{
//#ifdef __TESTING_MODE__
//   constructCount++;
//#endif
//
//   x = x_in;
//   fvals[0] = f_in;
//   grad[0] = grad_in;
//}

SurfPoint::SurfPoint(const vector<double> & x_in,const vector<double> & f_in)
//   : grad(1) 
{
#ifdef __TESTING_MODE__
   constructCount++;
#endif

   x = x_in;
   fvals = f_in; 
 //  grad.clear(); 
 //  grad.resize(fvals.size());
   responseCount = fvals.size();
}

//SurfPoint::SurfPoint(const vector<double> & x_in, 
//                     const vector<double> & f_in, 
//                     const vector<vector<double> > & grad_in) 
//{
//#ifdef __TESTING_MODE__
//   constructCount++;
//#endif
//   
//   x = x_in;
//   fvals = f_in; 
//   grad = grad_in; 
//   responseCount = fvals.size();
//}

//SurfPoint::SurfPoint(const double * x_in,const double * f_in,const double * grad_in, 
//                     int dim, int respCnt) 
//   : x(dim), fvals(respCnt), grad(respCnt)
//{
//#ifdef __TESTING_MODE__
//   constructCount++;
//#endif
// 
//   // fill the x values
//   for(int i=0;i<dim;i++)
//      x[i] = x_in[i];
//
//   // now fill the grad and fvals values
//   for(int j=0;j<respCnt;j++)
//   {
//      if(grad_in!=0)
//      {
//         for(int i=0;i<dim;i++)
//            grad[j].push_back(grad_in[j*dim+i]); 
//      }
//      fvals[j] = f_in[j];
//   }
//
//   responseCount = fvals.size();
//}

SurfPoint::~SurfPoint() 
{
#ifdef __TESTING_MODE__
   destructCount++;
#endif
}

int SurfPoint::getDimension() const
{ return x.size(); }

int SurfPoint::getResponseCount() const
{ return responseCount; }

void SurfPoint::write(ostream & os, bool binary) const
{
  if (!binary) {
      int width = 15;
      //os.setf(ios::scientific);
      for (unsigned int i = 0; i < x.size(); i++) {
        os  << setw(width) << x[i] ;
      }
      for (unsigned int j = 0; j < fvals.size(); j++) {
        os << setw(width) << fvals[j];
      }
  } else {
      for (unsigned int i = 0; i < x.size(); i++) {
        os.write((char*)&x[i],sizeof(double)) ;
      }
      for (unsigned int j = 0; j < fvals.size(); j++) {
        os.write((char*)&fvals[j],sizeof(double));
      }
  }
      
  //os << endl;  // omit newline to allow for printing more stuff on the same line
  //os.unsetf(ios::scientific);
}
SurfPoint& SurfPoint::operator=(const SurfPoint& pt) 
{
  // Assignment operator. Copies 'x',grad and 'fvals' from 'pt'.

  if (this == (&pt))  
    return (*this);
  
  x = pt.x;
  fvals = pt.fvals;
//  grad = pt.grad;
  responseCount = pt.responseCount;
  
  return (*this);
}

//bool SurfPoint::operator==(const SurfPoint& pt) const
//{ 
////   bool xequal,fequal,gradequal;
////   xequal = equal(x.begin(),x.end(),pt.x.begin());
////   fequal = equal(fvals.begin(),fvals.end(),pt.fvals.begin());
////   gradequal = equal(grad.begin(),grad.end(),pt.grad.begin(),gradCmp);
////   return xequal && fequal && gradequal;
// 
//   return  x==pt.x && fvals==pt.fvals && grad==pt.grad;
//}
//
//bool SurfPoint::operator<(const SurfPoint& pt) const
//{ return fvals<pt.fvals && grad<pt.grad && fvals<pt.fvals; }
//
//bool SurfPoint::operator>(const SurfPoint& pt) const
//{
//  bool notLessThan = !(*this < pt);
//  bool notEqualTo  = !(*this == pt);
//
//  return (notLessThan && notEqualTo); 
//}
//
//void SurfPoint::setX(const vector<double>& xvec)
//{ x = xvec; }

const vector<double> & SurfPoint::getX() const
{ return x; }

void SurfPoint::setF(int response, double fval)
{ 
    if (response < 0 || response >= responseCount) {
        cerr << "Invalid response index (" << response << ") supplied; no update was made" << endl;
	return;
    }
    fvals[response] = fval; 
}

//void SurfPoint::setResponses(const std::vector<double> & f_in)
//{
//   fvals = f_in;
//
//   // build up the gradients, so that the next ones
//   // are in the right place
//   grad.clear();
//   grad.resize(f_in.size()); 
//
//   responseCount = fvals.size();
//}
//
//void SurfPoint::setResponses(const std::vector<double> & f_in,
//                            const std::vector<std::vector<double> > & grad_in)
//{
//
//   fvals = f_in;
//   grad = grad_in;
//
//   responseCount = fvals.size();
//}
//
//
//void SurfPoint::setResponses(const double * f_in,
//                             const double * grad_in,int respCnt)
//{
//   int dim = getDimension();
// 
//   // pre-allocate space, its a little quicker
//   fvals.resize(respCnt);
//   grad.resize(respCnt);
//
//   for(int i=0;i<respCnt;i++)
//   {
//      if(grad_in!=0)
//      {
//         for(int j=0;j<dim;j++)
//            grad[i].push_back(grad_in[i*dim+j]);
//      }
//      fvals[i] = f_in[i];
//   }
//
//   responseCount = fvals.size();
//}
//
//void SurfPoint::addResponse(double f_in, const std::vector<double> & grad_in)
//{
//   if(responseCount==0)  // this is already allocated
//   {
//      fvals[0] = f_in;
//      grad[0] = grad_in;
//   }
//   else                  // tack on the points
//   {
//      fvals.push_back(f_in);
//      grad.push_back(grad_in);
//   }
//
//   responseCount = fvals.size();
//}
//
//void SurfPoint::addResponse(double f_in,const double * grad_in)
//{
//   vector<double> moreGrad(0);
//
//   if(grad_in!=0)
//   {
//      for(int i=0;i<getDimension();i++)
//         moreGrad.push_back(grad_in[i]);
//   }
//
//   addResponse(f_in,moreGrad);
//
///*   if(responseCount==0)
// *  {
// *     fvals[0] = f_in;
// *     grad[0] = moreGrad;
// *  }
// *  else
// *  {
// *     fvals.push_back(f_in);
// *     grad.push_back(moreGrad);
// *  }
// *  
// *  responseCount = fvals.size();
// */
//}

int SurfPoint::addResponse(double val)
{
    fvals.resize(++responseCount);
    fvals[responseCount - 1] = val;
    return responseCount - 1;
}

double SurfPoint::getF(int i) const
{ return fvals[i]; }

//void SurfPoint::setGrad(const vector<double> & g,int index)
//{ grad[index] = g; }
//
//void SurfPoint::getGrad(vector<double> & g,int index) const
//{ g = grad[index]; }
//
//const std::vector<double> & SurfPoint::getGrad(int i) const
//{ return grad[i]; }

ostream& SurfPoint::print(ostream& stream) const
{
  // Output the contents of this object to the current stream with no
  // line feeds. Return a pointer to the stream. This is used in
  // conjunction with ostream to avoid making operator<< a friend of
  // this object.

  stream.setf(ios::scientific);

  stream << "f=" << fvals[0];
  stream << " x=[ ";
  for (unsigned int i = 0; i < x.size(); i++)
    stream  << x[i] << " ";
  stream << "]";

  stream.unsetf(ios::scientific);

  return stream;
}

ostream& operator<<(ostream& stream, const SurfPoint & pt) 
{
  return pt.print(stream);
}
