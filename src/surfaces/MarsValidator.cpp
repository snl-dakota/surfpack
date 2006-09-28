/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#endif
#include "surfpack.h"
#include "SurfData.h"
#include "MarsValidator.h"

using std::cerr;
using std::cout;
using std::endl;
using std::istream;
using std::numeric_limits;
using std::ostream;
using std::ostream_iterator;
using std::ostringstream;
using std::set;
using std::setw;
using std::string;
using std::vector;

void printMatrixCpp(double* mat, unsigned rows, unsigned columns, ostream& os)
{
  for (unsigned r = 0; r < rows; r++) {
    for (unsigned c = 0; c < columns; c++) {
      os << setw(15) << mat[r + c*rows];
    }
    os << endl;
  }
}

void printIntMatrixCpp(int* mat, unsigned rows, unsigned columns, ostream& os)
{
  for (unsigned r = 0; r < rows; r++) {
    for (unsigned c = 0; c < columns; c++) {
      os << setw(15) << mat[r + c*rows];
    }
    os << endl;
  }
}

//_____________________________________________________________________________
// Mars Basis Function class 
//_____________________________________________________________________________
string MarsBasis::asString() const
{
  ostringstream os;
  os << ((sign < 0) ? "-" : " ") << "(v" << var_index << " - " 
     << knot << ")";
  if (multiplier) os << " * " << multiplier->asString();
  return os.str();
}
MarsBasis::MarsBasis(double sign_in, double knot_in, unsigned var_index_in,
  const Basis* multiplier_in)
  : Basis(), sign(sign_in), knot(knot_in), var_index(var_index_in), 
  multiplier(multiplier_in)
{

}

Basis* MarsBasis::clone() const
{
  return new MarsBasis(sign,knot,var_index,multiplier);
}

double MarsBasis::eval(const vector<double>& x) const
{
  assert(var_index < x.size());
  double val = sign*(x[var_index] - knot);
  if (val <= 0.0) {
    return 0.0;
  } else {
    assert(this != multiplier);
    return multiplier ? multiplier->eval(x)*val : val;
  } 
}

UnityBasis::UnityBasis() : Basis() {}

string UnityBasis::asString() const
{
  return "1.0";
}

double UnityBasis::eval(const vector<double>& x) const
{
  return 1.0;
}

Basis* UnityBasis::clone() const
{
  return new UnityBasis;
}
//_____________________________________________________________________________
// Data members 
//_____________________________________________________________________________

const string MarsCppSurface::name = "marsc";

//_____________________________________________________________________________
// Creation, Destruction, Initialization
//_____________________________________________________________________________

MarsCppSurface::MarsCppSurface(SurfData* sd) : Surface(sd)
{
  init();
}

MarsCppSurface::MarsCppSurface(const string filename) : Surface(0)
{
  init();
  read(filename);
}
MarsCppSurface::~MarsCppSurface()
{
}

void MarsCppSurface::init()
{
  max_bases = 3;
  max_interactions = 2;
  interpolation = 1;
}
//_____________________________________________________________________________
// Overloaded Operators 
//_____________________________________________________________________________

//_____________________________________________________________________________
// Queries
//_____________________________________________________________________________

const string MarsCppSurface::surfaceName() const
{
  return name;
}

unsigned MarsCppSurface::minPointsRequired() const
{
  if (xsize <= 0) {
    throw string(
      "Dimensionality of data needed to determine number of required samples."
    );
  } else {
    return xsize+1; 
  }
}

double MarsCppSurface::evaluate(const vector<double>& x)
{
  assert (coeffs.size() == bfs.size());
  double sum = 0;
  for (unsigned i = 0; i < bfs.size(); i++) {
    sum += coeffs[i]*bfs[i]->eval(x);
  } 
  return sum;
}

//_____________________________________________________________________________
// Commands 
//_____________________________________________________________________________

double lsv (SurfpackMatrix< double >& x_matrix, vector<double>& coeffs, 
  const vector<double>& y_vector)
{
  vector<double> result;
  surfpack::matrixVectorMult(result,x_matrix,coeffs);
  double sum = 0;
  double diff;
  for (unsigned i = 0;i < y_vector.size(); i++) {
    diff = y_vector[i] - result[i];
    sum += diff*diff;
  }
  return sum; 
}

void MarsCppSurface::build(SurfData& data)
{
  vector< set< double > > knot_candidates(data.xSize());
  for (unsigned i = 0; i < data.size(); i++) {
    for (unsigned j = 0; j < data.xSize(); j++) {
      knot_candidates[j].insert(data[i].X()[j]);
    }
  }
  for (unsigned k = 0; k < knot_candidates.size(); k++) {
    cout << "knot candidates for dimensions " << k << endl;
    copy(knot_candidates[k].begin(),knot_candidates[k].end(),
	ostream_iterator<double>(cout,"\n"));
  }
  
  SurfpackMatrix< double > x_matrix;
  SurfpackMatrix< double > x_matrix_copy;
  vector< double > y_vector(data.size());
  vector< double > y_vector_copy;
  for (unsigned i = 0; i < data.size(); i++) y_vector[i] = data[i].F();
  copy(y_vector.begin(),y_vector.end(),ostream_iterator<double>(cout,"\n"));
  bfs.clear();
  bfs.push_back(new UnityBasis);

  while (bfs.size() < 10 ) {
    unsigned bfsize = bfs.size();
    // Initialize the best values to garbage
    double bestlsv = numeric_limits<double>::max();
    Basis* best1 = 0;
    Basis* best2 = 0;
    // Iterate through each combination of parent basis, variable, and knot 
    for (unsigned parenti = 0; parenti < bfsize; parenti++) {
      Basis* parent = bfs[parenti];
      cout << "Parent basis (" << parenti << "/" << (bfsize-1) << "): " 
	   << parent->asString() << endl;
      for (unsigned var = 0; var < data.xSize(); var++) {
        cout << "   var: " << var << endl;
        for(set< double >::iterator itr = knot_candidates[var].begin();
	    itr != knot_candidates[var].end(); itr++) {
          double knot = *itr;
          cout << "        knot: " << knot << endl;
          // delete candidates from previous iteration, if any
          for(unsigned i = bfsize; i < bfs.size(); i++) {
            delete bfs[i];
            bfs[i] = 0;
          }
          bfs.erase(bfs.begin()+bfsize,bfs.end());
          assert(!bfs.empty());

          // create the new candidates
          bfs.push_back(new MarsBasis(1.0,knot,var,parent));
          bfs.push_back(new MarsBasis(-1.0,knot,var,parent));

	  // compute the least squares values
          // This can go in its own function later
          x_matrix.reshape(data.size(),bfs.size());
          for (unsigned r = 0; r < data.size(); r++) {
            for (unsigned c = 0; c < bfs.size(); c++) {
              x_matrix[r][c] = bfs[c]->eval(data[r].X());
            }
          }
          y_vector_copy = y_vector;
          x_matrix_copy = x_matrix;
          surfpack::linearSystemLeastSquares(x_matrix_copy,this->coeffs,
            y_vector_copy); 
          double newlsv = lsv(x_matrix,coeffs,y_vector);
          cout << "       newlsv: " << newlsv << endl;
          if (newlsv < bestlsv) {
            bestlsv = newlsv;
            delete best1;
            delete best2;
            best1 = new MarsBasis(1.0,knot,var,parent);
            best2 = new MarsBasis(-1.0,knot,var,parent);
            cout << "       New bests " << endl;
            cout << "       " << best1->asString() << endl;
            cout << "       " << best2->asString() << endl;
          }
        } // for each knot value
      } // for each variable
    } // for each parent basis 
          // delete candidates from previous iteration, if any
          for(unsigned i = bfsize; i < bfs.size(); i++) {
            delete bfs[i];
            bfs[i] = 0;
          }
          bfs.erase(bfs.begin()+bfsize,bfs.end());
    bfs.push_back(best1->clone());
    bfs.push_back(best2->clone());
    delete best1;
    delete best2;
  } // while not too many basis functions
  

  // Now solve one last time with the right bases
  // compute the least squares values
  // This can go in its own function later
          x_matrix.reshape(data.size(),bfs.size());
          for (unsigned r = 0; r < data.size(); r++) {
            for (unsigned c = 0; c < bfs.size(); c++) {
              x_matrix[r][c] = bfs[c]->eval(data[r].X());
            }
          }
          y_vector_copy = y_vector;
          x_matrix_copy = x_matrix;
          surfpack::linearSystemLeastSquares(x_matrix_copy,this->coeffs,
            y_vector_copy); 
          double finallsv = lsv(x_matrix,coeffs,y_vector);
          cout << "       finallsv: " << finallsv << endl;
          cout << "       finallsv/n: " << finallsv/data.size() << endl;

}

void MarsCppSurface::config(const Arg& arg)
{
  string argname = arg.name;
  if (argname == "max_bases") {
    max_bases = arg.getRVal()->getInteger();
  } else if (argname == "max_interactions") {
    max_interactions = arg.getRVal()->getInteger();
  } else if (argname == "interpolation") {
    if (arg.getRVal()->getStringLiteral() == "linear") {
      interpolation = 1;
    } else if (arg.getRVal()->getStringLiteral() == "cubic") {
      interpolation = 2;
    } else {
      cerr << "Expected value for interpolation: 'linear' or 'cubic'" << endl;
    } 
  } else {
    Surface::config(arg);
  }
}
/// Create a surface of the same type as 'this.'  This objects data should
/// be replaced with the dataItr passed in, but all other attributes should
/// be the same (e.g., a second-order polynomial should return another 
/// second-order polynomial.  Surfaces returned by this method can be used
/// to compute the PRESS statistic.
MarsCppSurface* MarsCppSurface::makeSimilarWithNewData(SurfData* surfData)
{
  return new MarsCppSurface(surfData);
}

//_____________________________________________________________________________
// Helper methods 
//_____________________________________________________________________________

//_____________________________________________________________________________
// I/O 
//_____________________________________________________________________________

void MarsCppSurface::writeBinary(ostream& os)
{
}

void MarsCppSurface::writeText(ostream& os)
{
}

void MarsCppSurface::readBinary(istream& is)
{

}

void MarsCppSurface::readText(istream& is)
{
}

//_____________________________________________________________________________
// Testing 
//_____________________________________________________________________________

