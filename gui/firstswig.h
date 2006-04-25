/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2001, Sandia National Laboratories.

    Surfpack: A Software Library of Multidimensional Surface Fitting Methods

    Surfpack is distributed under the DAKOTA GNU General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

class FirstClass
{
public:
  FirstClass(int x_);
  ~FirstClass();
  void printVal(char* msg);
  double shiftArray(double* vals, int size);
  double myevaluate(double* vals, int size);
private:
  int x;
};
