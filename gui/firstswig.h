
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
