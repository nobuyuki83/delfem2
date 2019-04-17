#include <iostream>
#include <math.h>

#include "delfem2/mathexpeval.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc,char* argv[])
{
  CMathExpressionEvaluator e;
  e.SetKey("x", 3.0);
  e.SetExp("x+3.0");
  e.SetKey("x", 5.0);
  std::cout << e.Eval() << std::endl;
  return 0;
}
