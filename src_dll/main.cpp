#include "exported.h"

// main.cpp
extern "C"
{
  EXPORTED int add(int a, int b)
  {
      return a + b;
  }
}
