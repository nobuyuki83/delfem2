// exported.h
#pragma once

// Define EXPORTED for any platform
#ifdef _WIN32
# ifdef WIN_EXPORT
#   define EXPORTED  __declspec( dllexport )
# else
#   define EXPORTED  __declspec( dllimport )
# endif
#else
# define EXPORTED
#endif



// main.cpp
extern "C"
{
  EXPORTED int add(int a, int b)
  {
      return a + b;
  }
}
