#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#if defined(__APPLE__) && defined(__MACH__) // Mac
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/glu.h>
#elif defined(WIN32) // windows
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#else // linux
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "delfem2/mshio.h"
#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"
#include "delfem2/funcs_gl.h"
#include "delfem2/dyntri_v3.h"

#include "delfem2/mshtopoio_gl.h"

CTriangulationOutput Triangulation
(const std::vector<double>& aXY,
 double edge_length)
{
  CTriangulationOutput out;
  std::vector< std::vector<double> > aaXY;
  aaXY.push_back(aXY);
  GenerateTesselation2(out.me.aElem, out.me.aPos,
                       out.aPtrVtxInd, out.aVtxInd,
                       edge_length, true, aaXY);
  out.me.ndim = 2;
  out.me.elem_type = MESHELEM_TRI;
  return out;
}

