#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <set>
#include <stack>
#include "delfem2/mshmisc.h"
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vec3.h"
#include "delfem2/slice.h"
#include "delfem2/primitive.h"
#include "delfem2/color.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/gl2_funcs.h"
#include "delfem2/opengl/gl2_color.h"

// -------------------------

std::vector<double> aXY;
std::vector<unsigned int> aTri;
std::vector<double> aVal;
std::vector<delfem2::CSliceTriMesh::CSegInfo> aSeg;

// ---------------------------

void myGlutDisplay(void)
{
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  delfem2::opengl::DrawMeshTri2D_Edge(aXY.data(), aXY.size()/2, aTri.data(), aTri.size()/3);
  
  std::vector< std::pair<double,delfem2::CColor> > colorMap;
  delfem2::ColorMap_BlueCyanGreenYellowRed(colorMap,
                                           -1, 1);
  delfem2::opengl::DrawMeshTri2D_ScalarP1(aXY.data(), aXY.size()/2,
                                          aTri.data(), aTri.size()/3,
                                          aVal.data(), 1,
                                          colorMap);
  ::glLineWidth(5);
  ::glBegin(GL_LINES);
  ::glColor3d(0,0,0);
  for(int iseg=0;iseg<aSeg.size();++iseg){
    double pA[2], pB[2];
    aSeg[iseg].Pos2D(pA,pB,
                     aXY, aTri);
    ::glVertex2dv(pA);
    ::glVertex2dv(pB);
  }
  ::glEnd();
  
}

void Hoge()
{
  { // make mesh
    const int ndiv = 16;
    std::vector<unsigned int> aQuad;
    delfem2::MeshQuad2D_Grid(aXY, aQuad,
                             ndiv,ndiv);
    delfem2::convert2Tri_Quad(aTri,
                              aQuad);
    delfem2::Translate(aXY,
                       -ndiv*0.5, -ndiv*0.5);
    delfem2::Scale(aXY,
                   1.0/ndiv);
  }
  
  { // make value
    const unsigned int np = aXY.size()/2;
    aVal.resize(np);
    for(int ip=0;ip<np;++ip){
      double x0 = aXY[ip*2+0];
      double y0 = aXY[ip*2+1];
      aVal[ip] = sqrt(x0*x0+y0*y0)*4.0-1.5;
    }
  }
  

  
}

int main(int argc,char* argv[])
{
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode  = delfem2::CAMERA_ROT_TBALL;
  
  Hoge();
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    {
      static int iframe = 0;
      double thres = 0.9*sin(iframe*0.001);
      aSeg.clear();
      for(unsigned int itri=0;itri<aTri.size()/3;++itri){
        double v0 = aVal[ aTri[itri*3+0] ];
        double v1 = aVal[ aTri[itri*3+1] ];
        double v2 = aVal[ aTri[itri*3+2] ];
        if(   (v0-thres)*(v1-thres) >= 0
           && (v1-thres)*(v2-thres) >= 0
           && (v2-thres)*(v0-thres) >= 0 ){
          continue;
        }
        delfem2::CSliceTriMesh::CSegInfo seg;
        seg.Initialize(itri, aTri, aVal, thres);
        aSeg.push_back(seg);
      }
      iframe += 1;
    }
    // ------------------
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
