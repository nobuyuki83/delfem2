#include <cstdlib>
#include <vector>
#include <set>
#include "delfem2/dtri.h"

#include "delfem2/objfunc_v23.h"
#include "delfem2/objfunc_v23dtri.h"
#include "delfem2/dtri_v2.h"

// --------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/glold_v23dtricad.h"
#include "delfem2/opengl/glold_funcs.h"

namespace dfm2 = delfem2;

// ------------------------------------

std::vector<dfm2::CDynPntSur> aPo2D;
std::vector<dfm2::CDynTri> aETri;
std::vector<dfm2::CVector2> aVec2;
std::vector<unsigned int> aLine;
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aXYZt;
std::vector<double> aUVW; // deformed vertex velocity
std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)
const double mass_point = 0.01;
const double dt = 0.01;
const double gravity[3] = {0.0, 0.0, -10.0};
bool is_animation = false;

// -------------------------------------

void StepTime()
{
  dfm2::PBD_Pre3D(aXYZt,
                  dt, gravity, aXYZ, aUVW, aBCFlag);
  dfm2::PBD_TriStrain(aXYZt.data(),
                aXYZt.size()/3, aETri, aVec2);
  dfm2::PBD_Bend(aXYZt.data(),
           aXYZt.size()/3, aETri, aVec2);
  dfm2::PBD_Seam(aXYZt.data(),
                 aXYZt.size()/3, aLine.data(), aLine.size()/2);
  dfm2::PBD_Post(aXYZ, aUVW,
                 dt, aXYZt, aBCFlag);

}

// -------------------------------------

void myGlutDisplay(void)
{
  ::glPointSize(5);
  ::glLineWidth(1);
  {
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0.8, 0.8, 0.8);
    /*
    float color[4] = {200.0/256.0, 200.0/256.0, 200.0/256.0,1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
    ::glEnable(GL_DEPTH_TEST);
     */
//    DrawMeshTri3D_FaceNorm(aXYZ, aTri);
  }
  
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  delfem2::opengl::DrawMeshDynTri3D_Edge(aXYZ, aETri);
}

int main(int argc,char* argv[])
{
  {
    std::vector< std::vector<double> > aaXY;
    aaXY.resize(1);
    double xys[12] = {
      -0.5,-0.5,
      +0.5,-0.5,
      +0.5,+0.5,
      +0.1,+0.6,
      -0.1,+0.6,
      -0.5,+0.5,
    };
    aaXY[0].assign(xys,xys+12);
    GenMesh(aPo2D,aETri,aVec2,
            aaXY, 0.05, 0.05);
  }
  // -------------
  const int np = aPo2D.size();
  aXYZ.resize(np*3);
  for(int ip=0;ip<np;++ip){
    aXYZ[ip*3+0] = aVec2[ip].x();
    aXYZ[ip*3+1] = aVec2[ip].y();
    aXYZ[ip*3+2] = 0.0;
  }
  aXYZt = aXYZ;
  aUVW.resize(np*3,0.0);
  aBCFlag.resize(np,0);
  for(int ip=0;ip<np;++ip){
    if( aXYZ[ip*3+1]  > +0.59 ){
      aBCFlag[ip] = 1;
    }
  }
  aLine.clear();
  { // make aLine
    std::map<int,int> mapY2Ip;
    for(int ip=0;ip<np;++ip){
      if( aXYZ[ip*3+0]  > +0.49 ){
        double y0 = aXYZ[ip*3+1];
        int iy = (int)(y0/0.0132);
        mapY2Ip[iy] = ip;
      }
    }
    for(int ip=0;ip<np;++ip){
      if( aXYZ[ip*3+0]  < -0.49 ){
        double y1 = aXYZ[ip*3+1];
        int iy = (int)(y1/0.0132);
        assert( mapY2Ip.find(iy) != mapY2Ip.end() );
        int ip0 = mapY2Ip[iy];
        aLine.push_back(ip);
        aLine.push_back(ip0);
      }
    }
  }
  
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  delfem2::opengl::setSomeLighting();
  while (!glfwWindowShouldClose(viewer.window))
  {
    {
      static int iframe = 0;
      if( iframe % 5 == 0 ){
        StepTime();
      }
      iframe++;
    }
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
