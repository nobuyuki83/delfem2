/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/femsolidhyper.h"
#include "delfem2/femutil.h"
namespace dfm2 = delfem2;

#include <random>
#include <vector>
#include <cstdlib>
#include <GLFW/glfw3.h>

// --------------------------------------------------------------

/*
void Draw()
{
  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  {
    float color[4] = {200.0/256.0, 200.0/256.0, 200.0/256.0,1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
//    glShadeModel(GL_SMOOTH);
    glShadeModel(GL_FLAT);
 }
  ::glDisable(GL_LIGHTING);
  
  {
    ::glColor3d(0,0,0);
    delfem2::opengl::DrawMeshTet3D_EdgeDisp(aXYZ.data(),
                                            aTet.data(),aTet.size()/4,
                                            aMode.data(),
                                            0.1);
    ::glEnable(GL_LIGHTING);
    {
      float color[4] = {180.0/256.0, 180.0/256.0, 130.0/256.0,1.0f};
      ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
      ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
      glShadeModel(GL_FLAT);
    }
    
    delfem2::opengl::DrawMeshTet3D_FaceNorm(aXYZ.data(),
                                            aTet.data(), aTet.size()/4);

  }

  
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }
  else{              ::glDisable(GL_LIGHTING); }
}
*/

void Simulation(
    std::vector<double>& aXYZ,
    std::vector<unsigned int>& aHex,
    std::vector<double>& aDisp)
{
  double c1 = 1.0;
  double c2 = 1.0;
  for(unsigned int ih=0;ih<aHex.size()/8;++ih){
    double aP0[8][3]; delfem2::FetchData<8,3>(aP0,aHex.data()+ih*8,aXYZ.data());
    double aU[8][3]; delfem2::FetchData<8,3>(aU,aHex.data()+ih*8,aDisp.data());
    //
    double W = 0.0;
    double dW[8][3]; for(int i=0;i<8*3;++i){ (&dW[0][0])[i] = 0.0; }
    double ddW[8][8][3][3]; for(int i=0;i<8*8*3*3;++i){ (&ddW[0][0][0][0])[i] = 0.0; }
    double vol = 0.0;
    delfem2::AddWdWddW_Solid3HyperMooneyrivlin2Reduced_Hex(
        W, dW, ddW, vol,
        c1, c2, aP0, aU, 1);
    //
    std::cout << ih << " " << vol << std::endl;
  }
}

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aHex;
  dfm2::MeshHex3_Grid(
      aXYZ0, aHex,
      5, 5, 1, 0.2);

  std::vector<double> aDisp;
  aDisp.resize(aXYZ0.size()/3,0.0);

  Simulation(aXYZ0,aHex,aDisp);

  // ----------------------
  delfem2::glfw::CViewer3 viewer;
  viewer.camera.view_height = 1.5;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();

  delfem2::opengl::setSomeLighting();
  while(!glfwWindowShouldClose(viewer.window)){
    // -----
    ::glDisable(GL_LIGHTING);
    viewer.DrawBegin_oldGL();
    ::glColor3d(0,0,0);
    ::glPointSize(3);
    dfm2::opengl::DrawPoints3_Points(aXYZ0);
    //
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::DrawMeshHex3D_FaceNorm(aXYZ0.data(), aHex.data(), aHex.size() / 8);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
