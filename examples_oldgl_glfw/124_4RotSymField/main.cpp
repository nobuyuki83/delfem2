/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief implementation of 4 rotatoinal symetry field
 * @details implementation is based on "Wenzel Jakob, Marco Tarini, Daniele Panozzo, and Olga Sorkine-Hornung. Instant field-aligned meshes. Siggraph Asia 2015"
 */

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vec3.h"
#include "delfem2/4rotsym.h"
#include <GLFW/glfw3.h>
#include <cstdlib>

namespace dfm2 = delfem2;

// ------------------------------------------------
// TODO: Add random permutation version in the demo

int main()
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  delfem2::Read_Ply(std::string(PATH_INPUT_DIR)+"/bunny_1k.ply",
      aXYZ,aTri);
  delfem2::Normalize_Points3(aXYZ);
  std::vector<double> aNorm(aXYZ.size());
  dfm2::Normal_MeshTri3D(aNorm.data(),
      aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
  
  std::vector<double> aOdir;
  {
    const double minCoords[3] =  {-1., -1., -1.};
    const double maxCoords[3] =  {+1., +1., +1.};
    aOdir.resize(aXYZ.size());
    dfm2::Points_RandomUniform(aOdir.data(),
        aXYZ.size() / 3, 3, minCoords, maxCoords);
    dfm2::TangentVector_Points3(aOdir,
        aNorm);
  }
  
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
      aTri.data(), aTri.size()/3, 3, aXYZ.size()/3);
  
  // ------------------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  dfm2::opengl::setSomeLighting();
  unsigned int iframe = 0;
  while (true)
  {
    if( iframe == 0 ){
      const double minCoords[3] =  {-1., -1., -1.};
      const double maxCoords[3] =  {+1., +1., +1.};
      aOdir.resize(aXYZ.size());
      dfm2::Points_RandomUniform(aOdir.data(),
                                 aXYZ.size() / 3, 3, minCoords, maxCoords);
      dfm2::TangentVector_Points3(aOdir,aNorm);
    }
    if( iframe > 30 ){
      dfm2::Smooth4RotSym(aOdir,
          aNorm, psup_ind, psup);
    }
    // --------------------
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ.data(), aTri.data(), aTri.size()/3);
    {
      ::glDisable(GL_LIGHTING);
      double len = 0.03;
      ::glLineWidth(3);
      unsigned int np = aXYZ.size()/3;
      for(unsigned int ip=0;ip<np;++ip){
        const dfm2::CVec3d p = dfm2::CVec3d(aXYZ.data()+ip*3);
        const dfm2::CVec3d n = dfm2::CVec3d(aNorm.data()+ip*3).Normalize();
        const dfm2::CVec3d o = dfm2::CVec3d(aOdir.data()+ip*3).Normalize();
        const dfm2::CVec3d q = dfm2::Cross(n,o);
        ::glBegin(GL_LINES);
        ::glColor3d(0,0,0);
        dfm2::opengl::myGlVertex(p);
        dfm2::opengl::myGlVertex(p+len*n);
        ::glColor3d(0,0,1);
        dfm2::opengl::myGlVertex(p-len*o);
        dfm2::opengl::myGlVertex(p);
        ::glColor3d(1,0,0);
        dfm2::opengl::myGlVertex(p);
        dfm2::opengl::myGlVertex(p+len*o);
        dfm2::opengl::myGlVertex(p-len*q);
        dfm2::opengl::myGlVertex(p+len*q);
        ::glEnd();
      }
    }
    iframe = (iframe+1)%100;
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
