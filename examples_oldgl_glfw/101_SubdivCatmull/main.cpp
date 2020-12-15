/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @brief simple demo of subdivision surface
 */

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshsubdiv.h"
#include "delfem2/primitive.h"
#include <GLFW/glfw3.h>
#include <cstdlib>

namespace dfm2 = delfem2;

// --------------------------------------------------

int main(int argc,char* argv[])
{
  // --------------------------
  // begin computing subdivision

  std::vector< std::vector<unsigned int> > aaQuad;
  std::vector< std::vector<double> > aaXYZ;
  const unsigned int nlevel_subdiv = 3;
  
  aaXYZ.resize(nlevel_subdiv+1);
  aaQuad.resize(nlevel_subdiv+1);
  {
    const double bbmin[3] = {-1,-1,-1};
    const double bbmax[3] = {+1,+1,+1};
    delfem2::MeshQuad3_CubeVox(
        aaXYZ[0],aaQuad[0],
        bbmin, bbmax);
  }
  for(unsigned int il=0;il<nlevel_subdiv;++il){
    const std::vector<double>& aXYZ0 = aaXYZ[il];
    const std::vector<unsigned int>& aQuad0 = aaQuad[il];
    std::vector<unsigned int>& aQuad1 = aaQuad[il+1];
    std::vector<int> aEdgeFace0;
    std::vector<unsigned int> psupIndQuad0, psupQuad0;
    dfm2::SubdivTopo_MeshQuad(
        aQuad1,
        psupIndQuad0,psupQuad0, aEdgeFace0,
        aQuad0.data(), (unsigned int)(aQuad0.size()/4),
        (unsigned int)(aXYZ0.size()/3));
    std::vector<double>& aXYZ1 = aaXYZ[il+1];
    delfem2::SubdivisionPoints_QuadCatmullClark(
        aXYZ1,
        aQuad1,aEdgeFace0,psupIndQuad0,psupQuad0,
        aQuad0.data(),
        (unsigned int)aQuad0.size()/4,
        aXYZ0.data(),
        (unsigned int)aXYZ0.size()/3);
  }
  
  // end computing subdivision
  // ----------------------------

  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.camera.view_height = 2.0;
  delfem2::opengl::setSomeLighting();
  while (!glfwWindowShouldClose(viewer.window))
  {
    viewer.DrawBegin_oldGL();
    ::glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
    ::glEnable(GL_LIGHTING);
    delfem2::opengl::DrawMeshQuad3D_FaceNorm(aaXYZ[nlevel_subdiv],aaQuad[nlevel_subdiv]);
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0,0,0);
    delfem2::opengl::DrawMeshQuad3D_Edge(aaXYZ[nlevel_subdiv],aaQuad[nlevel_subdiv]);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}