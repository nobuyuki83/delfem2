/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @brief simple demo of subdivision surface
 */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"
#include "delfem2/primitive.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/gl2_funcs.h"

int main(int argc,char* argv[])
{
  CViewer_GLFW viewer;
  viewer.Init_GLold();
  
  // --------------------------
  // begin computing subdivision

  std::vector< std::vector<unsigned int> > aaQuad;
  std::vector< std::vector<double> > aaXYZ;
  const unsigned int nlevel_subdiv = 3;
  
  aaXYZ.resize(nlevel_subdiv+1);
  aaQuad.resize(nlevel_subdiv+1);
  MeshQuad3D_CubeVox(aaXYZ[0],aaQuad[0],
                     -1,+1,  -1,+1,  -1,+1);
  for(unsigned int il=0;il<nlevel_subdiv;++il){
    const std::vector<double>& aXYZ0 = aaXYZ[il];
    const std::vector<unsigned int>& aQuad0 = aaQuad[il];
    std::vector<unsigned int>& aQuad1 = aaQuad[il+1];
    std::vector<int> aEdgeFace0;
    std::vector<int> psupIndQuad0, psupQuad0;
    QuadSubdiv(aQuad1,
               psupIndQuad0,psupQuad0, aEdgeFace0,
               aQuad0.data(), aQuad0.size()/4,
               aXYZ0.size()/3);
    std::vector<double>& aXYZ1 = aaXYZ[il+1];
    SubdivisionPoints_QuadCatmullClark(aXYZ1,
                                       aQuad1,aEdgeFace0,psupIndQuad0,psupQuad0,
                                       aQuad0.data(), aQuad0.size()/4,
                                       aXYZ0.data(),  aXYZ0.size()/3);
  }
  
  // end computing subdivision
  // ----------------------------
  
  viewer.nav.camera.view_height = 2.0;
  setSomeLighting();
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    viewer.DrawBegin_Glold();
    
    ::glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
    ::glEnable(GL_LIGHTING);
    DrawMeshQuad3D_FaceNorm(aaXYZ[nlevel_subdiv],aaQuad[nlevel_subdiv]);
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0,0,0);
    DrawMeshQuad3D_Edge(aaXYZ[nlevel_subdiv],aaQuad[nlevel_subdiv]);
    
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


