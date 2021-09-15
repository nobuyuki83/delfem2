/*
 * Copyright (c) 2019-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <algorithm>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/points.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mshmisc.h"
#include "delfem2/srch_v3bvhmshtopo.h"
#include "delfem2/srchbv3sphere.h"
#include "delfem2/srchbvh.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"

namespace dfm2 = delfem2;

// ----------------------------------------

void myGlutDisplay(
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& aTri,
    const std::vector<unsigned int>& aIndElm,
    const double src0[3],
    const double dir0[3])
{
  ::glDisable(GL_LIGHTING);
  //
  ::glPointSize(2);
  ::glBegin(GL_POINTS);
  ::glColor3d(0,0,0);
  for(size_t ip=0;ip<aXYZ.size()/3;++ip){
    ::glVertex3d(aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]);
  }
  ::glEnd();
  //
  ::glBegin(GL_TRIANGLES);
  ::glColor3d(0,0,1);
  for(unsigned int it : aIndElm){
    ::glVertex3dv(aXYZ.data()+aTri[it*3+0]*3);
    ::glVertex3dv(aXYZ.data()+aTri[it*3+1]*3);
    ::glVertex3dv(aXYZ.data()+aTri[it*3+2]*3);
  }
  ::glEnd();
  //
  ::glPointSize(4);
  ::glBegin(GL_POINTS);
  ::glColor3d(1,0,0);
  ::glVertex3dv(src0);
  ::glEnd();
  //
  ::glColor3d(1,0,0);
  const double L = 100;
  ::glBegin(GL_LINES);
  ::glVertex3d(src0[0]+L*dir0[0],src0[1]+L*dir0[1],src0[2]+L*dir0[2]);
  ::glVertex3d(src0[0]-L*dir0[0],src0[1]-L*dir0[1],src0[2]-L*dir0[2]);
  ::glEnd();
}

int main()
{
  std::vector<double> vec_xyz; // 3d points
  std::vector<unsigned int> vec_tri;

  { // load input mesh
    delfem2::Read_Ply(
        vec_xyz, vec_tri,
        std::filesystem::path(PATH_INPUT_DIR) / "bunny_2k.ply");
    dfm2::Normalize_Points3(vec_xyz, 2.0);
  }

  std::vector<dfm2::CNodeBVH2> vec_node_bvh;
  std::vector<dfm2::CBV3_Sphere<double>> vec_bv;
  dfm2::ConstructBVHTriangleMeshMortonCode(
      vec_node_bvh, vec_bv, vec_xyz, vec_tri);
  
  dfm2::glfw::CViewer3 viewer;
  viewer.camera.view_height = 1.5;
  viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  delfem2::opengl::setSomeLighting();

  double cur_time = 0.0;
  while (!glfwWindowShouldClose(viewer.window))
  {
    cur_time += 0.001;
    const double src0[3] = {
        0.5*sin(cur_time*1),
        0.5*sin(cur_time*2),
        0.5*sin(cur_time*3) };
    const double dir0[3] = {
        0.5*sin(cur_time*5),
        0.5*sin(cur_time*6),
        0.5*sin(cur_time*7) };
    // -----------
    std::vector<unsigned int> vec_tri_index;
    dfm2::BVH_GetIndElem_Predicate(
        vec_tri_index,
        dfm2::CIsBV_IntersectLine< dfm2::CBV3_Sphere<double>, double>(src0,dir0),
        0, vec_node_bvh, vec_bv);
    // -----------
    viewer.DrawBegin_oldGL();
    myGlutDisplay(vec_xyz, vec_tri, vec_tri_index, src0, dir0);
    viewer.SwapBuffers();
    
    glfwPollEvents();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


