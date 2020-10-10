/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <algorithm>
#include "delfem2/points.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/bv.h"
#include "delfem2/bvh.h"
// ---------------------------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

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

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ; // 3d points
  std::vector<unsigned int> aTri;

  { // load input mesh
    delfem2::Read_Ply(std::string(PATH_INPUT_DIR) + "/bunny_2k.ply",
                      aXYZ, aTri);
    dfm2::Normalize_Points3(aXYZ, 2.0);
  }

  std::vector<dfm2::CNodeBVH2> aNodeBVH;
  std::vector<dfm2::CBV3_Sphere<double>> aAABB;
  {
    std::vector<double> aCent;
    double rad = dfm2::CentsMaxRad_MeshTri3(aCent,
        aXYZ,aTri);
    double min_xyz[3], max_xyz[3];
    delfem2::BoundingBox3_Points3(min_xyz,max_xyz,
        aCent.data(), aCent.size()/3);
    std::vector<unsigned int> aSortedId;
    std::vector<std::uint32_t> aSortedMc;
    dfm2::SortedMortenCode_Points3(aSortedId,aSortedMc,
                                   aCent,min_xyz,max_xyz);
    dfm2::BVHTopology_Morton(aNodeBVH,
                             aSortedId,aSortedMc);

    dfm2::CLeafVolumeMaker_Mesh<dfm2::CBV3_Sphere<double>,double> lvm(
        1.0e-10,
        aXYZ.data(), aXYZ.size()/3,
        aTri.data(), aTri.size()/3, 3);
    dfm2::BVH_BuildBVHGeometry(aAABB,
        0, aNodeBVH,
        lvm);
    { // assertion
      dfm2::Check_MortonCode_Sort(aSortedId, aSortedMc, aCent, min_xyz, max_xyz);
      dfm2::Check_MortonCode_RangeSplit(aSortedMc);
      dfm2::Check_BVH(aNodeBVH,aCent.size()/3);
    }
  }
  
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.nav.camera.view_height = 1.5;
  viewer.nav.camera.camera_rot_mode = dfm2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  viewer.Init_oldGL();
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
    std::vector<unsigned int> aIndElem;
    dfm2::BVH_GetIndElem_IntersectLine(aIndElem,
      src0, dir0,
      0, aNodeBVH, aAABB);
    // -----------
    viewer.DrawBegin_oldGL();
    myGlutDisplay(aXYZ,aTri,aIndElem,src0,dir0);
    viewer.DrawEnd_oldGL();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


