/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <string>
#include <cstdlib>
#include <set>
#include "delfem2/mshmisc.h"
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/clusterpoints.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// -----------------------------




void DrawConnectedPoints(
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& psup_ind,
    const std::vector<unsigned int>& psup)
{
  unsigned int np = aXYZ.size()/3;
  assert(psup_ind.size()==np+1);
  ::glBegin(GL_LINES);
  for(unsigned int ip=0;ip<np;++ip){
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      unsigned int jp = psup[ipsup];
      ::glVertex3dv(aXYZ.data()+ip*3);
      ::glVertex3dv(aXYZ.data()+jp*3);
    }
  }
  ::glEnd();
}

// -----------------------------

int main(int argc,char* argv[])
{
  class CPointData
  {
  public:
    std::vector<double> aXYZ;
    std::vector<double> aArea;
    std::vector<double> aNorm;
    std::vector<unsigned int> psup_ind, psup;
  };

  std::vector<CPointData> aPointData;

  aPointData.resize(1);
  {
    std::vector<unsigned int> aTri0;
    {
      CPointData &pd0 = aPointData[0];
      delfem2::Read_Ply(std::string(PATH_INPUT_DIR) + "/bunny_2k.ply",
                        pd0.aXYZ, aTri0);
      dfm2::Normalize_Points3(pd0.aXYZ, 2.0);
    }
    { // make normal
      CPointData &pd0 = aPointData[0];
      pd0.aNorm.resize(pd0.aXYZ.size());
      dfm2::Normal_MeshTri3D(pd0.aNorm.data(),
                             pd0.aXYZ.data(), pd0.aXYZ.size() / 3,
                             aTri0.data(), aTri0.size() / 3);
    }
    { // make area
      CPointData &pd0 = aPointData[0];
      pd0.aArea.resize(pd0.aXYZ.size() / 3);
      dfm2::MassPoint_Tri3D(pd0.aArea.data(),
                            1.0,
                            pd0.aXYZ.data(), pd0.aXYZ.size() / 3,
                            aTri0.data(), aTri0.size() / 3);
    }
    { // make psup
      CPointData &pd0 = aPointData[0];
      dfm2::JArray_PSuP_MeshElem(pd0.psup_ind, pd0.psup,
                                 aTri0.data(), aTri0.size() / 3, 3,
                                 pd0.aXYZ.size() / 3);
    }
  }

  for(unsigned int itr=0;itr<8;++itr) {
    aPointData.resize(aPointData.size() + 1);
    const CPointData& pd0 = aPointData[itr];
    CPointData& pd1 = aPointData[itr+1];
    std::vector<unsigned int> map01;
    dfm2::BinaryClusteringPoints(
        pd1.aXYZ, pd1.aArea, pd1.aNorm, map01,
        pd0.aXYZ, pd0.aArea, pd0.aNorm, pd0.psup_ind, pd0.psup);
    dfm2::BinaryClusteringPoints_FindConnection(pd1.psup_ind,pd1.psup,
        pd1.aXYZ.size()/3,map01,pd0.psup_ind,pd0.psup);
  }

  // -----------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.nav.camera.camera_rot_mode = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.5;
  while (!glfwWindowShouldClose(viewer.window))
  {
    for(const auto& pd: aPointData) {
//      const auto& pd = aPointData[5];
      for(unsigned int itr=0;itr<30;++itr) {
        viewer.DrawBegin_oldGL();
        ::glColor3d(0, 0, 0);
        DrawConnectedPoints(pd.aXYZ, pd.psup_ind, pd.psup);
        viewer.DrawEnd_oldGL();
      }
    }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
