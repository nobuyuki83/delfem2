/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <cmath>
#include <cstdlib>
#include "delfem2/primitive.h"
#include "delfem2/mshmisc.h"
#include "delfem2/rig_geo3.h"

// gl related includes
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/rigv3_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

#ifndef M_PI
#  define M_PI 3.141592
#endif

namespace dfm2 = delfem2;

// ---------------------------------------------------

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aElm;
  // -----
  dfm2::MeshTri3_Capsule(aXYZ0,aElm,
                         0.2, 1.0,
                         16, 5, 8);
  std::vector<dfm2::CRigBone> aBone;
  {
    dfm2::CRigBone b;
    b.ibone_parent = -1;
    b.invBindMat[7] = +0.5;
    b.transRelative[1] = -0.5;
    aBone.push_back(b);
  }
  {
    dfm2::CRigBone b;
    b.ibone_parent = 0;
    b.transRelative[1] = +0.5;
    aBone.push_back(b);
  }
  std::vector<double> aW;
  {
    const unsigned int np = aXYZ0.size()/3;
    const unsigned int nb = aBone.size();
    aW.resize(np*nb);
    for(unsigned int ip=0;ip<aXYZ0.size()/3;++ip) {
      const double* p0 = aXYZ0.data()+ip*3;
      double w_tot = 0;
      for(unsigned int ib=0;ib<nb;++ib){
        double pb[3] = {
            -aBone[ib].invBindMat[3],
            -aBone[ib].invBindMat[7],
            -aBone[ib].invBindMat[11]};
        double len = dfm2::Distance3(p0,pb);
        double wb = 1.0/(len+1.0e-10);
        aW[ip*nb+ib] = wb;
        w_tot += wb;
      }
      for(unsigned int ib=0;ib<nb;++ib) {
        aW[ip*nb+ib] /= w_tot;
      }
    }
  }
  // ------
  std::vector<double> aXYZ1 = aXYZ0;

  // ----------------
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  viewer.nav.camera.Rot_Camera(+0.0,+0.0);
  dfm2::opengl::setSomeLighting();

  int iframe = 0;
  while (true)
  {
    ++iframe;
    {
      dfm2::Quat_Bryant(aBone[1].quatRelativeRot,  0.,0.,0.8*sin(0.1*iframe));
      dfm2::UpdateBoneRotTrans(aBone);
      dfm2::Skinning_LBS(aXYZ1,
                         aXYZ0, aBone, aW);
    }
    // -----
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_DEPTH_TEST);
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1, aElm);
    ::glColor3d(0,0,0);
    dfm2::opengl::DrawMeshTri3D_Edge(aXYZ1,aElm);
    ::glDisable(GL_DEPTH_TEST);
    ::glDisable(GL_LIGHTING);
    ::glColor3d(1,0,0);
    dfm2::opengl::DrawBone(aBone,-1,0,0.02,0.2);
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
