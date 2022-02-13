/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/msh_primitive.h"
#include "delfem2/mshmisc.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/vec3_funcs.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/rigv3.h"
#include "delfem2/opengl/old/mshuni.h"

#ifndef M_PI
#  define M_PI 3.141592
#endif

namespace dfm2 = delfem2;

// ---------------------------------------------------

int main() {
  std::vector<double> vtx_xyz_ini;
  std::vector<unsigned int> tri_vtx;
  // -----
  dfm2::MeshTri3_Capsule(
      vtx_xyz_ini, tri_vtx,
      0.2, 1.0,
      16, 5, 8);
  std::vector<dfm2::CRigBone> bones;
  {
    dfm2::CRigBone b;
    b.ibone_parent = -1;
    b.invBindMat[7] = +0.5;
    b.transRelative[1] = -0.5;
    bones.push_back(b);
  }
  {
    dfm2::CRigBone b;
    b.ibone_parent = 0;
    b.transRelative[1] = +0.5;
    bones.push_back(b);
  }
  std::vector<double> skinning_weight;
  {
    const size_t np = vtx_xyz_ini.size() / 3;
    const size_t nb = bones.size();
    skinning_weight.resize(np * nb);
    for (unsigned int ip = 0; ip < np; ++ip) {
      const double *p0 = vtx_xyz_ini.data() + ip * 3;
      double w_tot = 0;
      for (unsigned int ib = 0; ib < nb; ++ib) {
        double pb[3] = {
            -bones[ib].invBindMat[3],
            -bones[ib].invBindMat[7],
            -bones[ib].invBindMat[11]};
        double len = dfm2::Distance3(p0, pb);
        double wb = 1.0 / (len + 1.0e-10);
        skinning_weight[ip * nb + ib] = wb;
        w_tot += wb;
      }
      for (unsigned int ib = 0; ib < nb; ++ib) {
        skinning_weight[ip * nb + ib] /= w_tot;
      }
    }
  }
  // ------
  std::vector<double> vtx_xyz_def = vtx_xyz_ini;

  // ----------------
  dfm2::glfw::CViewer3 viewer;
  //
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  dfm2::opengl::setSomeLighting();

  int iframe = 0;
  while (true) {
    ++iframe;
    {
      dfm2::Quat_Bryant(
          bones[1].quatRelativeRot,
          0., 0., 0.8 * sin(0.1 * iframe));
      dfm2::UpdateBoneRotTrans(bones);
      dfm2::Skinning_LBS(
          vtx_xyz_def,
          vtx_xyz_ini, bones, skinning_weight);
    }
    // -----
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_DEPTH_TEST);
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(vtx_xyz_def, tri_vtx);
    ::glColor3d(0, 0, 0);
    dfm2::opengl::DrawMeshTri3D_Edge(vtx_xyz_def, tri_vtx);
    ::glDisable(GL_DEPTH_TEST);
    ::glDisable(GL_LIGHTING);
    ::glColor3d(1, 0, 0);
    dfm2::opengl::DrawBone_Line(
        bones,
        -1, 0, 0.02, 0.2);
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { break; }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
