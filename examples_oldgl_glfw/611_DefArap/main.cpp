/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/defarap.h"
#include "delfem2/mat4.h"
#include "delfem2/thread.h"
#include "delfem2/msh_primitive.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/v3q.h"

namespace dfm2 = delfem2;

// -------------------------------------

void SetPositionAtFixedBoundary(
    std::vector<double> &aRhs,
    unsigned int iframe,
    const std::vector<double> &aXYZ0,
    const std::vector<int> &aBCFlag) {
  double A[16];
  {
    dfm2::CMat4d T0 = dfm2::Mat4_AffineTranslation(0., -0.8, 0.);
    dfm2::CMat4d R0 = dfm2::Mat4_AffineRotationRodriguez(0., +2.0 * sin(0.03 * iframe), 1.0 * sin(0.07 * iframe));
    dfm2::CMat4d T1 = dfm2::Mat4_AffineTranslation(0.2 * sin(0.03 * iframe), +0.5 + 0.1 * cos(0.05 * iframe), 0.);
    (T1 * R0 * T0).CopyTo(A);
  }
  const size_t np = aRhs.size() / 3;
  for (unsigned int ip = 0; ip < np; ++ip) {
    if (aBCFlag[ip * 3 + 0] == 0) { continue; }
    if (aBCFlag[ip * 3 + 0] == 1) {
      aRhs[ip * 3 + 0] = aXYZ0[ip * 3 + 0];
      aRhs[ip * 3 + 1] = aXYZ0[ip * 3 + 1];
      aRhs[ip * 3 + 2] = aXYZ0[ip * 3 + 2];
    }
    if (aBCFlag[ip * 3 + 0] == 2) {
      dfm2::CVec3d t0 = dfm2::Vec3_Mat4Vec3_Homography(A, aXYZ0.data() + ip * 3);
      t0.CopyTo(aRhs.data() + ip * 3);
    }
  }
}

void UpdateRotationsByMatchingCluster_SVD_Parallel(
    std::vector<double> &aQuat1,
    const std::vector<double> &aXYZ0,
    const std::vector<double> &aXYZ1,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup) {
  auto func_matchrot = [&aQuat1, &aXYZ0, & aXYZ1, &psup_ind, &psup](unsigned int ip) {
    dfm2::UpdateRotationsByMatchingCluster_SVD(
        aQuat1,
        ip, aXYZ0, aXYZ1, psup_ind, psup);
  };
  delfem2::parallel_for(aXYZ0.size() / 3, func_matchrot);
}

void myGlutDisplay_Mesh(
    const std::vector<double> &aXYZ0,
    const std::vector<double> &aXYZ1,
    const std::vector<unsigned int> &aTri) {
  ::glLineWidth(1);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1, 0, 0);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1, aTri);
  ::glColor3d(0.8, 0.8, 0.8);
  dfm2::opengl::DrawMeshTri3D_Edge(
      aXYZ0.data(), aXYZ0.size() / 3,
      aTri.data(), aTri.size() / 3);
  ::glColor3d(0, 0, 0);
  dfm2::opengl::DrawMeshTri3D_Edge(
      aXYZ1.data(), aXYZ1.size() / 3,
      aTri.data(), aTri.size() / 3);

}

void Draw_BCFlag(
    const std::vector<double> &aXYZ1,
    const std::vector<int> &aBCFlag) {  // draw bc as a point
  ::glPointSize(10);
  ::glBegin(GL_POINTS);
  for (unsigned int ip = 0; ip < aXYZ1.size() / 3; ++ip) {
    if (aBCFlag[ip * 3 + 0] == 0) { continue; }
    if (aBCFlag[ip * 3 + 0] == 1) { ::glColor3d(0, 0, 1); }
    if (aBCFlag[ip * 3 + 0] == 2) { ::glColor3d(0, 1, 0); }
    ::glVertex3dv(aXYZ1.data() + ip * 3);
  }
  ::glEnd();
}

// --------------------------------------------------

int main() {
  std::vector<unsigned int> tri_vtx;
  std::vector<double> vtx_xyz_ini;
  std::vector<int> dof_bcflag;
  {
    dfm2::MeshTri3D_CylinderClosed(
        vtx_xyz_ini, tri_vtx,
        0.2, 1.6,
        32, 32);
    const size_t np = vtx_xyz_ini.size() / 3;
    dof_bcflag.assign(np * 3, 0);
    for (unsigned int ip = 0; ip < np; ++ip) {
      double y0 = vtx_xyz_ini[ip * 3 + 1];
      if (y0 < -0.65) {
        dof_bcflag[ip * 3 + 0] = 1;
        dof_bcflag[ip * 3 + 1] = 1;
        dof_bcflag[ip * 3 + 2] = 1;
      }
      if (y0 > +0.65) {
        dof_bcflag[ip * 3 + 0] = 2;
        dof_bcflag[ip * 3 + 1] = 2;
        dof_bcflag[ip * 3 + 2] = 2;
      }
    }
  }

  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();

  {
    const unsigned int np = vtx_xyz_ini.size() / 3;
    std::vector<double> vtx_xyz_def = vtx_xyz_ini;
    std::vector<double> vtx_quaternion(np * 4);
    for (unsigned int ip = 0; ip < np; ++ip) {
      dfm2::Quat_Identity(vtx_quaternion.data() + 4 * ip);
    }
    dfm2::Deformer_Arap2 def;
    def.Init(vtx_xyz_ini, tri_vtx, vtx_quaternion, dof_bcflag);

    int iframe = 0;
    glfwSetWindowTitle(viewer.window, "hoge");
    for (; iframe < 200; ++iframe) {
      SetPositionAtFixedBoundary(
          vtx_xyz_def,
          iframe, vtx_xyz_ini, dof_bcflag);
      for(int itr=0;itr<3;++itr) {
        def.Deform(
            vtx_xyz_def, vtx_quaternion,
            vtx_xyz_ini, dof_bcflag);
        UpdateRotationsByMatchingCluster_SVD_Parallel(
            vtx_quaternion,
            vtx_xyz_ini, vtx_xyz_def,
            def.psup_ind, def.psup);
      }
      // --------------------
      viewer.DrawBegin_oldGL();
      myGlutDisplay_Mesh(vtx_xyz_ini, vtx_xyz_def, tri_vtx);
      dfm2::opengl::Draw_QuaternionsCoordinateAxes(vtx_xyz_def, vtx_quaternion, 0.04);
      Draw_BCFlag(vtx_xyz_def, dof_bcflag);
      viewer.SwapBuffers();
      glfwPollEvents();
      viewer.ExitIfClosed();
    }  // end linear displacement only
  }

  for (unsigned int itr = 0; itr < 2; ++itr) {
    const unsigned int np = vtx_xyz_ini.size() / 3;
    std::vector<double> vtx_xyz_def = vtx_xyz_ini;
    std::vector<double> vtx_quaternion(np * 4);
    for (unsigned int ip = 0; ip < np; ++ip) {
      dfm2::Quat_Identity(vtx_quaternion.data() + 4 * ip);
    }
    int iframe = 0;
    { // arap edge linear disponly
      dfm2::Deformer_Arap def0;
      def0.Init(vtx_xyz_ini, tri_vtx, false);
      if (itr == 0) {
        glfwSetWindowTitle(viewer.window, "(1) ARAP without preconditioner");
      } else {
        glfwSetWindowTitle(viewer.window, "(1) ARAP without preconditioner using thread");
      }
      for (; iframe < 200; ++iframe) {
        SetPositionAtFixedBoundary(
            vtx_xyz_def,
            iframe, vtx_xyz_ini, dof_bcflag);
        def0.Deform(
            vtx_xyz_def, vtx_quaternion,
            vtx_xyz_ini, dof_bcflag);
        if (itr == 0) {
          dfm2::UpdateQuaternions_Svd(
              vtx_quaternion,
              vtx_xyz_ini, vtx_xyz_def,
              def0.psup_ind, def0.psup);
        } else {
          UpdateRotationsByMatchingCluster_SVD_Parallel(
              vtx_quaternion,
              vtx_xyz_ini, vtx_xyz_def, def0.psup_ind, def0.psup);
        }
        // --------------------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(vtx_xyz_ini, vtx_xyz_def, tri_vtx);
        dfm2::opengl::Draw_QuaternionsCoordinateAxes(vtx_xyz_def, vtx_quaternion, 0.04);
        Draw_BCFlag(vtx_xyz_def, dof_bcflag);
        viewer.SwapBuffers();
        glfwPollEvents();
        viewer.ExitIfClosed();
      }
    }  // end linear displacement only
    {  // arap edge linear displacement only
      dfm2::Deformer_Arap def0;
      def0.Init(vtx_xyz_ini, tri_vtx, true);
      if (itr == 0) {
        glfwSetWindowTitle(viewer.window, "(2) ARAP with preconditioner");
      } else {
        glfwSetWindowTitle(viewer.window, "(2) ARAP with preconditioner using thread");
      }
      for (; iframe < 400; ++iframe) {
        SetPositionAtFixedBoundary(
            vtx_xyz_def,
            iframe, vtx_xyz_ini, dof_bcflag);
        def0.Deform(
            vtx_xyz_def, vtx_quaternion,
            vtx_xyz_ini, dof_bcflag);
        if (itr == 0) {
          dfm2::UpdateQuaternions_Svd(
              vtx_quaternion,
              vtx_xyz_ini, vtx_xyz_def, def0.psup_ind, def0.psup);
        } else {
          UpdateRotationsByMatchingCluster_SVD_Parallel(
              vtx_quaternion,
              vtx_xyz_ini, vtx_xyz_def, def0.psup_ind, def0.psup);
        }
        // --------------------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(vtx_xyz_ini, vtx_xyz_def, tri_vtx);
        dfm2::opengl::Draw_QuaternionsCoordinateAxes(vtx_xyz_def, vtx_quaternion, 0.04);
        Draw_BCFlag(vtx_xyz_def, dof_bcflag);
        viewer.SwapBuffers();
        glfwPollEvents();
        viewer.ExitIfClosed();
      }
    }  // end linear disponly
  }
}


