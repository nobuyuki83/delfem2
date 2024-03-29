/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <random>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/defarap.h"
#include "delfem2/mat4.h"
#include "delfem2/msh_primitive.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/mshuni.h"

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

void myGlutDisplay_Mesh(
    const std::vector<double> &aXYZ0,
    const std::vector<double> &aXYZ1,
    const std::vector<unsigned int> &aTri) {
  ::glLineWidth(1);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1, 0, 0);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1, aTri);
  ::glColor3d(0.8, 0.8, 0.8);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ0.data(), aXYZ0.size() / 3,
                                   aTri.data(), aTri.size() / 3);
  ::glColor3d(0, 0, 0);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ1.data(), aXYZ1.size() / 3,
                                   aTri.data(), aTri.size() / 3);

}

void Draw_BCFlag(
    const std::vector<double> &aXYZ1,
    const std::vector<int> &aBCFlag) { // draw bc as a point
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
        16, 16);
    {
      const auto np = static_cast<unsigned int>(vtx_xyz_ini.size() / 3);
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
  }
  std::vector<double> vtx_xyz_def = vtx_xyz_ini;

  // ---------------

  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  //
  for (unsigned int itr = 0; itr < 3; ++itr) {
    const double weight_bc = 100.0;
    int iframe = 0;
    { // arap edge linear disponly
      dfm2::CDef_ArapEdgeLinearDisponly def0(vtx_xyz_ini, tri_vtx, weight_bc, dof_bcflag);
      glfwSetWindowTitle(viewer.window, "(1) ARAP Edge Linear Disponly");
      for (; iframe < 50; ++iframe) {
        SetPositionAtFixedBoundary(
            vtx_xyz_def,
            iframe, vtx_xyz_ini, dof_bcflag);
        def0.Deform(vtx_xyz_def,
                    vtx_xyz_ini);
        // --------------------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(vtx_xyz_ini, vtx_xyz_def, tri_vtx);
        Draw_BCFlag(vtx_xyz_def, dof_bcflag);
        viewer.SwapBuffers();
        glfwPollEvents();
        viewer.ExitIfClosed();
      }
    } // end linear disponly
    // -------------------------------------------------------
    { // begin lienar disprot without preconditioner
      glfwSetWindowTitle(viewer.window, "(2) Arap Edge Linear Disprot without Prec");
      auto np = static_cast<unsigned int>(vtx_xyz_ini.size() / 3);
      dfm2::CDef_ArapEdge def1;
      def1.Init(vtx_xyz_ini, tri_vtx, weight_bc, dof_bcflag, false);
      std::vector<double> aQuat(np * 4); // array of quaternion
      for (; iframe < 100; ++iframe) {
        for (unsigned int ip = 0; ip < np; ++ip) { dfm2::Quat_Identity(aQuat.data() + ip * 4); } // Initialize
        SetPositionAtFixedBoundary(
            vtx_xyz_def,
            iframe, vtx_xyz_ini, dof_bcflag);
        def1.Deform(vtx_xyz_def, aQuat,
                    vtx_xyz_ini);
        // ------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(vtx_xyz_ini, vtx_xyz_def, tri_vtx);
        Draw_BCFlag(vtx_xyz_def, dof_bcflag);
        dfm2::opengl::Draw_QuaternionsCoordinateAxes(vtx_xyz_def, aQuat, 0.04);
        viewer.SwapBuffers();
        glfwPollEvents();
        viewer.ExitIfClosed();
      } // end of frame loop
    } // end linear disprot without preconditioner
    // -------------------------------
    { // begin lienar disprot with preconditioner
      glfwSetWindowTitle(viewer.window, "(3) Arap Edge Linear Disprot with Prec");
      const auto np = static_cast<unsigned int>(vtx_xyz_ini.size() / 3);
      dfm2::CDef_ArapEdge def1;
      def1.Init(vtx_xyz_ini, tri_vtx, weight_bc, dof_bcflag, true);
      std::vector<double> aQuat(np * 4);
      for (; iframe < 200; ++iframe) {
        for (unsigned int ip = 0; ip < np; ++ip) { dfm2::Quat_Identity(aQuat.data() + ip * 4); } // initialize
        for (unsigned int i = 0; i < np * 3; ++i) { // adding noise for debuggng purpose
          vtx_xyz_def[i] += 0.02 * (double) rand() / (RAND_MAX + 1.0) - 0.01;
        }
        SetPositionAtFixedBoundary(
            vtx_xyz_def,
            iframe, vtx_xyz_ini, dof_bcflag);
        def1.Deform(
            vtx_xyz_def, aQuat,
            vtx_xyz_ini);
        // ------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(vtx_xyz_ini, vtx_xyz_def, tri_vtx);
        Draw_BCFlag(vtx_xyz_def, dof_bcflag);
        dfm2::opengl::Draw_QuaternionsCoordinateAxes(vtx_xyz_def, aQuat, 0.04);
        viewer.SwapBuffers();
        glfwPollEvents();
        viewer.ExitIfClosed();
      } // end of frame loop
    } // end linear disprot with preconditioner
    // -------------------------------
    { // begin nonlienar disprot with preconditioner
      glfwSetWindowTitle(viewer.window, "(4) Arap Edge NonLinear Disprot with Prec");
      const auto np = static_cast<unsigned int>(vtx_xyz_ini.size() / 3);
      dfm2::CDef_ArapEdge def1;
      def1.Init(vtx_xyz_ini, tri_vtx, weight_bc, dof_bcflag, true);
      std::vector<double> vtx_quaternion(np * 4);
      for (unsigned int ip = 0; ip < np; ++ip) { dfm2::Quat_Identity(vtx_quaternion.data() + ip * 4); }
      for (; iframe < 400; ++iframe) {
        for (unsigned int i = 0; i < np * 3; ++i) { // adding noise for debuggng purpose
          vtx_xyz_def[i] += 0.02 * (double) rand() / (RAND_MAX + 1.0) - 0.01;
        }
        SetPositionAtFixedBoundary(
            vtx_xyz_def,
            iframe, vtx_xyz_ini, dof_bcflag);
        def1.Deform(
            vtx_xyz_def, vtx_quaternion,
            vtx_xyz_ini);
        // ------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(vtx_xyz_ini, vtx_xyz_def, tri_vtx);
        Draw_BCFlag(vtx_xyz_def, dof_bcflag);
        dfm2::opengl::Draw_QuaternionsCoordinateAxes(vtx_xyz_def, vtx_quaternion, 0.04);
        viewer.SwapBuffers();
        glfwPollEvents();
        viewer.ExitIfClosed();
      } // end of frame loop
    } // end linear disprot with preconditioner
  }
}


