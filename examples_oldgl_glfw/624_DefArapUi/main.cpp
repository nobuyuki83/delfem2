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

#include "delfem2/defarap.h"
#include "delfem2/gizmo_geo3.h"
#include "delfem2/msh_primitive.h"
#include "delfem2/vec3.h"
#include "delfem2/quat.h"
#include "delfem2/mat4.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/gizmo.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// ----------------------------------------

int main() {
  class MyViewer : public delfem2::glfw::CViewer3 {
   public:
    MyViewer() {
      dfm2::MeshTri3D_CylinderClosed(
          aXYZ0, aTri,
          0.2, 1.6,
          32, 32);
      const auto np = static_cast<unsigned int>(aXYZ0.size() / 3);
      aBCFlag.assign(np * 3, 0);
      for (unsigned int ip = 0; ip < np; ++ip) {
        double y0 = aXYZ0[ip * 3 + 1];
        if (y0 < -0.65) {
          aBCFlag[ip * 3 + 0] = 1;
          aBCFlag[ip * 3 + 1] = 1;
          aBCFlag[ip * 3 + 2] = 1;
        }
        if (y0 > +0.65) {
          aBCFlag[ip * 3 + 0] = 2;
          aBCFlag[ip * 3 + 1] = 2;
          aBCFlag[ip * 3 + 2] = 2;
        }
      }
      { // initialize gizmo
        unsigned int nbc = 2;
        std::vector<dfm2::CVec3d> aCntBC;
        aCntBC.assign(nbc, dfm2::CVec3d(0, 0, 0));
        std::vector<unsigned int> aW(nbc, 0);
        for (unsigned int ip = 0; ip < np; ++ip) {
          if (aBCFlag[ip * 3 + 0] <= 0) { continue; }
          unsigned int ibc = aBCFlag[ip * 3 + 0] - 1;
          aCntBC[ibc] += dfm2::CVec3d(aXYZ0.data() + ip * 3);
          aW[ibc] += 1;
        }
        for (unsigned int ibc = 0; ibc < nbc; ++ibc) {
          aCntBC[ibc] /= (double) aW[ibc];
        }
        giz1.pivot0 = aCntBC[1].cast<float>();
        giz1.gizmo_rot.pos = aCntBC[1].cast<float>();
        giz1.gizmo_trnsl.pos = aCntBC[1].cast<float>();
        giz1.gizmo_rot.size = 0.3f;
        giz1.gizmo_trnsl.size = 0.3f;
      }
      aXYZ1 = aXYZ0;
      aQuat1.resize(np * 4);
      for (unsigned int ip = 0; ip < np; ++ip) {
        dfm2::Quat_Identity(aQuat1.data() + 4 * ip);
      }
      def0.Init(aXYZ0, aTri, true);
    }
    void mouse_press(const float src[3], const float dir[3]) override {
      giz1.Pick(src, dir);
    }
    void mouse_drag(const float src0[3], const float src1[3], const float dir[3]) override {
      giz1.Drag(src0, src1, dir);
    }
    void key_release(
		[[maybe_unused]] int key, 
		[[maybe_unused]] int mods) override {
    }
    void key_press(int key, [[maybe_unused]] int mods) override {
      if (key == GLFW_KEY_R) { giz1.igizmo_mode = 1; }
      if (key == GLFW_KEY_G) { giz1.igizmo_mode = 0; }
    }
    //
    void Draw() {
      { // set boundary condition
        const dfm2::CMat4<double> aff1 = giz1.Affine().cast<double>();
        for (unsigned int ip = 0; ip < aXYZ0.size() / 3; ++ip) {
          if (aBCFlag[ip * 3 + 0] == 0) { continue; }
          if (aBCFlag[ip * 3 + 0] == 1) {
            aXYZ1[ip * 3 + 0] = aXYZ0[ip * 3 + 0];
            aXYZ1[ip * 3 + 1] = aXYZ0[ip * 3 + 1];
            aXYZ1[ip * 3 + 2] = aXYZ0[ip * 3 + 2];
          }
          if (aBCFlag[ip * 3 + 0] == 2) {
            dfm2::CVec3d t0 = dfm2::Vec3_Mat4Vec3_Homography(
                aff1.mat,
                aXYZ0.data() + ip * 3);
            t0.CopyTo(aXYZ1.data() + ip * 3);
          }
        }
        for (int itr = 0; itr < 2; ++itr) {
          dfm2::UpdateRotationsByMatchingCluster_Linear(
              aQuat1,
              aXYZ0, aXYZ1,
              def0.psup_ind, def0.psup);
        }
      }
      def0.Deform(
          aXYZ1, aQuat1,
          aXYZ0, aBCFlag);
      dfm2::UpdateQuaternions_Svd(
          aQuat1,
          aXYZ0, aXYZ1,
          def0.psup_ind, def0.psup);
      // -------------------------------
      DrawBegin_oldGL();
      delfem2::opengl::DrawAxis(1);
      { // mesh
        ::glEnable(GL_LIGHTING);
        ::glColor3d(0, 0, 0);
        delfem2::opengl::DrawMeshTri3D_Edge(
            aXYZ1.data(), aXYZ1.size() / 3,
            aTri.data(), aTri.size() / 3);
        delfem2::opengl::DrawMeshTri3D_FaceNorm(
            aXYZ1.data(),
            aTri.data(), aTri.size() / 3);
      }
      { // draw bc
        ::glDisable(GL_LIGHTING);
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
      /*
      for(int ibc=0;ibc<aCntBC.size();++ibc){
        dfm2::CVec3d p0 = aCntBC[ibc];
        dfm2::opengl::DrawSphereAt(32, 32, 0.1, p0.x(), p0.y(), p0.z());
        
      }
       */
      delfem2::opengl::Draw(giz1);
      SwapBuffers();
    }
   public:
    delfem2::CGizmo_Affine<float> giz1; // bcflag==2
    std::vector<unsigned int> aTri;
    std::vector<double> aXYZ0, aXYZ1;
    std::vector<double> aQuat1;
    std::vector<int> aBCFlag;
    delfem2::Deformer_Arap def0;
  } viewer;
  // --------------------
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  // --------------------
  while (!glfwWindowShouldClose(viewer.window)) {
    viewer.Draw();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


