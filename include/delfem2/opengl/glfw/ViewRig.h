/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_OPENGL_GLFW_VIEWRIG_H
#define DFM2_OPENGL_GLFW_VIEWRIG_H

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/rigv3_glold.h"
#include "delfem2/gizmo_geo3.h"

namespace delfem2 {
namespace opengl {

class CViewerGlfw_RiggedMesh : public CViewer_GLFW {
public:
  CViewerGlfw_RiggedMesh() {
    gizmo.SetSize(0.3);
  }

  void SetRiggedMesh(
      std::vector<double> &aXYZ0_,
      std::vector<unsigned int> &aTri_,
      std::vector<double> &aSkinningSparseW_,
      std::vector<unsigned int> &aSkinningSparseI_,
      std::vector<CRigBone> &aBone_) {
    aXYZ0 = aXYZ0_;
    aTri = aTri_;
    aSkinningSparseW = aSkinningSparseW_;
    aSkinningSparseI = aSkinningSparseI_;
    aBone = aBone_;
    aXYZ1 = aXYZ0;
  }

  void SetUndeformed() {
    for (auto &ib : aBone) {
      Quat_Identity(ib.quatRelativeRot);
    }
    aBone[0].transRelative[0] = 0.0;
    aBone[0].transRelative[1] = 0.0;
    aBone[0].transRelative[2] = 0.0;
    UpdateBoneRotTrans(aBone);
    SkinningSparse_LBS(
        aXYZ1,
        aXYZ0, aBone, aSkinningSparseW, aSkinningSparseI);
  }

  void Draw() {
    ::glEnable(GL_LIGHTING);
    ::glEnable(GL_DEPTH_TEST);
    DrawMeshTri3D_FaceNorm(
        aXYZ1.data(),
        aTri.data(), aTri.size() / 3);
    ::glDisable(GL_DEPTH_TEST);
    DrawBone(
        aBone,
        -1, -1,
        0.01, 1.0);
    ::delfem2::opengl::Draw(gizmo, aBone);
  }

  void mouse_drag(
      const float src0[3],
      const float src1[3],
      const float dir[3]) override
  {
    const bool is_edited = gizmo.Drag(
        aBone,
        src0, src1, dir);
    if (!is_edited) { return; }
    SkinningSparse_LBS(
        aXYZ1,
        aXYZ0, aBone, aSkinningSparseW, aSkinningSparseI);
  }

  void mouse_press(const float src[3], const float dir[3]) override {
    gizmo.Pick(src, dir, aBone);
    std::cout << "ib:" << gizmo.ipicked_bone << " " << aBone.size() << std::endl;
  }

  void key_press(int key, int mods) override {
    if (key == GLFW_KEY_G) { gizmo.SetMode(CGizmo_Rig<float>::MODE_EDIT::TRNSL); }
    if (key == GLFW_KEY_R) { gizmo.SetMode(CGizmo_Rig<float>::MODE_EDIT::ROT); }
    if (key == GLFW_KEY_S) {
      std::ofstream fout("pose.txt");
      for (const auto &bone: aBone) {
        const double *q = bone.quatRelativeRot;
        fout << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << std::endl;
      }
      {
        const double *t = aBone[0].transRelative;
        fout << t[0] << " " << t[1] << " " << t[2] << std::endl;
      }
    }
  }

public:
  std::vector<double> aXYZ0, aXYZ1;
  std::vector<unsigned int> aTri;
  std::vector<double> aSkinningSparseW;
  std::vector<unsigned int> aSkinningSparseI;
  std::vector<CRigBone> aBone;
  CGizmo_Rig<float> gizmo;
};

}
}


#endif