/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_GLFW_VIEWRIG_H
#define DFM2_GLFW_VIEWRIG_H

#include "delfem2/glfw/viewer3.h"
#include "delfem2/opengl/old/rigv3.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/gizmo_geo3.h"

namespace delfem2 {

/**
 * @brief gizmo for rig
 */
template<typename REAL>
class CGizmo_Rig {
public:
  enum MODE_EDIT {
    ROT, TRNSL
  };
public:
  void Pick(const float src[3], const float dir[3],
            const std::vector<CRigBone> &aBone) {
    if (mode_edit == ROT) {
      const CVec3d s0(src), d0(dir);
      if (ipicked_bone != -1) {
        gizmo_rot.Pick(true, src, dir, 0.03);
        if (gizmo_rot.ielem_picked != -1) {
          return;
        }
      }
      ipicked_bone = -1;
      for (unsigned int ib = 0; ib < aBone.size(); ++ib) {
        CVec3d p0 = aBone[ib].Pos();
        CVec3d p1 = nearest_Line_Point(p0, s0, d0);
        double len01 = (p0 - p1).norm();
        if (len01 < 0.03) {
          ipicked_bone = ib;
          gizmo_rot.quat[0] = aBone[ipicked_bone].quatRelativeRot[0];
          gizmo_rot.quat[1] = aBone[ipicked_bone].quatRelativeRot[1];
          gizmo_rot.quat[2] = aBone[ipicked_bone].quatRelativeRot[2];
          gizmo_rot.quat[3] = aBone[ipicked_bone].quatRelativeRot[3];
        }
      }
    } else if (mode_edit == TRNSL) {
      gizmo_trnsl.pos[0] = aBone[0].transRelative[0];
      gizmo_trnsl.pos[1] = aBone[0].transRelative[1];
      gizmo_trnsl.pos[2] = aBone[0].transRelative[2];
      gizmo_trnsl.Pick(true, src, dir, 0.03);
    }
  }

  bool Drag(std::vector<CRigBone> &aBone,
            const float src0[3], const float src1[3], const float dir[3]) {
    if (mode_edit == ROT) {
      if (ipicked_bone == -1) { return false; }
      assert(ipicked_bone >= 0 && ipicked_bone < (int) aBone.size());
      const delfem2::CVec3d hoge = aBone[ipicked_bone].Pos();
      gizmo_rot.pos = hoge.cast<float>();
      gizmo_rot.Drag(src0, src1, dir);
      {
        const int ibp = aBone[ipicked_bone].ibone_parent;
        CMat3d m3;
        if (ibp == -1) { m3.SetIdentity(); }
        else { m3.SetMat4(aBone[ibp].affmat3Global); }
        CQuatd qp;
        m3.GetQuat_RotMatrix(qp.p);
        CQuatd qg = CQuatf(gizmo_rot.quat).cast<double>();
        CQuatd qj = qp.conjugate() * qg;
        qj.CopyTo(aBone[ipicked_bone].quatRelativeRot);
      }
      UpdateBoneRotTrans(aBone);
      return true;
    } else if (mode_edit == TRNSL) {
      gizmo_trnsl.Drag(src0, src1, dir);
      aBone[0].transRelative[0] = gizmo_trnsl.pos.x;
      aBone[0].transRelative[1] = gizmo_trnsl.pos.y;
      aBone[0].transRelative[2] = gizmo_trnsl.pos.z;
      UpdateBoneRotTrans(aBone);
      return true;
    }
    return false;
  }

  // --------
  void SetSize(double size) {
    gizmo_rot.size = size;
    gizmo_trnsl.size = size;
  }

  void SetMode(MODE_EDIT mode) { mode_edit = mode; }

public:
  MODE_EDIT mode_edit = ROT; // 0:rot, 1:trnsl
  int ipicked_bone = -1;
  CGizmo_Transl<REAL> gizmo_trnsl;
  CGizmo_Rotation<REAL> gizmo_rot;
};

namespace glfw {

DFM2_INLINE void Draw(
    CGizmo_Rig<float>& giz,
    const std::vector<CRigBone>& aBone)
{
  if( giz.mode_edit == CGizmo_Rig<float>::MODE_EDIT::TRNSL ){ // translation
    giz.gizmo_trnsl.pos[0] = aBone[0].transRelative[0];
    giz.gizmo_trnsl.pos[1] = aBone[0].transRelative[1];
    giz.gizmo_trnsl.pos[2] = aBone[0].transRelative[2];
    opengl::Draw(giz.gizmo_trnsl);
  }
  else if( giz.mode_edit == CGizmo_Rig<float>::MODE_EDIT::ROT ){ // translation
    if( giz.ipicked_bone != -1 ){
      assert( giz.ipicked_bone >= 0 && giz.ipicked_bone < (int)aBone.size() );
      giz.gizmo_rot.pos = aBone[giz.ipicked_bone].Pos().cast<float>();
      { // set quaternion
        CMat3<double> m3;
        m3.SetMat4(aBone[giz.ipicked_bone].affmat3Global);
        CQuat<double> qj;
        m3.GetQuat_RotMatrix(qj.p);
        qj.CopyTo(giz.gizmo_rot.quat);
      }
      opengl::Draw(giz.gizmo_rot);
    }
  }
}

class CViewerGlfw_RiggedMesh : public CViewer3 {
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
    ::delfem2::opengl::DrawMeshTri3D_FaceNorm(
        aXYZ1.data(),
        aTri.data(), aTri.size() / 3);
    ::glDisable(GL_DEPTH_TEST);
    ::delfem2::opengl::DrawBone_Line(
        aBone,
        -1, -1,
        0.01, 1.0);
    ::delfem2::glfw::Draw(gizmo, aBone);
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
