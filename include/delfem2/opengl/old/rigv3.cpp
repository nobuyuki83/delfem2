/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/rigv3.h"
#include "delfem2/rig_geo3.h"

#if defined(_WIN32) // windows
  #define NOMINMAX   // to remove min,max macro
  #include <windows.h>
#endif

#if defined(__APPLE__) && defined(__MACH__) // mac
  #define GL_SILENCE_DEPRECATION
  #include <OpenGL/gl.h>
#else
  #include <GL/gl.h>
#endif

#ifndef M_PI
  #define M_PI 3.1415926535
#endif

// -------------------------------------------------------

DFM2_INLINE void delfem2::opengl::Draw_RigBone(
    int ibone,
    bool is_selected,
    [[maybe_unused]] int ielem_selected,
    const std::vector<CRigBone>& aBone,
    double rad_bone_sphere,
    [[maybe_unused]] double rad_rot_hndlr) {
  { // draw point
    if (is_selected) { ::glColor3d(0, 1, 1); }
    else { ::glColor3d(1, 0, 0); }
    const CVec3d pos = aBone[ibone].Pos();
    delfem2::opengl::DrawSphereAt(32, 32, rad_bone_sphere, pos.x, pos.y, pos.z);
  }
  /*
  if(is_selected){
    opengl::DrawHandlerRotation_Mat4(aBone[ibone].affmat3Global, rad_rot_hndlr, ielem_selected);
    int ibone_parent = aBone[ibone].ibone_parent;
    if( ibone_parent>=0&&ibone_parent<(int)aBone.size() ){
      const CVec3d pp(aBone[ibone_parent].Pos());
    }
    else{
    }
  }
   */
}

DFM2_INLINE void delfem2::opengl::DrawBone_Line(
    const std::vector<CRigBone>& aBone,
    int ibone_selected,
    int ielem_selected,
    double rad_bone_sphere,
    double rad_rot_hndlr)
{
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_2D);
  for(unsigned int iskel=0;iskel<aBone.size();++iskel){
    const bool is_selected = (int)iskel==ibone_selected;
    Draw_RigBone(iskel,
        is_selected,ielem_selected,aBone,
        rad_bone_sphere,rad_rot_hndlr);
  }
  for(unsigned int ibone=0;ibone<aBone.size();++ibone){
    const CRigBone& bone = aBone[ibone];
    const int ibone_p = aBone[ibone].ibone_parent;
    if( ibone_p < 0 || ibone_p >= (int)aBone.size() ){ continue; }
    const CRigBone& bone_p = aBone[ibone_p];
    bool is_selected_p = (ibone_p == ibone_selected);
    if(is_selected_p){ ::glColor3d(1.0,1.0,1.0); } // white if selected
    else{              ::glColor3d(0.0,0.0,0.0); } // black if not selected
    ::glBegin(GL_LINES);
    opengl::myGlVertex(bone.Pos());
    opengl::myGlVertex(bone_p.Pos());
    ::glEnd();
  }
}

DFM2_INLINE void delfem2::opengl::DrawBone_Octahedron(
    const std::vector<CRigBone>& aBone,
    unsigned int ibone_selected,
    [[maybe_unused]] unsigned int ielem_selected,
    [[maybe_unused]] double rad_bone_sphere,
    [[maybe_unused]] double rad_rot_hndlr)
{
  namespace dfm2 = delfem2;
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_2D);
  /*
  for(unsigned int iskel=0;iskel<aBone.size();++iskel){
    const bool is_selected = (iskel==ibone_selected);
    Draw_RigBone(
        iskel,
        is_selected,ielem_selected,aBone,
        rad_bone_sphere,rad_rot_hndlr);
  }
   */
  for(unsigned int ibone=0;ibone<aBone.size();++ibone){
    const CRigBone& bone = aBone[ibone];
    const int ibone_p = aBone[ibone].ibone_parent;
    if( ibone_p < 0 || ibone_p >= (int)aBone.size() ){ continue; }
    const CRigBone& bone_p = aBone[ibone_p];
    bool is_selected_p = (ibone_p == (int)ibone_selected);
    const CVec3d p0(bone_p.invBindMat[3], bone_p.invBindMat[7], bone_p.invBindMat[11]);
    const CVec3d p1(bone.invBindMat[3], bone.invBindMat[7], bone.invBindMat[11]);
    ::glPushMatrix();
    double At[16]; dfm2::Transpose_Mat4(At,bone_p.affmat3Global);
    ::glMultMatrixd(At);
    ::glEnable(GL_LIGHTING);
    delfem2::opengl::DrawArrowOcta_FaceNrm(
        dfm2::CVec3d(0,0,0),
        p0-p1,
        0.1,
        0.2);
    //
    ::glDisable(GL_LIGHTING);
    if(is_selected_p){ ::glColor3d(1.0,1.0,1.0); } // white if selected
    else{              ::glColor3d(0.0,0.0,0.0); } // black if not selected
    delfem2::opengl::DrawArrowOcta_Edge(
        dfm2::CVec3d(0,0,0),
        p0-p1,
        0.1,
        0.2);
    /*
    delfem2::opengl::DrawArrow(
        dfm2::CVec3d(0,0,0),
        p0-p1);
        */
    ::glPopMatrix();
  }
}

DFM2_INLINE void delfem2::opengl::DrawJoints(
    const std::vector<double>& aJntPos,
    const std::vector<int>& aIndBoneParent)
{
  for(unsigned int ib=0;ib<aJntPos.size()/3;++ib){
    const double* p = aJntPos.data()+ib*3;
    ::glColor3d(0,0,1);
    ::glDisable(GL_LIGHTING);
    ::glDisable(GL_DEPTH_TEST);
    opengl::DrawSphereAt(8, 8, 0.01, p[0], p[1], p[2]);
    int ibp = aIndBoneParent[ib];
    if( ibp == -1 ){ continue; }
    const double* pp = aJntPos.data()+ibp*3;
    ::glColor3d(0,0,0);
    ::glBegin(GL_LINES);
    ::glVertex3dv(p);
    ::glVertex3dv(pp);
    ::glEnd();
  }
}

