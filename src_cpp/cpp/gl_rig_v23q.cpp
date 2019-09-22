/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <set>
#include <map>
#include <cassert>

#if defined(__APPLE__) // Mac
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#elif defined(_WIN32) // windows
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "delfem2/mat3.h"
#include "delfem2/vec3.h"
#include "delfem2/quat.h"
#include "delfem2/funcs.h" // isFileExists
#include "delfem2/v23m3q.h"
#include "delfem2/msh.h"

#include "delfem2/rig_v3q.h"

#include "delfem2/gl2_funcs.h"
#include "delfem2/gl_tex.h"
#include "delfem2/gl_color.h"
#include "delfem2/gl_v23.h"
#include "delfem2/gl_rig_v23q.h"


#ifndef M_PI 
#define M_PI 3.1415926535
#endif

///////////////////////////////////////////////////////////////////////////////

void Draw_RigBone
(int ibone,
 bool is_selected,
 int ielem_selected,
 const std::vector<CRigBone>& aBone,
 double rad_bone_sphere,
 double rad_rot_hndlr)
{
  { // draw point
    if(is_selected){ ::glColor3d(0,1,1); }
    else{            ::glColor3d(1,0,0); }
    const CVector3 pos = aBone[ibone].Pos();
    DrawSphereAt(32, 32, rad_bone_sphere, pos.x,pos.y,pos.z);
  }
  if(is_selected){
    DrawHandlerRotation_Mat4(aBone[ibone].Mat, rad_rot_hndlr, ielem_selected);
    int ibone_parent = aBone[ibone].ibone_parent;
    if( ibone_parent>=0&&ibone_parent<(int)aBone.size() ){
      const CVector3 pp(aBone[ibone_parent].Pos());
    }
    else{
    }
  }
}

void DrawBone
(const std::vector<CRigBone>& aBone,
 int ibone_selected,
 int ielem_selected,
 double rad_bone_sphere,
 double rad_rot_hndlr)
{
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_2D);
  ::glPointSize(3);
  for( int iskel=0;iskel<(int)aBone.size();++iskel){
    Draw_RigBone(iskel,
                 (iskel==ibone_selected),ielem_selected,aBone,
                 rad_bone_sphere,rad_rot_hndlr);
  }
  // draw edges whilte
  for( int ibone=0;ibone<(int)aBone.size();++ibone){
    const CRigBone& bone = aBone[ibone];
    const int ibone_p = aBone[ibone].ibone_parent;
    if( ibone_p < 0 || ibone_p >= (int)aBone.size() ){ continue; }
    const CRigBone& bone_p = aBone[ibone_p];
    bool is_selected_p = (ibone_p == ibone_selected);
    if(is_selected_p){ ::glColor3d(1.0,1.0,1.0); }
    else{              ::glColor3d(0.0,0.0,0.0); }
    ::glBegin(GL_LINES);
    myGlVertex(bone.Pos());
    myGlVertex(bone_p.Pos());
    ::glEnd();
  }
}

