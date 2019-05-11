/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <set>
#include <map>
#include <cassert>

#if defined(__APPLE__) && defined(__MACH__) // Mac
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/glu.h>
#elif defined(_WIN32) // windows
#include <windows.h>
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#else
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "delfem2/mat3.h"
#include "delfem2/vec3.h"
#include "delfem2/quat.h"
#include "delfem2/funcs.h" // isFileExists
#include "delfem2/v23m3q.h"
#include "delfem2/msh.h"
#include "delfem2/tex.h"

#include "delfem2/funcs_gl.h"
#include "delfem2/color_gl.h"
#include "delfem2/v23q_gl.h"

#include "delfem2/rigmesh.h"

#ifndef M_PI 
#define M_PI 3.1415926535
#endif

static double myStod(const std::string& str){
  char* e;
  double d = std::strtod(str.c_str(),&e);
  return d;
}

// probably std::stroi is safer to use but it is only for C++11
static int myStoi(const std::string& str){
  char* e;
  long d = std::strtol(str.c_str(),&e,0);
  return (int)d;
}

bool isActive_AABB(const double aabb[6]){
  if( aabb[0]>aabb[1] ) return false;
  return true;
}

void copy_AABB(double aabb[6], const double aabb0[6]){
  if( aabb == aabb0 ) return;
  for(int i=0;i<6;++i){ aabb[i] = aabb0[i]; }
}

void myAdd_AABB(double aabb[6], const double aabb0[6], const double aabb1[6])
{
  if( !isActive_AABB(aabb0) && !isActive_AABB(aabb1) ){ aabb[0]=1; aabb[1]=-1; return; }
  if( !isActive_AABB(aabb0) ){ copy_AABB(aabb, aabb1); return; }
  if( !isActive_AABB(aabb1) ){ copy_AABB(aabb, aabb0); return; }
  aabb[0] = ( aabb0[0] < aabb1[0] ) ? aabb0[0] : aabb1[0];
  aabb[1] = ( aabb0[1] > aabb1[1] ) ? aabb0[1] : aabb1[1];
  aabb[2] = ( aabb0[2] < aabb1[2] ) ? aabb0[2] : aabb1[2];
  aabb[3] = ( aabb0[3] > aabb1[3] ) ? aabb0[3] : aabb1[3];
  aabb[4] = ( aabb0[4] < aabb1[4] ) ? aabb0[4] : aabb1[4];
  aabb[5] = ( aabb0[5] > aabb1[5] ) ? aabb0[5] : aabb1[5];
}

///////////////////////////////////////////////////////////////////////////////

void CBone_RigMsh::Draw
(bool is_selected,
 int ielem_selected,
 const std::vector<CBone_RigMsh>& aBone,
 double bone_rad) const
{
  if( !is_active ) return;
  { // draw point
    if(is_selected){ ::glColor3d(0,1,1); }
    else{            ::glColor3d(1,0,0); }
    DrawSphereAt(32, 32, bone_rad, pos[0],pos[1],pos[2]);
  }
  if(is_selected){
    const CVector3 p0(pos);
    DrawHandlerRotation(p0, quat, 15, ielem_selected);
    if( ibone_parent>=0&&ibone_parent<(int)aBone.size() && aBone[ibone_parent].is_active ){
      const CVector3 pp(aBone[ibone_parent].pos);
//      double len = (p0-pp).Length();
//      drawCircle((p0-pp).Normalize(), p0, len*1.5);
    }
    else{
    }
  }
}

int CBone_RigMsh::PickHandler
(const CVector3& org,
 const CVector3& dir,
 double rad_handlr,
 double tol) const
{
  const CVector3 p0(pos);
  return PickHandlerRotation(org,dir,
                             p0, quat, rad_handlr,
                             tol);
}




void CBone_RigMsh::Affine(const double A[16])
{
  Affine3D(pos_ini, A, pos_ini);
  Affine3D(pos, A, pos);
}

/*
 void Scale(double scale){
 pos_ini[0] *= scale;
 pos_ini[1] *= scale;
 pos_ini[2] *= scale;
 pos[0] *= scale;
 pos[1] *= scale;
 pos[2] *= scale;
 }
 void Translate(double dx, double dy, double dz){
 pos_ini[0] += dx;
 pos_ini[1] += dy;
 pos_ini[2] += dz;
 pos[0] += dx;
 pos[1] += dy;
 pos[2] += dz;
 }
 */

////////////////////////////////////////////////////////////////////////////
// from here std::vector<CBone_RigMsh>

void DrawBone
(const std::vector<CBone_RigMsh>& aBone,
 int ibone_selected,
 int ielem_selected,
 double bone_rad)
{
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_2D);
  ::glPointSize(3);
  for( int iskel=0;iskel<(int)aBone.size();++iskel){
    const CBone_RigMsh& bone = aBone[iskel];
    bone.Draw((iskel==ibone_selected),ielem_selected,aBone,bone_rad);
  }
  // draw edges whilte
  for( int ibone=0;ibone<(int)aBone.size();++ibone){
    const CBone_RigMsh& bone = aBone[ibone];
    if( !bone.is_active ) continue;
    const int ibone_p = aBone[ibone].ibone_parent;
    if( ibone_p < 0 || ibone_p >= (int)aBone.size() ){ continue; }
    const CBone_RigMsh& bone_p = aBone[ibone_p];
    if( !bone_p.is_active ) continue;
    bool is_selected_p = (ibone_p == ibone_selected);
    if(is_selected_p){ ::glColor3d(1.0,1.0,1.0); }
    else{              ::glColor3d(0.0,0.0,0.0); }
    ::glBegin(GL_LINES);
    ::glVertex3dv(bone.pos);
    ::glVertex3dv(bone_p.pos);
    ::glEnd();
  }
}


void ReadBVH
(std::vector<CBone_RigMsh>& aBone,
 std::vector<CChannel_RotTransBone_BVH>& aChannelRotTransBone,
 int& nframe,
 std::vector<double>& aValueRotTransBone,
 const std::string& path_bvh)
{
  std::ifstream fin;
  fin.open(path_bvh);
  if( !fin.is_open() ){
    std::cout << "cannot open file" << std::endl;
    return;
  }
  aBone.clear();
  aChannelRotTransBone.clear();
  ///////////////////////////////////////////
  std::string line;
  std::vector<int> stackIndBone;
  while(std::getline(fin,line)){
    if (line[line.size()-1] == '\n') line.erase(line.size()-1); // remove the newline code
    if (line[line.size()-1] == '\r') line.erase(line.size()-1); // remove the newline code
    line = Replace(line, '\t', ' ');
    std::vector<std::string> aToken = Split(line,' ');
    std::cout << aToken[0] << std::endl;
    if( aToken[0] == "HIERARCHY" ){
      assert(aBone.empty());
    }
    else if( aToken[0] == "ROOT" ){
      assert(aBone.size()==0);
      CBone_RigMsh br;
      assert( aToken.size() == 2 );
      br.name = aToken[1];
      aBone.push_back(br);
    }
    else if( aToken[0] == "{" ){
      stackIndBone.push_back(aBone.size()-1);
      if( stackIndBone.size() > 1 ){
        int ibp = stackIndBone[stackIndBone.size()-2];
        int ib = aBone.size()-1;
        aBone[ib].ibone_parent  = ibp;
      }
    }
    else if( aToken[0] == "}" ){
      stackIndBone.resize(stackIndBone.size()-1);
    }
    else if( aToken[0] == "OFFSET"){
      assert( aToken.size()==4 );
      int ib = aBone.size()-1;
      aBone[ib].pos_ini[0] = myStod(aToken[1]);
      aBone[ib].pos_ini[1] = myStod(aToken[2]);
      aBone[ib].pos_ini[2] = myStod(aToken[3]);
      if( stackIndBone.size() > 1 ){
        const int ibp = stackIndBone[stackIndBone.size()-2];
        assert(ibp<(int)aBone.size());
        aBone[ib].pos_ini[0] += aBone[ibp].pos_ini[0];
        aBone[ib].pos_ini[1] += aBone[ibp].pos_ini[1];
        aBone[ib].pos_ini[2] += aBone[ibp].pos_ini[2];
      }
    }
    else if( aToken[0] == "CHANNELS" ){
      assert(aToken.size()>=2);
      int nch = myStoi(aToken[1]);
      assert((int)aToken.size()==nch+2);
      const int ib = aBone.size()-1;
      for(int ich=0;ich<nch;++ich){
        const std::string& type_ch = aToken[ich+2];
        if(      type_ch == "Xposition" ){ aChannelRotTransBone.push_back( CChannel_RotTransBone_BVH(ib,0,false) ); }
        else if( type_ch == "Yposition" ){ aChannelRotTransBone.push_back( CChannel_RotTransBone_BVH(ib,1,false) ); }
        else if( type_ch == "Zposition" ){ aChannelRotTransBone.push_back( CChannel_RotTransBone_BVH(ib,2,false) ); }
        else if( type_ch == "Xrotation" ){ aChannelRotTransBone.push_back( CChannel_RotTransBone_BVH(ib,0,true) ); }
        else if( type_ch == "Yrotation" ){ aChannelRotTransBone.push_back( CChannel_RotTransBone_BVH(ib,1,true) ); }
        else if( type_ch == "Zrotation" ){ aChannelRotTransBone.push_back( CChannel_RotTransBone_BVH(ib,2,true) ); }
        else{
          std::cout << "ERROR-->undefiend type" << std::endl;
        }
      }
    }
    else if( aToken[0] == "JOINT" ){
      CBone_RigMsh br;
      assert( aToken.size() == 2 );
      br.name = aToken[1];
      aBone.push_back(br);
    }
    else if( aToken[0] == "End" ){
      assert(aToken[1] == "Site");
      CBone_RigMsh br;
      assert( aToken.size() == 2 );
      br.name = aToken[1];
      aBone.push_back(br);
    }
    else if( aToken[0] == "MOTION"){
      break;
    }
  }
  nframe = 0;
  {
    std::string stmp0;
    {
      std::getline(fin,line);
      std::stringstream ss(line);
      ss >> stmp0 >> nframe;
      std::cout << "frame: " << nframe << std::endl;
    }
    std::getline(fin,line);
    std::cout << "frametime: " << line << std::endl;
  }
  const int nchannel = aChannelRotTransBone.size();
  aValueRotTransBone.resize(nframe*nchannel);
  for(int iframe=0;iframe<nframe;++iframe){
    std::getline(fin,line);
    line = Replace(line, '\t', ' ');
    if (line[line.size()-1] == '\n') line.erase(line.size()-1); // remove the newline code
    if (line[line.size()-1] == '\r') line.erase(line.size()-1); // remove the newline code
    std::vector<std::string> aToken = Split(line,' ');
    std::cout << aToken.size() << " " << aChannelRotTransBone.size() << std::endl;
    assert(aToken.size()==aChannelRotTransBone.size());
    for(int ich=0;ich<nchannel;++ich){
      aValueRotTransBone[iframe*nchannel+ich] = myStod(aToken[ich]);
    }
  }
  ///////
  InitializeBone(aBone);
}

void InitializeBone
(std::vector<CBone_RigMsh>& aBone)
{
  for(unsigned int ib=0;ib<aBone.size();++ib){
    CBone_RigMsh& bone = aBone[ib];
    bone.pos[0] = bone.pos_ini[0];
    bone.pos[1] = bone.pos_ini[1];
    bone.pos[2] = bone.pos_ini[2];
    bone.quat[0] = 1.0;
    bone.quat[1] = 0.0;
    bone.quat[2] = 0.0;
    bone.quat[3] = 0.0;
    bone.quat_joint[0] = 1.0;
    bone.quat_joint[1] = 0.0;
    bone.quat_joint[2] = 0.0;
    bone.quat_joint[3] = 0.0;
  }
}

void UpdateBoneRotTrans
(std::vector<CBone_RigMsh>& aBone)
{
  for(unsigned int ibone=0;ibone<aBone.size();++ibone){
    const int ibone_p = aBone[ibone].ibone_parent;
    if( ibone_p < 0 || ibone_p >= (int)aBone.size() ){ // root bone
      QuatCopy(aBone[ibone].quat,aBone[ibone].quat_joint);
      continue;
    }
    const double* quat_parent = aBone[ibone_p].quat;
    CVector3 v0 = CVector3(aBone[ibone].pos_ini) - CVector3(aBone[ibone_p].pos_ini);
    CVector3 v1 = QuatVec(quat_parent, v0) + CVector3(aBone[ibone_p].pos);
    v1.CopyValueTo(aBone[ibone].pos);
    QuatMult(aBone[ibone].quat, aBone[ibone].quat_joint, quat_parent);
  }
}

void SetRotTransBVH
(std::vector<CBone_RigMsh>& aBone,
 const std::vector<CChannel_RotTransBone_BVH>& aChannelRotTransBone,
 const double *aVal)
{
  InitializeBone(aBone);
  const int nch = aChannelRotTransBone.size();
  for(int ich=0;ich<nch;++ich){
    const int ibone = aChannelRotTransBone[ich].ibone;
    const int iaxis = aChannelRotTransBone[ich].iaxis;
    const bool isrot = aChannelRotTransBone[ich].isrot;
    const double val = aVal[ich];
    assert(ibone<(int)aBone.size());
    assert(iaxis>=0&&iaxis<3);
    if( !isrot ){
      aBone[ibone].pos[iaxis] = val;
    }
    else{
      const double ar = -val*M_PI/180.0;
      double v0[3] = {0,0,0};
      v0[iaxis] = 1.0;
      double dq[4] = { cos(ar*0.5), v0[0]*sin(ar*0.5), v0[1]*sin(ar*0.5), v0[2]*sin(ar*0.5) };
      double qtmp[4]; QuatMult(qtmp, dq, aBone[ibone].quat_joint);
      QuatCopy(aBone[ibone].quat_joint,qtmp);
    }
  }
  UpdateBoneRotTrans(aBone);
}

void PickBone
(int& ibone_selected,
 int& ielem_selected,
 const std::vector<CBone_RigMsh>& aBone,
 const CVector3& src,
 const CVector3& dir,
 double rad_hndlr,
 double tol)
{
  if( ibone_selected>=0 && ibone_selected<(int)aBone.size() ){
    const CBone_RigMsh& bone = aBone[ibone_selected];
    ielem_selected = bone.PickHandler(src,dir,rad_hndlr,tol);
  }
  if( ielem_selected == -1 ){
    ibone_selected = -1;
    for(int ibone=0;ibone<(int)aBone.size();++ibone){
      CVector3 pos(aBone[ibone].pos);
      double distance = Distance(nearest_Line_Point(pos, src, dir),pos);
      if( distance < tol ){
        ibone_selected = ibone;
        break;
      }
    }
  }
}


void CBoneGoal::GetGoalPos
(double* pos_trg,
 const double* org_rot,
 const double* pos_cur) const
{
  if( itype == 0 ){
    pos_trg[0] = pos[0];
    pos_trg[1] = pos[1];
    pos_trg[2] = pos[2];
  }
  else if( itype == 1 ){
    CVector3 s0(pos);
    CVector3 d0(dir);
    double t0,t1;
    {
      double rad = Distance3D(pos_cur,org_rot);
      double a = d0*d0;
      double b = 2*(s0-org_rot)*d0;
      double c = (s0-org_rot)*(s0-org_rot)-rad*rad;
      double t01 = -0.5*b/a;
      double det = b*b-4*a*c;
      t0=t01;
      t1=t01;
      if( det > 0 ){
        t0 -= 0.5*sqrt(det)/a;
        t1 += 0.5*sqrt(det)/a;
      }
    }
    CVector3 p0 = s0+t0*d0;
    CVector3 p1 = s0+t1*d0;
    if( (pos_cur-p0).Length()<(pos_cur-p1).Length() ){
      pos_trg[0] = p0[0];
      pos_trg[1] = p0[1];
      pos_trg[2] = p0[2];
    }
    else{
      pos_trg[0] = p1[0];
      pos_trg[1] = p1[1];
      pos_trg[2] = p1[2];
    }
  }
}


void BoneOptimization
(std::vector<CBone_RigMsh>& aBone,
 const std::vector<CBoneGoal>& aBoneGoal)
{
  if( aBoneGoal.empty() ) return;
  ///
  std::vector< std::vector<int> > src2trg(aBone.size());
  for(int itrg=0;itrg<(int)aBoneGoal.size();++itrg){
    int ib0 = aBoneGoal[itrg].ibone;
    //  src2trg[ib0].push_back(itrg); // this should be added if we set goal apart from the bone
    for(;;){
      int ib1 = aBone[ib0].ibone_parent;
      if( ib1 < 0 ){ break; }
      assert( ib1 < (int)aBoneGoal.size() );
      src2trg[ib1].push_back(itrg);
      ib0 = ib1;
    }
  }
  
  /*
   std::cout << "#########" << std::endl;
   for(int ibone=0;ibone<aBone.size();++ibone){
   std::cout << "src2trg: " << ibone << " --> ";
   for(int itrg=0;itrg<src2trg[ibone].size();++itrg){
   std::cout << src2trg[ibone][itrg] << " ";
   }
   std::cout << std::endl;
   }
   */
  for(unsigned int ibone=0;ibone<aBone.size();++ibone){
    //  for(int ibone=aBone.size()-1;ibone>=0;--ibone){
    if( src2trg[ibone].empty() ){ continue; }
    //    std::cout << "optimize bone:" << ibone << std::endl;
    {
      const CVector3 org(aBone[ibone].pos);
      CMatrix3 A;
      double len = 0.0;
      for(unsigned int iitrg=0;iitrg<src2trg[ibone].size();++iitrg){
        const int itrg = src2trg[ibone][iitrg];
        const int ib0 = aBoneGoal[itrg].ibone;
        double pos_goal[3];
        aBoneGoal[itrg].GetGoalPos(pos_goal,
                                   aBone[ibone].pos, aBone[ib0].pos);
        CVector3 a = CVector3(aBone[ib0].pos) - org;
        CVector3 b = pos_goal-org;
        A += OuterProduct(a,b);
        len += a.Length() + b.Length();
      }
      A += CMatrix3::Identity()*len*len*1.0e-3;
      CMatrix3 R; GetRotPolarDecomp(R.mat, A.mat, 20);
      double quat[4]; R.GetQuat_RotMatrix(quat);
      double quat0[4]; QuatMult(quat0, aBone[ibone].quat_joint,quat);
      QuatCopy(aBone[ibone].quat_joint,quat0);
      UpdateBoneRotTrans(aBone);
    }
    {
      double off[3] = {0,0,0};
      double w = 0;
      for(unsigned int itrg=0;itrg<aBoneGoal.size();++itrg){
        const int ib0 = aBoneGoal[itrg].ibone;
        double pos_goal[3];
        aBoneGoal[itrg].GetGoalPos(pos_goal,
                                   aBone[ibone].pos, aBone[ib0].pos);
        off[0] += aBone[ib0].pos[0] - pos_goal[0];
        off[1] += aBone[ib0].pos[1] - pos_goal[1];
        off[2] += aBone[ib0].pos[2] - pos_goal[2];
        w += 1.0;
      }
      off[0] /= w;
      off[1] /= w;
      off[2] /= w;
      for(unsigned int jbone=0;jbone<aBone.size();++jbone){
        aBone[jbone].pos[0] -= off[0];
        aBone[jbone].pos[1] -= off[1];
        aBone[jbone].pos[2] -= off[2];
      }
    }
  }
}

void DrawBoneTarget
(const std::vector<CBoneGoal>& aBoneGoal,
 const std::vector<CBone_RigMsh>& aBone)
{
  double bone_rad = 6;
  ::glDisable(GL_LIGHTING);
  for(unsigned int itrg=0;itrg<aBoneGoal.size();++itrg){
    const CBoneGoal& bg = aBoneGoal[itrg];
    if( bg.itype == 0 ){ ::glColor3d(0,0,1); }
    if( bg.itype == 1 ){ ::glColor3d(0,0,0); }
    const double* pos = bg.pos;
    DrawSphereAt(32, 32, bone_rad, pos[0],pos[1],pos[2]);
    if( bg.itype == 1 ){
      CVector3 dir(bg.dir);
      CVector3 src(bg.pos);
      ::glBegin(GL_LINES);
      myGlVertex(src+1000*dir);
      myGlVertex(src-1000*dir);
      ::glEnd();
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

void CSkin_RigMsh::Finalize(const int npoint)
{
  // making point 2 bone indexing
  aPoint2BoneInd.assign(npoint+1,0);
  for(unsigned int ibone=0;ibone<aBone.size();++ibone){
    const CBoneWeight_RigMsh& bone = aBone[ibone];
    for(unsigned int iip=0;iip<bone.aIpWeight.size();++iip){
      int ip0 = bone.aIpWeight[iip].first;
      assert( ip0 >= 0 && ip0 < npoint );
      aPoint2BoneInd[ip0+1] += 1;
    }
  }
  for(int ip=0;ip<npoint;++ip){
    aPoint2BoneInd[ip+1] += aPoint2BoneInd[ip];
  }
  const int np2bi = aPoint2BoneInd[npoint];
  aPoint2Bone.resize(np2bi);
  aPoint2BoneWeight.resize(np2bi);
  for(unsigned int ibone=0;ibone<aBone.size();++ibone){
    const CBoneWeight_RigMsh& bone = aBone[ibone];
    for(unsigned int iip=0;iip<bone.aIpWeight.size();++iip){
      int ip0 = bone.aIpWeight[iip].first;
      int ip2bi = aPoint2BoneInd[ip0];
      aPoint2BoneWeight[ip2bi] = bone.aIpWeight[iip].second;
      aPoint2Bone[ip2bi] = ibone;
      aPoint2BoneInd[ip0] += 1;
    }
  }
  for(int ipoint=npoint-1;ipoint>=0;--ipoint){
    aPoint2BoneInd[ipoint+1] = aPoint2BoneInd[ipoint];
  }
  aPoint2BoneInd[0] = 0;
  ////
  // weight normalization
  assert( aPoint2BoneInd[npoint] == np2bi );
  for(int ipoint=0;ipoint<npoint;++ipoint){
    double weight = 0.0;
    for(int ipobo=aPoint2BoneInd[ipoint];ipobo<aPoint2BoneInd[ipoint+1];++ipobo){
      weight += aPoint2BoneWeight[ipobo];
    }
    double weight_inv = 1.0/weight;
    for(int ipobo=aPoint2BoneInd[ipoint];ipobo<aPoint2BoneInd[ipoint+1];++ipobo){
      aPoint2BoneWeight[ipobo] *= weight_inv;
    }
  }
}

/*
void CSkin_RigMsh::SetBoneHierarchy(const std::vector< std::pair<std::string,std::string> >& aSkeletonName)
{ // put hierachy in the
  std::map<std::string,int> mapName2Indb;
  for(int ibone=0;ibone<aBone.size();++ibone){
    mapName2Indb.insert( std::make_pair(aBone[ibone].name,ibone) );
  }
  for(int iskel=0;iskel<aSkeletonName.size();++iskel){
    std::map<std::string,int>::iterator itr = mapName2Indb.find(aSkeletonName[iskel].first);
    std::map<std::string,int>::iterator itr_p = mapName2Indb.find(aSkeletonName[iskel].second);
    int ibone = (itr==mapName2Indb.end()) ? -1:itr->second;
    int ibone_p = (itr_p==mapName2Indb.end()) ? -1:itr_p->second;
    if( ibone >= 0 && ibone < aBone.size() ){
      aBone[ibone].ibone_parent = ibone_p;
    }
  }
}
 */

void CSkin_RigMsh::SetSkeleton
(std::vector<double>& aXYZ,
 const std::vector<CBone_RigMsh>& aGBone,
 const std::vector<double>& aXYZ_ini)
{
  assert(aXYZ.size()==aXYZ_ini.size());
  std::vector<int> LB2GB(aBone.size(),-1);
  for(unsigned int ilb=0;ilb<aBone.size();++ilb){
    for(unsigned int igb=0;igb<aGBone.size();++igb){
      if( aGBone[igb].name == aBone[ilb].name ){ LB2GB[ilb] = igb; break; }
    }
  }
  const int np = aXYZ.size()/3;
  assert( (int)aPoint2BoneInd.size() == np+1 );
  aXYZ.assign(np*3,0.0);
  for(int ip=0;ip<np;++ip){
    for(int iilb=aPoint2BoneInd[ip];iilb<aPoint2BoneInd[ip+1];++iilb){
      const int ilb = aPoint2Bone[iilb];
      const double weight = aPoint2BoneWeight[iilb];
      const int igb = LB2GB[ilb];
      const double* pos_b_ini = aGBone[igb].pos_ini;
      const double* pos_b = aGBone[igb].pos;
      double v0[3] = {
        aXYZ_ini[ip*3+0]-pos_b_ini[0],
        aXYZ_ini[ip*3+1]-pos_b_ini[1],
        aXYZ_ini[ip*3+2]-pos_b_ini[2]};
      double v1[3]; QuatVec(v1, aGBone[igb].quat, v0);
      aXYZ[ip*3+0] += (pos_b[0]+v1[0])*weight;
      aXYZ[ip*3+1] += (pos_b[1]+v1[1])*weight;
      aXYZ[ip*3+2] += (pos_b[2]+v1[2])*weight;
    }
  }
}

//////////////////////////////////


void CMesh_RigMsh::DrawLayerWithTex(int ilayer,const CTexManager& tex_manager, bool is_ini) const
{
  assert(aXYZ_ini.size()==aXYZ.size());
  const std::vector<double>& aXYZ0 = (is_ini)?aXYZ_ini:aXYZ;
  
  if( ilayer < 0 || ilayer >= (int)aLayer.size() ){ return; }
  const CLayer_RigMsh& layer = aLayer[ilayer];
  if(layer.material_mapping_mode == "ALL_SAME" || layer.material_mapping_mode == "" ){
    if( this->isTextureWithUVSetName(layer.uv_setname) ){
      const CTextureInfo_RigMsh& tex = this->getTextureWithUVSetName(layer.uv_setname);
      tex_manager.BindTexturePath(tex.full_path);
      ::glDisable(GL_LIGHTING);
      glColor3d(1,1,1);
      const std::vector<double>& aUV = layer.aUV;
      DrawMeshElem3D_FaceNorm(aXYZ0,aElemInd,aElem,aUV);
    }
    else{
      ::glEnable(GL_LIGHTING);
      DrawMeshElem3D_FaceNorm(aXYZ0,aElemInd,aElem);
    }
  }
  else if(layer.material_mapping_mode == "BY_FACE" ){
    assert( layer.material_mapping_face.size() == aElemInd.size()-1 );
    const int nmat = layer.material_mapping_mat2face.size();
    assert( nmat == (int)aMaterial.size() );
    for(int imat=0;imat<nmat;++imat){
      const std::vector<int>& part = layer.material_mapping_mat2face[imat];
      const CMaterial_RigMsh& mat = aMaterial[imat];
      if( mat.aTexture_Diffuse.size() == 1 ){
        const CTextureInfo_RigMsh& tex = mat.aTexture_Diffuse[0];
        tex_manager.BindTexturePath(tex.full_path);
//        std::cout << " " << tex.full_path << " " << part.size() << std::endl;
        ::glDisable(GL_LIGHTING);
        ::glEnable(GL_TEXTURE_2D);
        glColor3d(1,1,1);
        const std::vector<double>& aUV = layer.aUV;
        DrawMeshElemPart3D_FaceNorm_TexPoEl(aXYZ0,aElemInd,aElem,part,aUV);
      }
      else{
//        std::cout << "hogehoge " << mat.aTexture_Diffuse.size() << "   " << layer.uv_setname << " " << imat << " " << part.size() << std::endl;
        ::glDisable(GL_TEXTURE_2D);
        ::glEnable(GL_LIGHTING);
        std::vector<double> tmp;
        DrawMeshElemPart3D_FaceNorm_TexPoEl(aXYZ0,aElemInd,aElem,part,tmp);
      }
    }
  }
  else{ // unknown mapping mode
    ::glEnable(GL_LIGHTING);
    DrawMeshElem3D_FaceNorm(aXYZ0,aElemInd,aElem);
  }
}

void CMesh_RigMsh::Affine(const double A[16])
{
  for(int ip=0;ip<(int)aXYZ_ini.size()/3;++ip){
    double* p =  aXYZ_ini.data()+ip*3;
    Affine3D(p, A, p);
  }
  for(int ip=0;ip<(int)aXYZ.size()/3;++ip){
    double* p =  aXYZ.data()+ip*3;
    Affine3D(p, A, p);
  }
}

////////////////////////////////////////////////////////////////

void CRigMsh::Affine(const double A[16])
{
  for(unsigned int im=0;im<aMesh.size();++im){ aMesh[im].Affine(A); }
  for(unsigned int ib=0;ib<aBone.size();++ib){ aBone[ib].Affine(A); }
}

std::vector<double> CRigMsh::MinMaxXYZ() const
{
  double mm[6] = {+1,-1, 0,0, 0,0};
  for(int imesh=0;imesh<(int)this->aMesh.size();++imesh){
    const CMesh_RigMsh& mesh = aMesh[imesh];
    double mm0[6]; ::MinMaxXYZ(mm0,mesh.aXYZ_ini);
    myAdd_AABB(mm, mm, mm0);
  }
  return std::vector<double>(mm,mm+6);
}

void CRigMsh::Draw(const CTexManager& tex_manager) const{
  if( is_draw_weight ){
    std::vector< std::pair<double,CColor> > colorMap;
    colorMap.push_back( std::make_pair(0.0,CColor(color_bone_weight_back) ) );
    colorMap.push_back( std::make_pair(1.0,CColor(color_bone_weight) ) );
    for(int imesh=0;imesh<(int)this->aMesh.size();++imesh){
      const CMesh_RigMsh& mesh = aMesh[imesh];
      if( mesh.aWeight.size() != mesh.aXYZ_ini.size()/3 ){ return; }
      glEnable(GL_DEPTH_TEST);
      glDisable(GL_LIGHTING);
      DrawMeshElem3D_Scalar_Vtx(mesh.aXYZ_ini, mesh.aElemInd, mesh.aElem, mesh.aWeight.data(),colorMap);
    }
  }
  else{
    for(int imesh=0;imesh<(int)this->aMesh.size();++imesh){
//      std::cout << imesh << std::endl;
      glEnable(GL_DEPTH_TEST);
      glDisable(GL_LIGHTING);
      glDisable(GL_TEXTURE_2D);
      aMesh[imesh].DrawLayerWithTex(aMesh[imesh].ilayer_active,tex_manager,false);
//      break;
    }
  }
  
  if( is_draw_bone ){
    glDisable(GL_DEPTH_TEST);
    DrawBone(aBone,ibone_selected,ielem_selected,draw_rep_length*0.005);
  }
  glEnable(GL_DEPTH_TEST);
}

void CRigMsh::Drag(double spx, double spy, double dsx, double dsy)
{
  float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
  float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
  CVector2 sp0(spx,spy);
  CVector2 sp1(spx+dsx, spy+dsy);
  if( ibone_selected>=0 && ibone_selected<(int)aBone.size() ){
    CBone_RigMsh& bone = aBone[ibone_selected];
    bool is_updated = DragHandlerRot(bone.quat_joint,
                                     ielem_selected, sp0, sp1, bone.pos,
                                     mMV, mPj);
    if( is_updated ){
      this->UpdateBonePos();
    }
    /*
    if( ielem_selected>=0 && ielem_selected< 3){
      CVector3 v(0,0,0); v[ielem_selected] = 1;
      CVector3 axis = QuatVec(quat,v).Normalize();
      double ar = -DragCircle(sp0,sp1, CVector3(bone.pos), axis, mMV, mPj);
      double dq[4] = { cos(ar*0.5), v.x*sin(ar*0.5), v.y*sin(ar*0.5), v.z*sin(ar*0.5) };
      double qtmp[4]; QuatMult(qtmp, dq, bone.quat_joint);
      QuatCopy(bone.quat_joint,qtmp);
    }
     */
  }
}

void CRigMsh::Pick(double spx, double spy)
{
  float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
  float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
  CVector3 org = screenUnProjection(CVector3(spx,spy,0), mMV, mPj);
  CVector3 depth = screenDepthDirection(org,mMV,mPj);
  PickBone(ibone_selected, ielem_selected,
           aBone,
           org,depth,
           15,
           5);
//  std::cout << ibone_selected << " " << ielem_selected << std::endl;
}

void CMesh_RigMsh::SetSleketon(const std::vector<CBone_RigMsh>& aGBone){
  if( iskin_active < 0 || iskin_active >= (int)aSkin.size() ) return;
  CSkin_RigMsh& skin = aSkin[iskin_active];
  skin.SetSkeleton(aXYZ,aGBone,aXYZ_ini);
}

std::vector<std::string> CRigMsh::GetArrayTexPath() const {
  std::set<std::string> setTexPath;
  for(unsigned int imesh=0;imesh<aMesh.size();++imesh){
    for(unsigned int imaterial=0;imaterial<aMesh[imesh].aMaterial.size();++imaterial){
      const CMaterial_RigMsh& material = aMesh[imesh].aMaterial[imaterial];
      for(unsigned int itex=0;itex<material.aTexture_Diffuse.size();++itex){
        const CTextureInfo_RigMsh& mat_tex = material.aTexture_Diffuse[itex];
        std::string path = mat_tex.full_path;
        setTexPath.insert(path);
      }
    }
  }
  std::vector<std::string> aPath(setTexPath.begin(),setTexPath.end());
  return aPath;
}

void CRigMsh::UpdateBonePos()
{
  UpdateBoneRotTrans(aBone);
  for(unsigned int imesh=0;imesh<aMesh.size();++imesh){
    aMesh[imesh].SetSleketon(aBone);
  }
}

void CRigMsh::FixTexPath(const std::string& path_fbx)
{
  for(unsigned int imesh=0;imesh<aMesh.size();++imesh){
    CMesh_RigMsh& mesh = aMesh[imesh];
    for(unsigned int imat=0;imat<mesh.aMaterial.size();++imat){
      CMaterial_RigMsh& mat = mesh.aMaterial[imat];
      for(unsigned int itex=0;itex<mat.aTexture_Diffuse.size();++itex){
        std::string path_tex = mat.aTexture_Diffuse[itex].full_path;
//        std::cout << imesh << " " << imat << " " << path_tex << std::endl;
        if( isFileExists(path_tex) ){ continue; }
        std::string path_dir = getPathDir(path_fbx);
        if( path_tex.find('/') != std::string::npos ){
          int iloc = path_tex.find_last_of('/');
          std::string bn(path_tex.begin()+iloc+1,path_tex.end());
          mat.aTexture_Diffuse[itex].full_path = path_dir + "/" + bn;
        }
        else if( path_tex.find('\\') != std::string::npos ){
          int iloc = path_tex.find_last_of('\\');
          std::string bn(path_tex.begin()+iloc+1,path_tex.end());
          mat.aTexture_Diffuse[itex].full_path = path_dir + "/" + bn;
        }
        else{ // this is base name
          mat.aTexture_Diffuse[itex].full_path = path_dir + "/" + path_tex;
//          std::cout << mat.aTexture_Diffuse[itex].full_path << std::endl;
        }
        path_tex = mat.aTexture_Diffuse[itex].full_path;
        if( !isFileExists(path_tex) ){
          std::cout << "Error!-> broken tex info cannot fixed: " << path_tex << std::endl;
        }
      }
    }
  }
}

void CRigMsh::PrintInfo() const
{
  std::cout << "  number of mesh: " << aMesh.size() << std::endl;
  
  for(int imesh=0;imesh<(int)aMesh.size();++imesh){
    const CMesh_RigMsh& mesh = aMesh[imesh];
    std::cout << "  mesh: " << imesh << " " << mesh.aXYZ_ini.size()/3 << " " << mesh.aElemInd.size()-1 << " " << mesh.aElem.size() << std::endl;
    std::cout << "  layer info" << std::endl;
    for( int ilayer =0;ilayer<(int)mesh.aLayer.size();++ilayer){
      const CLayer_RigMsh& layer = mesh.aLayer[ilayer];
      std::cout << "    layer: " << ilayer << std::endl;
      std::cout << "      uv_setname: " << layer.uv_setname << std::endl;
      std::cout << "      mapping mode: " << layer.material_mapping_mode << std::endl;
      std::cout << "      mapping material face size: " << layer.material_mapping_face.size() << std::endl;
      std::cout << "      mapping material number: " << layer.material_mapping_mat2face.size() << std::endl;
    }
    
    std::cout << "  material info" << std::endl;
    for(unsigned int imaterial=0;imaterial<mesh.aMaterial.size();++imaterial){
      const CMaterial_RigMsh& material = mesh.aMaterial[imaterial];
      for(unsigned int itex=0;itex<material.aTexture_Diffuse.size();++itex){
        const CTextureInfo_RigMsh& mat_tex = material.aTexture_Diffuse[itex];
        std::cout << "    diffuse: " << imaterial << " " << itex << " " << mat_tex.uv_setname << " " << mat_tex.full_path << std::endl;
      }
    }
    std::cout << "  skin info" << std::endl;
    for(unsigned  int iskin = 0; iskin<mesh.aSkin.size();++iskin){
      const CSkin_RigMsh& skin = mesh.aSkin[iskin];
      for(unsigned int ibone=0;ibone<skin.aBone.size();++ibone){
        const CBoneWeight_RigMsh& bone = skin.aBone[ibone];
        std::cout << "    bone: " << ibone << " in skin: " << iskin << " has points: " << bone.aIpWeight.size() << " name: " << bone.name << std::endl;
      }
    }
  }
}



