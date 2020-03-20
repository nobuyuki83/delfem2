/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <set>
#include <map>
#include <cassert>
#include <sstream>
#include <fstream>
#include "delfem2/vec2.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/mat4.h"
#include "delfem2/quat.h"
// -----
#include "delfem2/v23m3q.h"
#include "delfem2/rig_v3q.h"

#ifndef M_PI 
#define M_PI 3.1415926535
#endif

namespace dfm2 = delfem2;

// ------------------------------------------------------------

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
    return aabb[0] <= aabb[1];
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

static void CalcInvMat(double* a, const int n, int& info )
{
  double tmp1;
  
  info = 0;
  int i,j,k;
  for(i=0;i<n;i++){
    if( fabs(a[i*n+i]) < 1.0e-30 ){
      info = 1;
      return;
    }
    if( a[i*n+i] < 0.0 ){
      info--;
    }
    tmp1 = 1.0 / a[i*n+i];
    a[i*n+i] = 1.0;
    for(k=0;k<n;k++){
      a[i*n+k] *= tmp1;
    }
    for(j=0;j<n;j++){
      if( j!=i ){
        tmp1 = a[j*n+i];
        a[j*n+i] = 0.0;
        for(k=0;k<n;k++){
          a[j*n+k] -= tmp1*a[i*n+k];
        }
      }
    }
  }
}

static std::string MyReplace
 (const std::string& str,
  const char cf,
  const char ct)
{
  const int n = str.size();
  //
  std::string ss(str);
  for(int i=0;i<n;++i){
    if( ss[i] != cf ){ continue; }
    ss[i] = ct;
  }
  return ss;
}

static std::vector<std::string> MySplit
 (const std::string& str,
  char delimiter)
{
  std::vector<std::string> aToken;
  aToken.clear();
  std::stringstream data(str);
  std::string line;
  while(std::getline(data,line,delimiter)){
    if( line.empty() ){ continue; }
    aToken.push_back(line);
  }
  return aToken;
}


// ------------------------------------------------------------

int dfm2::CRigBone::PickHandler
(const dfm2::CVec3d& org,
 const dfm2::CVec3d& dir,
 double rad_handlr,
 double tol) const
{
  return dfm2::PickHandlerRotation_Mat4(org,dir,
                                        affmat3Global, rad_handlr,
                                        tol);
}

void dfm2::CRigBone::SetRotationBryant
(double rx, double ry, double rz)
{
  dfm2::Quat_Bryant(quatRelativeRot, rx, ry, rz);
}

void dfm2::CRigBone::DeformSkin(double pos2[3],
                                const double pos0[3]){
  const double pos0a[4] = {pos0[0], pos0[1], pos0[2], 1.0};
  double pos1a[4]; MatVec4(pos1a,invBindMat,pos0a);
  double pos2a[4]; MatVec4(pos2a,affmat3Global,pos1a);
  pos2[0] = pos2a[0];
  pos2[1] = pos2a[1];
  pos2[2] = pos2a[2];
}

void dfm2::CRigBone::SetTranslation
(double tx, double ty, double tz)
{
  this->transRelative[0] = tx;
  this->transRelative[1] = ty;
  this->transRelative[2] = tz;
}

void dfm2::UpdateBoneRotTrans
(std::vector<dfm2::CRigBone>& aBone)
{
  for(std::size_t ibone=0;ibone<aBone.size();++ibone){
    const int ibone_p = aBone[ibone].ibone_parent;
    if( ibone_p < 0 || ibone_p >= (int)aBone.size() ){ // root bone
      dfm2::Mat4_ScaleRotTrans(aBone[ibone].affmat3Global,
                               aBone[ibone].scale,
                               aBone[ibone].quatRelativeRot,
                               aBone[ibone].transRelative);
      continue;
    }
    double M01[16];
    dfm2::Mat4_ScaleRotTrans(M01,
                             aBone[ibone].scale,
                             aBone[ibone].quatRelativeRot,
                             aBone[ibone].transRelative);
    dfm2::MatMat4(aBone[ibone].affmat3Global,
                  aBone[ibone_p].affmat3Global,M01);
  }
}



void dfm2::UpdateRigSkin
(double* aXYZ,
 const double* aXYZ0,
 unsigned int nXYZ,
 const unsigned int* aTri,
 unsigned int nTri,
 const std::vector<dfm2::CRigBone>& aBone,
 const double* aRigWeight,
 const unsigned int* aRigJoint)
{
  for(unsigned int ip=0;ip<nXYZ;++ip){
    double pos0[4] = {aXYZ0[ip*3+0],aXYZ0[ip*3+1],aXYZ0[ip*3+2],1.0};
    double pos1[3] = {0,0,0};
    double sum_w = 0.0;
    for(int iij=0;iij<4;++iij){
      double w = aRigWeight[ip*4+iij];
      if( w < 1.0e-30 ){ continue; }
      unsigned int ij = aRigJoint[ip*4+iij];
      sum_w += w;
      assert (ij<aBone.size());
      double pos0a[4]; dfm2::MatVec4(pos0a,aBone[ij].invBindMat,pos0);
      double pos0b[4]; dfm2::MatVec4(pos0b,aBone[ij].affmat3Global,pos0a);
      pos1[0] += w*pos0b[0];
      pos1[1] += w*pos0b[1];
      pos1[2] += w*pos0b[2];
    }
    assert( fabs(sum_w)>1.0e-10 );
    pos1[0] /= sum_w;
    pos1[1] /= sum_w;
    pos1[2] /= sum_w;
    aXYZ[ip*3+0] = pos1[0];
    aXYZ[ip*3+1] = pos1[1];
    aXYZ[ip*3+2] = pos1[2];
  }
}


// ------------------------------------
// from here BioVisionHierarchy

void dfm2::Read_BioVisionHierarchy
(std::vector<dfm2::CRigBone>& aBone,
 std::vector<dfm2::CChannel_BioVisionHierarchy>& aChannelRotTransBone,
 int& nframe,
 std::vector<double>& aValueRotTransBone,
 const std::string& path_bvh)
{
  std::ifstream fin;
  fin.open(path_bvh.c_str());
  if( !fin.is_open() ){
    std::cout << "cannot open file" << std::endl;
    return;
  }
  aBone.clear();
  aChannelRotTransBone.clear();
  //
  std::string line;
  std::vector<int> stackIndBone;
  while(std::getline(fin,line)){
    if (line[line.size()-1] == '\n') line.erase(line.size()-1); // remove the newline code
    if (line[line.size()-1] == '\r') line.erase(line.size()-1); // remove the newline code
    line = MyReplace(line, '\t', ' ');
    std::vector<std::string> aToken = MySplit(line,' ');
//    std::cout << aToken[0] << std::endl;
    if( aToken[0] == "HIERARCHY" ){
      assert(aBone.empty());
    }
    else if( aToken[0] == "ROOT" ){
      assert(aBone.size()==0);
      dfm2::CRigBone br;
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
      double org_x = myStod(aToken[1]);
      double org_y = myStod(aToken[2]);
      double org_z = myStod(aToken[3]);
      aBone[ib].invBindMat[ 3] = -org_x;
      aBone[ib].invBindMat[ 7] = -org_y;
      aBone[ib].invBindMat[11] = -org_z;
      if( stackIndBone.size() > 1 ){
        const int ibp = stackIndBone[stackIndBone.size()-2];
        assert(ibp<(int)aBone.size());
        aBone[ib].invBindMat[ 3] += aBone[ibp].invBindMat[ 3];
        aBone[ib].invBindMat[ 7] += aBone[ibp].invBindMat[ 7];
        aBone[ib].invBindMat[11] += aBone[ibp].invBindMat[11];
      }
    }
    else if( aToken[0] == "CHANNELS" ){
      assert(aToken.size()>=2);
      int nch = myStoi(aToken[1]);
      assert((int)aToken.size()==nch+2);
      assert( !aBone.empty() );
      const std::size_t ib = aBone.size()-1;
      for(int ich=0;ich<nch;++ich){
        const std::string& type_ch = aToken[ich+2];
        if(      type_ch == "Xposition" ){ aChannelRotTransBone.emplace_back(ib,0,false ); }
        else if( type_ch == "Yposition" ){ aChannelRotTransBone.emplace_back(ib,1,false ); }
        else if( type_ch == "Zposition" ){ aChannelRotTransBone.emplace_back(ib,2,false ); }
        else if( type_ch == "Xrotation" ){ aChannelRotTransBone.emplace_back(ib,0,true ); }
        else if( type_ch == "Yrotation" ){ aChannelRotTransBone.emplace_back(ib,1,true ); }
        else if( type_ch == "Zrotation" ){ aChannelRotTransBone.emplace_back(ib,2,true ); }
        else{
          std::cout << "ERROR-->undefiend type" << std::endl;
        }
      }
    }
    else if( aToken[0] == "JOINT" ){
      dfm2::CRigBone br;
      assert( aToken.size() == 2 );
      br.name = aToken[1];
      aBone.push_back(br);
    }
    else if( aToken[0] == "End" ){
      assert(aToken[1] == "Site");
      dfm2::CRigBone br;
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
//      std::cout << "frame: " << nframe << std::endl;
    }
    std::getline(fin,line);
//    std::cout << "frametime: " << line << std::endl;
  }
  const int nchannel = aChannelRotTransBone.size();
  aValueRotTransBone.resize(nframe*nchannel);
  for(int iframe=0;iframe<nframe;++iframe){
    std::getline(fin,line);
    line = MyReplace(line, '\t', ' ');
    if (line[line.size()-1] == '\n') line.erase(line.size()-1); // remove the newline code
    if (line[line.size()-1] == '\r') line.erase(line.size()-1); // remove the newline code
    std::vector<std::string> aToken = MySplit(line,' ');
//    std::cout << aToken.size() << " " << aChannelRotTransBone.size() << std::endl;
    assert(aToken.size()==aChannelRotTransBone.size());
    for(int ich=0;ich<nchannel;++ich){
      aValueRotTransBone[iframe*nchannel+ich] = myStod(aToken[ich]);
    }
  }
  // ---------------
  for(std::size_t ibone=0;ibone<aBone.size();++ibone){
    dfm2::CRigBone& bone = aBone[ibone];
    bone.scale = 1.0;
    bone.quatRelativeRot[0] = 1.0;
    bone.quatRelativeRot[1] = 0.0;
    bone.quatRelativeRot[2] = 0.0;
    bone.quatRelativeRot[3] = 0.0;
    bone.transRelative[0] = 0.0;
    bone.transRelative[1] = 0.0;
    bone.transRelative[2] = 0.0;
    if( bone.ibone_parent != -1 ){
      const dfm2::CRigBone& bone_p = aBone[bone.ibone_parent];
      bone.transRelative[0] = (-bone.invBindMat[ 3])-(-bone_p.invBindMat[ 3]);
      bone.transRelative[1] = (-bone.invBindMat[ 7])-(-bone_p.invBindMat[ 7]);
      bone.transRelative[2] = (-bone.invBindMat[11])-(-bone_p.invBindMat[11]);
    }
  }
  for(auto & bone : aBone){
    for(int i=0;i<16;++i){ bone.affmat3Global[i] = bone.invBindMat[i]; }
    int info; CalcInvMat(bone.affmat3Global, 4, info);
  }
}


void dfm2::SetPose_BioVisionHierarchy
(std::vector<dfm2::CRigBone>& aBone,
 const std::vector<dfm2::CChannel_BioVisionHierarchy>& aChannelRotTransBone,
 const double *aVal)
{
  for(auto & bone : aBone){
    bone.quatRelativeRot[0] = 1.0;
    bone.quatRelativeRot[1] = 0.0;
    bone.quatRelativeRot[2] = 0.0;
    bone.quatRelativeRot[3] = 0.0;
  }
  const int nch = aChannelRotTransBone.size();
  for(int ich=0;ich<nch;++ich){
    const int ibone = aChannelRotTransBone[ich].ibone;
    const int iaxis = aChannelRotTransBone[ich].iaxis;
    const bool isrot = aChannelRotTransBone[ich].isrot;
    const double val = aVal[ich];
    assert(ibone<(int)aBone.size());
    assert(iaxis>=0&&iaxis<3);
    if( !isrot ){
      aBone[ibone].transRelative[iaxis] = val;
    }
    else{
      const double ar = val*M_PI/180.0;
      double v0[3] = {0,0,0};
      v0[iaxis] = 1.0;
      double dq[4] = { cos(ar*0.5), v0[0]*sin(ar*0.5), v0[1]*sin(ar*0.5), v0[2]*sin(ar*0.5) };
      double qtmp[4]; dfm2::QuatQuat(qtmp,
                                     aBone[ibone].quatRelativeRot, dq);
      dfm2::Copy_Quat(aBone[ibone].quatRelativeRot,qtmp);
    }
  }
  dfm2::UpdateBoneRotTrans(aBone);
}

// ----------------------------------

void dfm2::Smpl2Rig(
              std::vector<dfm2::CRigBone>& aBone,
              const std::vector<int>& aIndBoneParent,
              const std::vector<double>& aXYZ0,
              const std::vector<double>& aJntRgrs)
{
  const unsigned int nbone = aIndBoneParent.size();
  std::vector<double> aJntPos0;
  {
    const unsigned int nP = aXYZ0.size()/3;
    const unsigned int nBone = aIndBoneParent.size();
    aJntPos0.assign(nBone*3, 0.0);
    for(unsigned int ib=0;ib<nBone;++ib){
      aJntPos0[ib*3+0] = 0;
      aJntPos0[ib*3+1] = 0;
      aJntPos0[ib*3+2] = 0;
      for(unsigned int ip=0;ip<nP;++ip){
        aJntPos0[ib*3+0] += aJntRgrs[ip*nBone+ib]*aXYZ0[ip*3+0];
        aJntPos0[ib*3+1] += aJntRgrs[ip*nBone+ib]*aXYZ0[ip*3+1];
        aJntPos0[ib*3+2] += aJntRgrs[ip*nBone+ib]*aXYZ0[ip*3+2];
      }
    }
  }
  aBone.resize(nbone);
  for(unsigned int ib=0;ib<nbone;++ib){
    int ibp = aIndBoneParent[ib];
    aBone[ib].ibone_parent = ibp;
    aBone[ib].invBindMat[ 3] = -aJntPos0[ib*3+0];
    aBone[ib].invBindMat[ 7] = -aJntPos0[ib*3+1];
    aBone[ib].invBindMat[11] = -aJntPos0[ib*3+2];
    if( ibp != -1 ){
      aBone[ib].transRelative[0] = +aJntPos0[ib*3+0] - aJntPos0[ibp*3+0];
      aBone[ib].transRelative[1] = +aJntPos0[ib*3+1] - aJntPos0[ibp*3+1];
      aBone[ib].transRelative[2] = +aJntPos0[ib*3+2] - aJntPos0[ibp*3+2];
    }
    else{
      aBone[ib].transRelative[0] = +aJntPos0[ib*3+0];
      aBone[ib].transRelative[1] = +aJntPos0[ib*3+1];
      aBone[ib].transRelative[2] = +aJntPos0[ib*3+2];
    }
  }
  dfm2::UpdateBoneRotTrans(aBone);
}

// -------------------------------------------


void PickBone
(int& ibone_selected,
 int& ielem_selected,
 const std::vector<dfm2::CRigBone>& aBone,
 const dfm2::CVec3d& src,
 const dfm2::CVec3d& dir,
 double rad_hndlr,
 double tol)
{
  if( ibone_selected>=0 && ibone_selected<(int)aBone.size() ){
    const dfm2::CRigBone& bone = aBone[ibone_selected];
    ielem_selected = bone.PickHandler(src,dir,rad_hndlr,tol);
  }
  else{
    ielem_selected = -1;
  }
  if( ielem_selected == -1 ){
    ibone_selected = -1;
    for(int ibone=0;ibone<(int)aBone.size();++ibone){
      delfem2::CVec3d pos(aBone[ibone].Pos());
      double distance = Distance(nearest_Line_Point(pos, src, dir),pos);
      if( distance < tol ){
        ibone_selected = ibone;
        break;
      }
    }
  }
}



void dfm2::SetMat4AffineBone_FromJointRelativeRotation
 (std::vector<double>& aMat4AffineBone,
  const double trans_root[3],
  const std::vector<double>& aQuatRelativeRot,
  const std::vector<int>& aIndBoneParent,
  const std::vector<double>& aJntPos0)
{
  const unsigned int nBone = aIndBoneParent.size();
  assert( nBone >= 1 );
  assert( aMat4AffineBone.size() == nBone*16 );
  dfm2::Mat4_ScaleRotTrans(aMat4AffineBone.data(),
                           1.0, aQuatRelativeRot.data(), trans_root);
  for(unsigned int ibone=1;ibone<nBone;++ibone){
    int ibp = aIndBoneParent[ibone];
    assert( ibp >= 0 && ibp < (int)nBone );
    // inv binding mat
    double p1[3] = {aJntPos0[ibone*3+0], aJntPos0[ibone*3+1], aJntPos0[ibone*3+2]};
    dfm2::CMat4<double> M0, M1, M2;
    M0.Set_AffineTranslate(-p1[0], -p1[1], -p1[2]);
    dfm2::Mat4_Quat(M1.mat,
                    aQuatRelativeRot.data()+ibone*4);
    M2.Set_AffineTranslate(+p1[0], +p1[1], +p1[2]);
    dfm2::CMat4<double> M3 = M1.MatMat(M0);
    dfm2::CMat4<double> M4 = M2.MatMat(M3);
    dfm2::MatMat4(aMat4AffineBone.data()+ibone*16,
                  aMat4AffineBone.data()+ibp*16,
                  M4.mat);
  }
}



/*
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
*/
/*
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
        A += Mat3_OuterProduct(a,b);
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
 */
/*
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
 */

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

