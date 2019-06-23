/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef RIGMESH_H
#define RIGMESH_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>

class CVector3;

class CBone_RigMsh
{
public:
    CBone_RigMsh(){
      is_active = true;
      quat_joint[0] = 1;
      quat_joint[1] = 0;
      quat_joint[2] = 0;
      quat_joint[3] = 0;
      quat[0] = 1;
      quat[1] = 0;
      quat[2] = 0;
      quat[3] = 0;
      ibone_parent = -1;
    }
  void Draw(bool is_selected,
            int ielem_select,
            const std::vector<CBone_RigMsh>& aBone,
            double len) const;
  int PickHandler(const CVector3& org, const CVector3& dir,
                  double rad_handlr,
                  double tol) const;
  void Affine(const double a[16]);
public:
  std::string name;
  int ibone_parent;
  bool is_active;
  /////
  double pos_ini[3];
  double quat_joint[4];
  /////
  double pos[3];
  double quat[4]; // joint rotation at this bone_point. Rotation propageate to the child skeleton
};

class CChannel_RotTransBone_BVH
{
public:
  CChannel_RotTransBone_BVH(int ib, int ia, bool br){
    this->ibone = ib;
    this->iaxis = ia;
    this->isrot = br;
  }
public:
  int ibone;
  int iaxis;
  bool isrot;
};

void DrawBone(const std::vector<CBone_RigMsh>& aBone,
              int ibone_selected,
              int ielem_selected,
              double bone_rad);

void ReadBVH(std::vector<CBone_RigMsh>& aBone,
             std::vector<CChannel_RotTransBone_BVH>& aChannelRotTransBone,
             int& nframe,
             std::vector<double>& aRotTransBone,
             const std::string& path_bvh);

void InitializeBone(std::vector<CBone_RigMsh>& aBone);

void SetRotTransBVH(std::vector<CBone_RigMsh>& aBone,
                    const std::vector<CChannel_RotTransBone_BVH>& aChannelRotTransBone,
                    const double *aVal);

void UpdateBoneRotTrans(std::vector<CBone_RigMsh>& aBone);

void PickBone(int& ibone_selected,
              int& ielem_selected,
              const std::vector<CBone_RigMsh>& aBone,
              const CVector3& src,
              const CVector3& dir,
              double rad_hndlr,
              double tol);

class CBoneGoal
{
public:
  CBoneGoal(){
    itype = 0; // 0:position, 1:on line
  }
  void GetGoalPos(double* pos_trg,
                  const double* org_rot,
                  const double* pos_cur) const;
public:
  int ibone;
  ////
  int itype;
  double pos[3];
  double dir[3];
};

void BoneOptimization(std::vector<CBone_RigMsh>& aBone,
                      const std::vector<CBoneGoal>& aBoneGoal);

void DrawBoneTarget(const std::vector<CBoneGoal>& aBoneGoal,
                    const std::vector<CBone_RigMsh>& aBone);

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

class CTexManager; // implemented in the funcs_gl.h

class CTextureInfo_RigMsh
{
public:
  std::string uv_setname;
  std::string name;
  std::string full_path;
};

class CMaterial_RigMsh
{
public:
  double RGB_Diffuse[3];
  std::vector<CTextureInfo_RigMsh> aTexture_Diffuse;
};

class CBoneWeight_RigMsh
{
public:
  std::string name;
  std::vector< std::pair<int,double> > aIpWeight;
};

class CSkin_RigMsh
{
public:
  int findBone(const std::string& name_bone) const {
    int ibone = 0;
    for(;ibone<(int)aBone.size();++ibone){
      if( aBone[ibone].name == name_bone ){ return ibone; }
    }
    return -1;
  }
  void computeWeight(std::vector<double>& aWeight,
                     const std::string& name_bone_active, int npoint) const
  {
    aWeight.assign(npoint,0);
    int ibone_active = findBone(name_bone_active);
    if( ibone_active < 0 || ibone_active > (int)aBone.size() ){ return; }
    for(int ipoint=0;ipoint<npoint;++ipoint){
      for(int ipobo=aPoint2BoneInd[ipoint];ipobo<aPoint2BoneInd[ipoint+1];++ipobo){
        int ibone0 = aPoint2Bone[ipobo];
        if( ibone0 != ibone_active ) continue;
        double weight = aPoint2BoneWeight[ipobo];
        aWeight[ipoint] = weight;
      }
    }
  }
  void Finalize(const int npoint);
  void SetSkeleton(std::vector<double>& aXYZ,
                   const std::vector<CBone_RigMsh>& aSkeleton,
                   const std::vector<double>& aXYZ_ini);
public:
  std::vector<CBoneWeight_RigMsh> aBone;
  std::vector<int> aPoint2BoneInd;
  std::vector<int> aPoint2Bone;
  std::vector<double> aPoint2BoneWeight;
};


class CLayer_RigMsh
{
public:
    CLayer_RigMsh(){
      material_mapping_mode = "";
    }
  void Finalize(){
    int nmat = 0;
    for(int iface=0;iface<(int)material_mapping_face.size();++iface){
      int imat = material_mapping_face[iface];
      if( imat >= nmat ){
        nmat = imat+1;
      }
    }
    material_mapping_mat2face.resize(nmat);
    for(int iface=0;iface<(int)material_mapping_face.size();++iface){
      int imat = material_mapping_face[iface];
      material_mapping_mat2face[imat].push_back(iface);
    }
  }
public:
  std::string uv_setname;
  std::string uv_bindmode;
  std::vector<double> aUV;
  /////
  std::string material_mapping_mode;
  std::vector<int> material_mapping_face;
  std::vector< std::vector<int> > material_mapping_mat2face;
};

class CMesh_RigMsh
{
public:
    CMesh_RigMsh(){
      ilayer_active = 0;
      iskin_active = 0;
    }
  bool isTextureWithUVSetName(const std::string& str) const{
    for(int imat=0;imat<(int)aMaterial.size();++imat){
      const CMaterial_RigMsh& mat = aMaterial[imat];
      for(int itex=0;itex<(int)mat.aTexture_Diffuse.size();++itex){
        const CTextureInfo_RigMsh& tex = mat.aTexture_Diffuse[itex];
        if( tex.uv_setname == "default" ){ return true; }
        if( tex.uv_setname == str ){ return true; }
      }
    }
    return false;
  }
  const CTextureInfo_RigMsh& getTextureWithUVSetName(const std::string& str) const{
    for(int imat=0;imat<(int)aMaterial.size();++imat){
      const CMaterial_RigMsh& mat = aMaterial[imat];
      for(int itex=0;itex<(int)mat.aTexture_Diffuse.size();++itex){
        const CTextureInfo_RigMsh& tex = mat.aTexture_Diffuse[itex];
        if( tex.uv_setname == "default" ){ return tex; }
        if( tex.uv_setname == str ){ return tex; }
      }
    }
    assert(0);
    abort();
  }
  void DrawLayerWithTex(int ilayer,const CTexManager& tex_manager, bool is_ini) const;
  void SetActiveBone(const std::string& name_bone_active){
    if( iskin_active < 0 || iskin_active >= (int)aSkin.size() ){ return; }
    const CSkin_RigMsh& skin = aSkin[iskin_active];
    skin.computeWeight(aWeight,name_bone_active,aXYZ_ini.size()/3);
  }
  void SetSleketon(const std::vector<CBone_RigMsh>& aSkeleton); // initialize aWeight and aXYZ
  void Affine(const double a[16]);
public:
  std::vector<double> aXYZ_ini;
  std::vector<unsigned int> aElemInd,aElem;
  std::vector<CLayer_RigMsh> aLayer;
  std::vector<CMaterial_RigMsh> aMaterial;
  std::vector<CSkin_RigMsh> aSkin;
  int ilayer_active;
  int iskin_active;
  /////
  std::vector<double> aXYZ;
  std::vector<double> aWeight; // weigh of a bone for ilayer_active, iskin_active, name_bone_active
};


class CRigMsh
{
public:
  CRigMsh(){
    ibone_selected = 1;
    ielem_selected = 0;
    draw_rep_length = -1;
    is_draw_bone = true;
    is_draw_weight = false;
    ////
    color_bone_weight.assign(4, 0.0);
    color_bone_weight[0] = 1.0;
    color_bone_weight[3] = 1.0;
    color_bone_weight_back.assign(4, 0.8);
    color_bone_weight_back[3] = 1.0;
  }
  void Draw(const CTexManager& tex_manager) const;
  std::vector<double> MinMaxXYZ() const;
  ///////////
  void Affine(const double A[16]);
  /*
    for(unsigned int im=0;im<aMesh.size();++im){ aMesh[im].Scale(s); }
    for(unsigned int ib=0;ib<aBone.size();++ib){ aBone[ib].Scale(s); }
    draw_rep_length *= s;
  }
  void Translate(double dx, double dy, double dz){
    for(unsigned int im=0;im<aMesh.size();++im){ aMesh[im].Translate(dx,dy,dz); }
    for(unsigned int ib=0;ib<aBone.size();++ib){ aBone[ib].Translate(dx,dy,dz); }
  }
   */
  void Rotate(double dx, double dy, double dz);
  ///////////
  void Clear(){
    ibone_selected = -1;
    ielem_selected = -1;
    is_draw_bone = true;
    aMesh.clear();
    aBone.clear();
  }
  void Initialize(){
    const std::vector<double>& mm = this->MinMaxXYZ();
    double lx = mm[1]-mm[0];
    double ly = mm[3]-mm[2];
    double lz = mm[5]-mm[4];
    draw_rep_length = sqrt(lx*lx+ly*ly+lz*lz);
  }
  void UpdateBonePos();
  std::vector<std::string> GetArrayTexPath() const;
  void Pick(double spx, double spy);
  void Drag(double spx, double spy, double dsx, double dsy);
  void FixTexPath(const std::string& path_fbx);
  void PrintInfo() const;
  void DisplayBoneWeightOnMesh(int ibone){
    if( ibone < 0 || ibone >= (int)aBone.size() ){
      is_draw_weight = false;
      return;
    }
    ///////
    is_draw_weight = true;
    for(int imesh=0;imesh<(int)this->aMesh.size();++imesh){
      CMesh_RigMsh& mesh = aMesh[imesh];
      mesh.SetActiveBone(aBone[ibone].name);
    }
  }
private:
  bool is_draw_weight;
public:
  int ielem_selected;
  int ibone_selected;
  ////
  double draw_rep_length; // diagonal length of the bounding box，set in the Initialize()
  bool is_draw_bone;
  std::vector<CMesh_RigMsh> aMesh;
  std::vector<CBone_RigMsh> aBone;
  ////
  std::vector<double> color_bone_weight;
  std::vector<double> color_bone_weight_back;
};

#endif /* rigmesh_hpp */
