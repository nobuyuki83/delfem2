/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "gtest/gtest.h"

#include "delfem2/quat.h"
#include "delfem2/geo3_v23m34q.h"

// Define these only in *one* .cc file.
//#define TINYGLTF_IMPLEMENTATION
//#define STB_IMAGE_IMPLEMENTATION
//#define STB_IMAGE_WRITE_IMPLEMENTATION
// #define TINYGLTF_NOEXCEPTION // optional. disable exception handling.
#include "tinygltf/tiny_gltf.h"

#include "delfem2/rig_geo3.h"
#include "delfem2/tinygltf/io_gltf.h"

namespace dfm2 = delfem2;

// ----------------------------

TEST(gltf,formatcheck)
{
  tinygltf::TinyGLTF loader;
  for(int ifile=0;ifile<4;++ifile){
    std::string path_glb = std::string(PATH_INPUT_DIR);
    if(      ifile == 0 ){ path_glb += "/Duck.glb"; }
    else if( ifile == 1 ){ path_glb += "/RiggedSimple.glb"; }
    else if( ifile == 2 ){ path_glb += "/RiggedFigure.glb"; }
    else if( ifile == 3 ){ path_glb += "/CesiumMan.glb"; }
    tinygltf::Model model;
    std::string err;
    std::string warn;
    bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, path_glb); // for binary glTF(.glb)
    EXPECT_TRUE(warn.empty());
    EXPECT_TRUE(err.empty());
    EXPECT_TRUE(ret);
    for(int in=0;in<model.nodes.size();++in){
      const tinygltf::Node& node = model.nodes[in];
      EXPECT_TRUE( node.rotation.empty() || node.rotation.size() == 4 );
      EXPECT_TRUE( node.translation.empty() || node.translation.size() == 3 );
      EXPECT_TRUE( node.scale.empty() || node.scale.size() == 3 );
      EXPECT_TRUE( node.matrix.empty() || node.matrix.size() == 16 );
      EXPECT_FALSE( !node.matrix.empty() && !node.scale.empty() );      // if there is matrix, no scale
      EXPECT_FALSE( !node.matrix.empty() && !node.translation.empty() );     // if there is matrix, no translation
      EXPECT_FALSE( !node.matrix.empty() && !node.rotation.empty() );     // if there is matrix, no rotation
      if( node.skin != -1 ){
        EXPECT_TRUE( node.skin>=0 && node.skin<model.skins.size() );
      }
      if( node.mesh != -1 ){
        EXPECT_TRUE( node.mesh>=0 && node.mesh<model.meshes.size() );
      }
      for(int ic : node.children){
        EXPECT_TRUE( ic>=0 && ic<model.nodes.size() );
      }
    }
    for(int is=0;is<model.skins.size();++is){
      const tinygltf::Skin& skin = model.skins[is];
      for(int inode : skin.joints){
        EXPECT_TRUE( inode>=0 && inode < model.nodes.size() );
      }
      EXPECT_TRUE( skin.inverseBindMatrices>=0 && skin.inverseBindMatrices<model.accessors.size() );
      const tinygltf::Accessor& acc = model.accessors[skin.inverseBindMatrices];
      EXPECT_EQ( acc.type, TINYGLTF_TYPE_MAT4 );
      EXPECT_EQ( acc.componentType, TINYGLTF_COMPONENT_TYPE_FLOAT );
      EXPECT_EQ( acc.count, skin.joints.size() );
    }
  }
}

TEST(gltf,io_gltf_skin_sensitivity)
{
  std::random_device rd;
  std::mt19937 rndeng(rd());
  std::uniform_real_distribution<double> dist_01(+0,+1);
  //
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aTri;
  std::vector<dfm2::CRigBone> aBone;
  std::vector<double> aW;
  {
    std::vector<double> aRigWeight; // [np, 4]
    std::vector<unsigned int> aRigJoint; // [np, 4]
    {
      //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Duck.glb";
      //    std::string path_glb = std::string(PATH_INPUT_DIR)+"/Monster.glb";
      
      //      std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedSimple.glb";
      //    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedFigure.glb";
      std::string path_glb = std::string(PATH_INPUT_DIR)+"/CesiumMan.glb";
      dfm2::CGLTF gltf;
      gltf.Read(path_glb);
  //    gltf.Print();
      gltf.GetMeshInfo(aXYZ0, aTri, aRigWeight, aRigJoint, 0,0);
      gltf.GetBone(aBone, 0);
    }
    {
      for(auto & bone : aBone){
        dfm2::Quat_Identity(bone.quatRelativeRot);
      }
      UpdateBoneRotTrans(aBone);
      std::vector<double> aXYZ = aXYZ0;
      Skinning_LBS_LocalWeight(aXYZ.data(),
                               aXYZ0.data(), aXYZ0.size()/3,
                               aTri.data(), aTri.size()/3,
                               aBone, aRigWeight.data(), aRigJoint.data());
    }
    const unsigned int np = aXYZ0.size()/3;
    const unsigned int nb = aBone.size();
    aW.assign(np*nb,0.0);
    for(int ip=0;ip<np;++ip){
      for(int iib=0;iib<4;++iib){
        unsigned int ib = aRigJoint[ip*4+iib];
        aW[ip*nb+ib] = aRigWeight[ip*4+iib];
      }
    }
  }
  
  const unsigned int nb = aBone.size();
  assert( aW.size() == aXYZ0.size()/3*nb );
  
  // ------------
  std::vector<double> Lx, Ly, Lz;  // [ nsns, nb*4 ]
  for(int ibs=0;ibs<aBone.size();++ibs){
    for(int idims=0;idims<3;++idims){
      dfm2::Rig_SensitivityBoneTransform_Eigen(Lx,Ly,Lz,
                                               ibs,idims,true,
                                               aBone);
    }
  }
  for(int idims=0;idims<3;++idims){
    dfm2::Rig_SensitivityBoneTransform_Eigen(Lx,Ly,Lz,
                                             0,idims,false,
                                             aBone);
  }
  // ---------------
  
  const double eps = 1.0e-4;
  
  { // Check Sensitivity Skin
    for(auto & bone : aBone){
      dfm2::Quat_Identity(bone.quatRelativeRot);
    }
    UpdateBoneRotTrans(aBone);
    std::vector<double> aRefPos; // [ np, nBone*4 ]
    Rig_SkinReferncePositionsBoneWeighted(aRefPos,
                                          aBone,aXYZ0,aW);
    const unsigned int nsns = Lx.size()/(nb*4);
    assert( nsns==(nb+1)*3 );
    for(int isns=0;isns<nsns;++isns){
      unsigned int ib_s = isns/3;
      bool is_rot = true;
      unsigned int idim_s = isns - ib_s*3;
      if( ib_s == nb ){ ib_s = 0; is_rot = false; }
      std::vector<dfm2::CRigBone> aBone2 = aBone;
      if( is_rot ){
        dfm2::CQuatd dq = dfm2::Quat_CartesianAngle(eps*dfm2::CVec3d::Axis(idim_s));
        dfm2::CQuatd q0 = dq*dfm2::CQuatd(aBone2[ib_s].quatRelativeRot);
        q0.CopyTo(aBone2[ib_s].quatRelativeRot);
      }
      else{
        aBone2[ib_s].transRelative[idim_s] += eps;
      }
      std::vector<double> aXYZ1;
      dfm2::UpdateBoneRotTrans(aBone);
      dfm2::Skinning_LBS(aXYZ1,
                         aXYZ0, aBone, aW);
      // ----------------
      std::vector<double> aXYZ2;
      dfm2::UpdateBoneRotTrans(aBone2);
      dfm2::Skinning_LBS(aXYZ2,
                         aXYZ0, aBone2, aW);
      const unsigned int np = aXYZ0.size()/3;
      for(unsigned int ip=0;ip<np;++ip){
        const double val0[3] = {
          (aXYZ2[ip*3+0] - aXYZ1[ip*3+0])/eps,
          (aXYZ2[ip*3+1] - aXYZ1[ip*3+1])/eps,
          (aXYZ2[ip*3+2] - aXYZ1[ip*3+2])/eps };
        double val1[3] =  { 0, 0, 0 };
        for(unsigned int j=0;j<nb*4;++j){
          val1[0] += aRefPos[ip*(nb*4)+j]*Lx[isns*(nb*4)+j];
          val1[1] += aRefPos[ip*(nb*4)+j]*Ly[isns*(nb*4)+j];
          val1[2] += aRefPos[ip*(nb*4)+j]*Lz[isns*(nb*4)+j];
        }
        for(int i=0;i<3;++i){
          EXPECT_NEAR(val0[i],val1[i], 1.0e-3*(fabs(val1[i])+1.0) );
        }
      }
    }
  }
  
  { // check sensitivity target
    for(auto & bone : aBone){
      dfm2::Quat_Identity(bone.quatRelativeRot);
    }
    dfm2::UpdateBoneRotTrans(aBone);
    // ----------------------
    std::vector<dfm2::CTarget> aTarget;
    srand(0);
    for(int itr=0;itr<5;++itr){
      dfm2::CTarget t;
      t.ib = (unsigned int)(aBone.size()*(rand()/(RAND_MAX+1.0)));
      t.pos = aBone[t.ib].Pos() + dfm2::CVec3d::Random(dist_01,rndeng);
      aTarget.push_back(t);
    }
    
    std::vector<double> aO0; // [nC]
    std::vector<double> adO0; // [nC, nb*3 ]
    for(auto & it : aTarget){
      dfm2::Rig_WdW_Target_Eigen(aO0,adO0,
                                 aBone,it,Lx,Ly,Lz);
    }
    
    const unsigned int nsns = Lx.size()/(nb*4);
    assert( nsns==(nb+1)*3 );
    for(unsigned int isns=0;isns<nsns;++isns){
      unsigned int ib_s = isns/3;
      bool is_rot = true;
      unsigned int idim_s = isns - ib_s*3;
      if( ib_s == nb ){ ib_s = 0; is_rot = false; }
      std::vector<dfm2::CRigBone> aBone1 = aBone;
      if( is_rot ){
        dfm2::CQuatd dq = dfm2::Quat_CartesianAngle(eps*dfm2::CVec3d::Axis(idim_s));
        dfm2::CQuatd q0 = dq*dfm2::CQuatd(aBone1[ib_s].quatRelativeRot);
        q0.CopyTo(aBone1[ib_s].quatRelativeRot);
      }
      else{
        aBone1[ib_s].transRelative[idim_s] += eps;
      }
      dfm2::UpdateBoneRotTrans(aBone1);
      // -------------
      std::vector<double> aO1; // [nC]
      std::vector<double> adO1; // [nC, nb*3 ]
      for(auto & it : aTarget){
        dfm2::Rig_WdW_Target_Eigen(aO1,adO1,
            aBone1,it,Lx,Ly,Lz);
      }
      // -------------
      for(int io=0;io<aO0.size();++io){
        EXPECT_NEAR((aO1[io]-aO0[io])/eps, adO0[io*nsns+isns], 0.5*eps*(fabs(adO0[io*nsns+isns])+1.0));
      }
    }
  }
}
