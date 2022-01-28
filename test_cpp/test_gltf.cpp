/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "gtest/gtest.h"

#include "delfem2/tinygltf/io_gltf.h"
#include "delfem2/rigopt.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/quat.h"
#include "delfem2/sampling.h"

// Define these only in *one* .cc file.
//#define TINYGLTF_IMPLEMENTATION
//#define STB_IMAGE_IMPLEMENTATION
//#define STB_IMAGE_WRITE_IMPLEMENTATION
// #define TINYGLTF_NOEXCEPTION // optional. disable exception handling.
#include "tinygltf/tiny_gltf.h"


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

TEST(gltf, skin_sensitivity) {
  std::random_device rd;
  std::mt19937 rndeng(rd());
  std::uniform_real_distribution<double> dist_01(+0, +1);
  //
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aTri;
  std::vector<dfm2::CRigBone> aBone0;
  std::vector<double> aSkinSparseW; // [np, 4]
  std::vector<unsigned int> aSkinSparseI; // [np, 4]
  {
    std::string path_glb = std::string(PATH_INPUT_DIR) + "/CesiumMan.glb";
    dfm2::CGLTF gltf;
    gltf.Read(path_glb);
    //    gltf.Print();
    gltf.GetMeshInfo(aXYZ0, aTri, aSkinSparseW, aSkinSparseI, 0, 0);
    gltf.GetBone(aBone0, 0);
    dfm2::SetCurrentBoneRotationAsDefault(aBone0);
  }
  const size_t nb = aBone0.size();
  const size_t np = aXYZ0.size() / 3;
  assert( aSkinSparseI.size() == aSkinSparseW.size() );
  assert( aSkinSparseI.size() % np == 0 );

  // --------

  const double eps = 1.0e-4;

  for(unsigned int itr=0;itr<3;++itr){ // Check Sensitivity Skin
    // set random joint rotation to Bone1
    std::vector<dfm2::CRigBone> aBone1 = aBone0;
    for (auto &bone : aBone1) {
      dfm2::CQuat<double>::Random(0.2).CopyTo(bone.quatRelativeRot);
//      dfm2::Quat_Identity(bone.quatRelativeRot);
    }
    UpdateBoneRotTrans(aBone1);
    dfm2::UpdateBoneRotTrans(aBone1);
    std::vector<double> aXYZ1;
    dfm2::SkinningSparse_LBS(aXYZ1,
        aXYZ0, aBone1, aSkinSparseW, aSkinSparseI);
    // compute Bon rotation sensitivity
    std::vector<double> L;  // [ nsns, nb*4 ]
    Rig_Sensitivity_Skeleton(
        L,
        aBone1);
    //
    const size_t nsns = L.size() / (nb * 12);
    assert(nsns == (nb + 1) * 3);
    for (int isns = 0; isns < nsns; ++isns) {
      std::vector<dfm2::CRigBone> aBone2 = aBone1;
      { // aBone2 is perturbed aBone1
        unsigned int ib_s = isns / 3;
        bool is_rot = true;
        unsigned int idim_s = isns - ib_s * 3;
        if (ib_s == nb) {
          ib_s = 0;
          is_rot = false;
        }
        if (is_rot) {
          dfm2::CQuatd dq = dfm2::Quat_CartesianAngle(eps * dfm2::CVec3d::Axis(idim_s));
          dfm2::CQuatd q0 = dq * dfm2::CQuatd(aBone2[ib_s].quatRelativeRot);
          q0.CopyTo(aBone2[ib_s].quatRelativeRot);
        } else {
          aBone2[ib_s].transRelative[idim_s] += eps;
        }
      }
      // ----------------
      dfm2::UpdateBoneRotTrans(aBone2);
      std::vector<double> aXYZ2;
      dfm2::SkinningSparse_LBS(aXYZ2,
          aXYZ0, aBone2, aSkinSparseW, aSkinSparseI);
      for (unsigned int ip = 0; ip < np; ++ip) {
        const double val0[3] = {
            (aXYZ2[ip * 3 + 0] - aXYZ1[ip * 3 + 0]) / eps,
            (aXYZ2[ip * 3 + 1] - aXYZ1[ip * 3 + 1]) / eps,
            (aXYZ2[ip * 3 + 2] - aXYZ1[ip * 3 + 2]) / eps};
        double val1[3] = {0, 0, 0};
        {
          std::vector<double> dP;
          Sensitivity_RigSkinPoint(dP,
              ip, aXYZ0.data()+ip*3, aBone1, L,
              static_cast<unsigned int>(aSkinSparseW.size()/np),
              aSkinSparseW, aSkinSparseI);
          val1[0] = dP[0*nsns+isns];
          val1[1] = dP[1*nsns+isns];
          val1[2] = dP[2*nsns+isns];
        }
        for (int i = 0; i < 3; ++i) {
          EXPECT_NEAR(val0[i], val1[i], 1.0e-3 * (fabs(val1[i]) + 1.0));
        }
      }
    }
  }
}


TEST(gltf, constraint_sensitivity )
{
  std::vector<dfm2::CRigBone> aBone;
  {
    std::string path_glb = std::string(PATH_INPUT_DIR) + "/CesiumMan.glb";
    dfm2::CGLTF gltf;
    gltf.Read(path_glb);
    //    gltf.Print();
    gltf.GetBone(aBone, 0);
    dfm2::SetCurrentBoneRotationAsDefault(aBone);
  }
  // ------------
  const size_t nb = aBone.size();
  std::vector<double> L;  // [ nsns, nb*4 ]
  Rig_Sensitivity_Skeleton(
      L,
      aBone);
  assert( L.size()==(nb+1)*3*nb*12 );

  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(+0, +1);
  std::uniform_int_distribution<int> rand_ib(+0, static_cast<int>(nb-1));
  const double eps = 1.0e-4;

  { // check sensitivity target
    for(auto & bone : aBone){
      dfm2::Quat_Identity(bone.quatRelativeRot);
    }
    dfm2::UpdateBoneRotTrans(aBone);
    // ----------------------
    std::vector<dfm2::CTarget> aTarget;
    for(int itr=0;itr<5;++itr){
      dfm2::CTarget t;
      t.ib = rand_ib(rndeng);
      t.pos = dfm2::CVec3d(dfm2::RandomVec3(dist_01,rndeng)) + aBone[t.ib].RootPosition();
      aTarget.push_back(t);
    }
    
    std::vector<double> aO0; // [nC]
    std::vector<double> adO0; // [nC, nb*3 ]
    for(auto & it : aTarget){
      dfm2::Rig_WdW_Target(aO0,adO0,
          aBone,it,L);
    }
    
    const size_t nsns = L.size()/(nb*12);
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
        dfm2::Rig_WdW_Target(aO1,adO1,
            aBone1,it,L);
      }
      // -------------
      for(int io=0;io<aO0.size();++io){
        EXPECT_NEAR(
            (aO1[io]-aO0[io])/eps,
            adO0[io*nsns+isns],
            0.5*eps*(fabs(adO0[io*nsns+isns])+1.0));
      }
    }
  }
}


TEST(gltf,set_default_rotation) {
  dfm2::CGLTF gltf;
  gltf.Read(std::string(PATH_INPUT_DIR) + "/CesiumMan.glb");
  // -----
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aTri;
  std::vector<double> aSkinSparseW;
  std::vector<unsigned int> aSkinSparseI;
  std::vector<dfm2::CRigBone> aBone;
  gltf.GetMeshInfo(
      aXYZ0, aTri, aSkinSparseW, aSkinSparseI,
      0,0);
  gltf.GetBone(aBone, 0);
  // there is the offset rotation in the rig bone
  std::vector<double> aXYZ1(aXYZ0.size());
  dfm2::Skinning_LBS_LocalWeight(
      aXYZ1.data(),
      aXYZ0.data(), aXYZ0.size()/3,
      aBone, aSkinSparseW.data(), aSkinSparseI.data());

  // remove offset rotation and make sure default mesh doesn't change
  dfm2::SetCurrentBoneRotationAsDefault(aBone);
  for(auto & bone : aBone){
    dfm2::Quat_Identity(bone.quatRelativeRot);
  }
  dfm2::UpdateBoneRotTrans(aBone);
  std::vector<double> aXYZ2(aXYZ0.size());
  dfm2::Skinning_LBS_LocalWeight(
      aXYZ2.data(),
      aXYZ0.data(), aXYZ0.size()/3,
      aBone, aSkinSparseW.data(), aSkinSparseI.data());
  //
  double d0 = 0.0, d1 = 0.0;
  for(unsigned int i=0;i<aXYZ0.size();++i){
    d0 += (aXYZ1[i]-aXYZ2[i])*(aXYZ1[i]-aXYZ2[i]);
    d1 += aXYZ1[i]*aXYZ1[i];
  }
  EXPECT_LE(d0,d1*1.0e-5);
}
