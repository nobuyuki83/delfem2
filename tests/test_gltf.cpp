/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>

#include "gtest/gtest.h"


// Define these only in *one* .cc file.
#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
// #define TINYGLTF_NOEXCEPTION // optional. disable exception handling.
#include "delfem2/../../external/tinygltf/tiny_gltf.h"

TEST(gltf,formatcheck)
{
  tinygltf::TinyGLTF loader;
  for(int ifile=0;ifile<5;++ifile){
    std::string path_glb = std::string(PATH_INPUT_DIR);
    if(      ifile == 0 ){ path_glb += "/Duck.glb"; }
    else if( ifile == 1 ){ path_glb += "/RiggedSimple.glb"; }
    else if( ifile == 2 ){ path_glb += "/RiggedFigure.glb"; }
    else if( ifile == 3 ){ path_glb += "/Monster.glb"; }
    else if( ifile == 4 ){ path_glb += "/CesiumMan.glb"; }
    tinygltf::Model model;
    std::string err;
    std::string warn;
    bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, path_glb); // for binary glTF(.glb)
    EXPECT_TRUE(warn.empty());
    EXPECT_TRUE(err.empty());
    EXPECT_TRUE(ret);
    for(int in=0;in<model.nodes.size();++in){
      const tinygltf::Node& node = model.nodes[in];
      EXPECT_TRUE( node.rotation.size() == 0 || node.rotation.size() == 4 );
      EXPECT_TRUE( node.translation.size() == 0 || node.translation.size() == 3 );
      EXPECT_TRUE( node.scale.size() == 0 || node.scale.size() == 3 );
      EXPECT_TRUE( node.matrix.size() == 0 || node.matrix.size() == 16 );
      EXPECT_FALSE( node.matrix.size() > 0 && node.scale.size() > 0 );      // if there is matrix, no scale
      EXPECT_FALSE( node.matrix.size() > 0 && node.translation.size() > 0 );     // if there is matrix, no translation
      EXPECT_FALSE( node.matrix.size() > 0 && node.rotation.size() > 0 );     // if there is matrix, no rotation
      if( node.skin != -1 ){
        EXPECT_TRUE( node.skin>=0 && node.skin<model.skins.size() );
      }
      if( node.mesh != -1 ){
        EXPECT_TRUE( node.mesh>=0 && node.mesh<model.meshes.size() );
      }
      for(int ic=0;ic<node.children.size();++ic){
        EXPECT_TRUE( node.children[ic]>=0 && node.children[ic]<model.nodes.size() );
      }
    }
    for(int is=0;is<model.skins.size();++is){
      const tinygltf::Skin& skin = model.skins[is];
      for(int ij=0;ij<skin.joints.size();++ij){
        int inode = skin.joints[ij];
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

