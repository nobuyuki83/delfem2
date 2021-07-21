#include <vector>
#include <iostream>
#include <map>

// Define these only in *one* .cc file.
#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
// #define TINYGLTF_NOEXCEPTION // optional. disable exception handling.
#include "tinygltf/tiny_gltf.h"
#include "io_gltf.h"

namespace dfm2 = delfem2;

// -----------------------------------------------------

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


bool GetArray_UInt(
    std::vector<unsigned int>& res,
    const tinygltf::Model& model,
    int iacc)
{
  const tinygltf::Accessor& acc = model.accessors[iacc];
  const int ibv = acc.bufferView;
  const tinygltf::BufferView& bv = model.bufferViews[ibv];
  const int ibuff = bv.buffer;
  const tinygltf::Buffer& buff = model.buffers[ibuff];
  const size_t ncnt = acc.count;
  unsigned int nelem = 0;
  if( acc.type == TINYGLTF_TYPE_SCALAR){ nelem = 1; }
  else if( acc.type == TINYGLTF_TYPE_VEC3 ){ nelem=3; }
  else if( acc.type == TINYGLTF_TYPE_VEC4 ){ nelem=4; }
  else{ std::cout << "Error!->unknown type: " << acc.type << std::endl; assert(0); abort(); }
  size_t ntot = ncnt*nelem;
  if( acc.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT ){ // unsigned short
    if( bv.byteStride != 0 && bv.byteStride != nelem*sizeof(unsigned short) ){
      std::cout << "Error!-->unsuppoted not packed" << std::endl;
      assert(0);
      abort();
    }
    assert( bv.byteLength >= ntot*sizeof(unsigned short) );
    const unsigned short* pdata = (unsigned short*)(buff.data.data() + bv.byteOffset + acc.byteOffset);
    res.assign(pdata,pdata+ntot);
    return true;
  }
  else{
    assert(0);
    abort();
  }
  return false;
}

bool GetArray_Double(
    std::vector<double>& res,
    const tinygltf::Model& model,
    int iacc)
{
  const tinygltf::Accessor& acc = model.accessors[iacc];
  const int ibv = acc.bufferView;
  const tinygltf::BufferView& bv = model.bufferViews[ibv];
  const int ibuff = bv.buffer;
  const tinygltf::Buffer& buff = model.buffers[ibuff];
  const size_t ncnt = acc.count;
  unsigned int nelem = 0;
  if( acc.type == TINYGLTF_TYPE_SCALAR){ nelem = 1; }
  else if( acc.type == TINYGLTF_TYPE_VEC3 ){ nelem=3; }
  else if( acc.type == TINYGLTF_TYPE_VEC4 ){ nelem=4; }
  else if( acc.type == TINYGLTF_TYPE_MAT4 ){ nelem=16; }
  else{ std::cout << "Error!->unknown type: " << acc.type << std::endl; assert(0); abort(); }
  size_t ntot = ncnt*nelem;
  if( acc.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT ){ // signed short
    if( bv.byteStride != 0 && bv.byteStride != nelem*sizeof(float) ){
      std::cout << "Error!-->unsuppoted not packed" << std::endl;
      assert(0);
      abort();
    }
    assert( bv.byteLength >= ntot*sizeof(float) );
    const float* pdata = (float*)(buff.data.data() + bv.byteOffset + acc.byteOffset);
    res.assign(pdata,pdata+ntot);
    return true;
  }
  return false;
}

void dfm2::Print(const tinygltf::Model& model){
  for(size_t ib=0;ib<model.buffers.size();++ib){
    std::cout << "buffer_idx: " << ib << " name: " << model.buffers[ib].name << std::endl;
    std::cout << "buffer_idx: " << ib << " uri: " << model.buffers[ib].uri << std::endl;
    std::cout << "buffer_idx: " << ib << " size: " << model.buffers[ib].data.size() << std::endl;
    //    const tinygltf::Buffer& buff = model.buffers[ib];
  }
  for(size_t ibv=0;ibv<model.bufferViews.size();++ibv){
    const tinygltf::BufferView& bv = model.bufferViews[ibv];
    std::cout << std::endl;
    std::cout << "bufferView_idx " << ibv << " name: " << bv.name << std::endl;
    std::cout << "bufferView_idx " << ibv << " offset: " << bv.byteOffset << std::endl;
    std::cout << "bufferView_idx " << ibv << " length: " << bv.byteLength << std::endl;
    std::cout << "bufferView_idx " << ibv << " buffer_idx: " << bv.buffer << std::endl;
  }
  for(size_t iac=0;iac<model.accessors.size();++iac){
    const tinygltf::Accessor& ac = model.accessors[iac];
    std::cout << std::endl;
    std::cout << "accessor_idx" << iac << " bufferView_idx: " << ac.bufferView << std::endl;
    std::cout << "accessor_idx" << iac << " componentType: " << ac.componentType << std::endl;
    std::cout << "accessor_idx" << iac << " type: " << ac.type << std::endl;
    std::cout << "accessor_idx" << iac << " byteStride: " << ac.ByteStride(model.bufferViews[ac.bufferView]) << std::endl;
    std::cout << "accessor_idx" << iac << " byteOffset: " << ac.byteOffset << std::endl;
  }
  for(size_t im=0;im<model.meshes.size();++im){
    const tinygltf::Mesh& mesh = model.meshes[im];
    for(size_t ipri=0;ipri<mesh.primitives.size();++ipri){
      std::cout << std::endl;
      const tinygltf::Primitive& primitive = mesh.primitives[ipri];
      for(const auto & attribute : primitive.attributes){
        std::cout << "mesh" << im << " primitive:" << ipri << " att: " << attribute.first << " acc:" << attribute.second << std::endl;
      }
      std::cout << "mesh" << im << " wieghtsize" << mesh.weights.size() << std::endl;
      {
        auto itr = primitive.attributes.find(std::string("POSITION"));
        if( itr != primitive.attributes.end() ){
          std::cout << "mesh" << im << " primitive" << ipri << " position acc: " << itr->second << std::endl;
        }
      }
      {
        auto itr = primitive.attributes.find(std::string("NORMAL"));
        if( itr != primitive.attributes.end() ){
          std::cout << "mesh" << im << " primitive" << ipri << " normal acc: " << itr->second << std::endl;
        }
      }
      std::cout << "mesh" << im << " primitive" << ipri << " index: " << primitive.indices << std::endl;
    }
  }
  for(size_t in=0;in<model.nodes.size();++in){
    const tinygltf::Node& node = model.nodes[in];
    std::cout << std::endl;
    std::cout << "node_idx: " << in << " name: " << node.name << std::endl;
    std::cout << "node_idx: " << in << " skin_idx: " << node.skin << std::endl;
    std::cout << "node_idx: " << in << " mesh_idx: " << node.mesh << std::endl;
    std::cout << "node_idx: " << in << " weightsize: " << node.weights.size() << std::endl;
    if( node.rotation.size() == 4 ){
      std::cout << "node_idx: " << in << " rot:" << node.rotation[0] << " " << node.rotation[1] << " " << node.rotation[2] << " " << node.rotation[3] << std::endl;
    }
    if( node.translation.size() == 3 ){
      std::cout << "node_idx: " << in << " trans:" << node.translation[0] << " " << node.translation[1] << " " << node.translation[2] << std::endl;
    }
    assert( node.scale.empty() || node.scale.size() == 3 );
    if( node.scale.size() == 3 ){
      std::cout << "node_idx: " << in << " scale:" << node.scale[0] << " " << node.scale[1] << " " << node.scale[2] << std::endl;
    }
    assert( node.matrix.empty() || node.matrix.size() == 16 );
    if( node.matrix.size() == 16 ){
      for(int i=0;i<16;++i){
        std::cout << "    " << i << " " << node.matrix[i] << std::endl;
      }
    }
    //
    std::cout << "node_idx: " << in << " child: ";
    for(int ic : node.children){ std::cout << ic << " "; }
    std::cout << std::endl;
  }
  /*
   for(int ia=0;ia<model.animations.size();++ia){
   const tinygltf::Animation& an = model.animations[ia];
   for(int ic=0;ic<an.channels.size();++ic){
   const tinygltf::AnimationChannel& ac = an.channels[ic];
   std::cout << "animation: " << ia << " " << "channel: " << ic << " sampler: " <<  ac.sampler << std::endl;
   std::cout << "animation: " << ia << " " << "channel: " << ic << " path: " <<  ac.target_path << std::endl;
   std::cout << "animation: " << ia << " " << "channel: " << ic << " node: " <<  ac.target_node << std::endl;
   }
   for(int is=0;is<an.samplers.size();++is){
   const tinygltf::AnimationSampler& as = an.samplers[is];
   std::cout << "animation: " << ia << " " << "sampler: " << is << " inut acc" << as.input << std::endl;
   std::cout << "animation: " << ia << " " << "sampler: " << is << " output" << as.output << std::endl;
   const tinygltf::Accessor& acc = model.accessors[as.input];
   std::cout << "   " << acc.minValues[0] << " " << acc.maxValues[0] << " " << acc.count << std::endl;
   std::vector<double> aTime; GetArray_Double(aTime, model, as.input);
   assert(aTime.size()==acc.count);
   std::vector<double> aVal; GetArray_Double(aVal, model, as.output);
   assert(aVal.size()%acc.count==0);
   //      std::cout << "      " << aVal.size() << std::endl;
   //      for(int ival=0;ival<aVal.size();++ival){
   //        std::cout << "    " << ival << " " << aVal[ival] << std::endl;
   /      }
   }
   }
   */
  for(size_t is=0;is<model.skins.size();++is){
    const tinygltf::Skin& skin = model.skins[is];
    std::cout << std::endl;
    std::cout << "skin_idx: " << is << "node_idx of bone root: " << skin.skeleton << std::endl;
    //
    std::cout << "skin_idx: " << is << " joint_idxs: ";
    for(int joint : skin.joints){
      std::cout << joint << " ";
    }
    std::cout << std::endl;
    //
    std::cout << "skin_idx" << is << " inverseBindMatrices: " << skin.inverseBindMatrices << std::endl;
    std::vector<double> M; GetArray_Double(M, model, skin.inverseBindMatrices);
    assert( M.size()%16 == 0 );
    assert( M.size()/16 == model.skins[is].joints.size() );
    const size_t nj = model.skins[is].joints.size();
    for(int ij=0;ij<nj;++ij){
      for(int i=0;i<16;++i){ std::cout << "   " << ij << " " << i << " " << M[ij*16+i] << std::endl; }
    }
  }
}

void dfm2::GetMeshInfo(
    std::vector<double>& aXYZ,
    std::vector<unsigned int>& aTri,
    std::vector<double>& aRigWeight,
    std::vector<unsigned int>& aRigJoint,
    const tinygltf::Model& model,
    int imsh,
    int iprimitive)
{
  aXYZ.clear();
  aTri.clear();
  aRigJoint.clear();
  aRigWeight.clear();
  //
  const tinygltf::Primitive& primitive = model.meshes[imsh].primitives[iprimitive];
  GetArray_UInt(aTri,
                model, primitive.indices);
  /*
   for(int it=0;it<aTri.size()/3;++it){
   std::cout << it << " --> " << aTri[it*3+0] << " " << aTri[it*3+1] << " " << aTri[it*3+2] << std::endl;
   }
   */
  {
    auto itr = primitive.attributes.find(std::string("POSITION"));
    if( itr != primitive.attributes.end() ){
      GetArray_Double(aXYZ,
                      model, itr->second);
//      std::cout << "has position: " << aXYZ.size()/3 << std::endl;
      /*
       for(int ip=0;ip<aXYZ.size()/3;++ip){
       std::cout << ip << " --> " << aXYZ[ip*3+0] << " " << aXYZ[ip*3+1] << " " << aXYZ[ip*3+2] << std::endl;
       }
       */
    }
  }
  {
    auto itr = primitive.attributes.find(std::string("NORMAL"));
    if( itr != primitive.attributes.end() ){
    }
  }
  {
    auto itr = primitive.attributes.find(std::string("WEIGHTS_0"));
    if( itr != primitive.attributes.end() ){
      GetArray_Double(aRigWeight,
                      model, itr->second);
//      std::cout << "has rig weight: " << aRigWeight.size()/4 << std::endl;
      //      assert( aRigWeight.size()/4 == aXYZ.size()/3 );
      /*
       for(int ir=0;ir<aRigWeight.size()/4;++ir){
       std::cout << ir << " ";
       std::cout << aRigWeight[ir*2+0] << " ";
       std::cout << aRigWeight[ir*2+1] << " ";
       std::cout << aRigWeight[ir*2+2] << " ";
       std::cout << aRigWeight[ir*2+3] << std::endl;
       }
       */
    }
  }
  {
    auto itr = primitive.attributes.find(std::string("JOINTS_0"));
    if( itr != primitive.attributes.end() ){
      GetArray_UInt(aRigJoint,
                    model, itr->second);
//      std::cout << "has rig joint: " << aRigJoint.size()/4 << std::endl;
      /*
       for(int ir=0;ir<aRigJoint.size()/4;++ir){
       std::cout << ir << " --> ";
       std::cout << aRigJoint[ir*4+0] << " ";
       std::cout << aRigJoint[ir*4+1] << " ";
       std::cout << aRigJoint[ir*4+2] << " ";
       std::cout << aRigJoint[ir*4+3] << std::endl;
       }
       //      assert( aRigJoint.size()/4 == aXYZ.size()/3 );
       */
    }
  }
}

void dfm2::GetBoneBinding(
    std::vector<dfm2::CRigBone>& aBone,
    const tinygltf::Model& model)
{
  const tinygltf::Skin& skin = model.skins[0];
  std::vector<double> M; GetArray_Double(M, model, skin.inverseBindMatrices);
  assert( M.size() == aBone.size()*16 );
  for(size_t ij=0;ij<M.size()/16;++ij){
    for(int i=0;i<4;++i){
      for(int j=0;j<4;++j){
        aBone[ij].invBindMat[i*4+j] = M[ij*16+j*4+i];
      }
    }
    for(int i=0;i<16;++i){ aBone[ij].affmat3Global[i] = aBone[ij].invBindMat[i]; }
    int info; CalcInvMat(aBone[ij].affmat3Global, 4, info);
  }
}

void dfm2::SetBone(
    std::vector<dfm2::CRigBone>& aBone,
    const tinygltf::Model& model,
    unsigned int inode0,
    int ibone_p,
    const std::vector<unsigned int>& mapNode2Bone)
{
  assert(inode0 < model.nodes.size() );
  const tinygltf::Node& node = model.nodes[inode0];
  const unsigned int ibone0 = mapNode2Bone[inode0];
  assert( ibone0 < aBone.size() );
  aBone[ibone0].ibone_parent = ibone_p;
  if( node.translation.size() == 3 ) { // there is a chance node.translation is omitted
    aBone[ibone0].transRelative[0] = node.translation[0];
    aBone[ibone0].transRelative[1] = node.translation[1];
    aBone[ibone0].transRelative[2] = node.translation[2];
  }
  if( node.rotation.size() == 4 ) { // there is a chance node.rotation is omitted
    aBone[ibone0].quatRelativeRot[0] = node.rotation[0];
    aBone[ibone0].quatRelativeRot[1] = node.rotation[1];
    aBone[ibone0].quatRelativeRot[2] = node.rotation[2];
    aBone[ibone0].quatRelativeRot[3] = node.rotation[3];
  }
  if( !node.scale.empty() ){
    aBone[ibone0].scale = node.scale[0];
  }
  aBone[ibone0].name = node.name;
  // ----------------------------
  for(int inode_ch : node.children){
    SetBone(aBone,
        model, inode_ch, ibone0, mapNode2Bone);
  }
}

void delfem2::GetBone(
    std::vector<dfm2::CRigBone>& aBone,
    const tinygltf::Model& model,
    unsigned int iskin)
{
  assert( iskin < model.skins.size() );
  aBone.resize( model.skins[iskin].joints.size() );
  std::vector<unsigned int> mapNode2Bone( model.nodes.size(), UINT_MAX);
  for(unsigned int ij=0;ij<model.skins[iskin].joints.size();++ij){
    const unsigned int inode = model.skins[iskin].joints[ij];
    assert( inode < model.nodes.size() );
    mapNode2Bone[inode] = ij;
  }
  unsigned int inode_root = model.skins[iskin].skeleton;
  if( inode_root == UINT_MAX && !model.skins[iskin].joints.empty() ){
    inode_root = model.skins[iskin].joints[0];
  }
  assert( inode_root < model.nodes.size() );

  SetBone(
      aBone,
      model, inode_root, -1, mapNode2Bone);
  GetBoneBinding(
      aBone,
      model);
}


// --------------------------------
// implementation of GLTF class

bool dfm2::CGLTF::Read(
    const std::string& fpath)
{
  std::string err;
  std::string warn;
  tinygltf::TinyGLTF loader;
  pModel = new tinygltf::Model;
  bool ret = loader.LoadBinaryFromFile(pModel, &err, &warn,
                                       fpath); // for binary glTF(.glb)
  if (!warn.empty()) { printf("Warn: %s\n", warn.c_str()); }
  if (!err.empty()) { printf("Err: %s\n", err.c_str()); }
  if (!ret) { printf("Failed to parse glTF\n"); return false; }
  return true;
}

void dfm2::CGLTF::Print() const
{
  dfm2::Print(*pModel);
}

void dfm2::CGLTF::GetMeshInfo(
    std::vector<double>& aXYZ0,
    std::vector<unsigned int>& aTri,
    std::vector<double>& aRigWeight,
    std::vector<unsigned int>& aRigJoint,
    int imesh,
    int iprimitive) const
{
  dfm2::GetMeshInfo(
      aXYZ0, aTri, aRigWeight, aRigJoint,
      *pModel, imesh, iprimitive);
}

void dfm2::CGLTF::GetBone(
    std::vector<dfm2::CRigBone>& aBone,
    unsigned int iskin) const
{
  ::delfem2::GetBone(aBone,*pModel,iskin);
}
