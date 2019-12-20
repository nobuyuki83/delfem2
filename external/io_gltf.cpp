#include <vector>
#include <string>
#include <iostream>
#include <map>

// Define these only in *one* .cc file.
#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
// #define TINYGLTF_NOEXCEPTION // optional. disable exception handling.
#include "delfem2/../../external/tinygltf/tiny_gltf.h"
#include "delfem2/../../external/io_gltf.h"


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


bool GetArray_UInt
(std::vector<unsigned int>& res,
 const tinygltf::Model& model, int iacc)
{
  const tinygltf::Accessor& acc = model.accessors[iacc];
  const int ibv = acc.bufferView;
  const tinygltf::BufferView& bv = model.bufferViews[ibv];
  const int ibuff = bv.buffer;
  const tinygltf::Buffer& buff = model.buffers[ibuff];
  const unsigned int ncnt = acc.count;
  unsigned int nelem = 0;
  if( acc.type == TINYGLTF_TYPE_SCALAR){ nelem = 1; }
  else if( acc.type == TINYGLTF_TYPE_VEC3 ){ nelem=3; }
  else if( acc.type == TINYGLTF_TYPE_VEC4 ){ nelem=4; }
  else{ std::cout << "Error!->unknown type: " << acc.type << std::endl; assert(0); abort(); }
  unsigned int ntot = ncnt*nelem;
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

bool GetArray_Double
(std::vector<double>& res,
 const tinygltf::Model& model, int iacc)
{
  const tinygltf::Accessor& acc = model.accessors[iacc];
  const int ibv = acc.bufferView;
  const tinygltf::BufferView& bv = model.bufferViews[ibv];
  const int ibuff = bv.buffer;
  const tinygltf::Buffer& buff = model.buffers[ibuff];
  const unsigned int ncnt = acc.count;
  unsigned int nelem = 0;
  if( acc.type == TINYGLTF_TYPE_SCALAR){ nelem = 1; }
  else if( acc.type == TINYGLTF_TYPE_VEC3 ){ nelem=3; }
  else if( acc.type == TINYGLTF_TYPE_VEC4 ){ nelem=4; }
  else if( acc.type == TINYGLTF_TYPE_MAT4 ){ nelem=16; }
  else{ std::cout << "Error!->unknown type: " << acc.type << std::endl; assert(0); abort(); }
  unsigned int ntot = ncnt*nelem;
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

void Print(const tinygltf::Model& model){
  for(size_t ib=0;ib<model.buffers.size();++ib){
    std::cout << "buffer: " << ib << " name: " << model.buffers[ib].name << std::endl;
    std::cout << "buffer: " << ib << " size: " << model.buffers[ib].data.size() << std::endl;
    //    const tinygltf::Buffer& buff = model.buffers[ib];
  }
  for(size_t ibv=0;ibv<model.bufferViews.size();++ibv){
    const tinygltf::BufferView& bv = model.bufferViews[ibv];
    std::cout << std::endl;
    std::cout << "buffer view " << ibv << " name: " << bv.name << std::endl;
    std::cout << "buffer view " << ibv << " offset: " << bv.byteOffset << std::endl;
    std::cout << "buffer view " << ibv << " length: " << bv.byteLength << std::endl;
    std::cout << "buffer view " << ibv << " buffer: " << bv.buffer << std::endl;
  }
  for(size_t iac=0;iac<model.accessors.size();++iac){
    const tinygltf::Accessor& ac = model.accessors[iac];
    std::cout << std::endl;
    std::cout << "accessor" << iac << " bufview: " << ac.bufferView << std::endl;
    std::cout << "accessor" << iac << " componentType: " << ac.componentType << std::endl;
    std::cout << "accessor" << iac << " type: " << ac.type << std::endl;
    std::cout << "accessor" << iac << " byteStride: " << ac.ByteStride(model.bufferViews[ac.bufferView]) << std::endl;
    std::cout << "accessor" << iac << " byteOffset: " << ac.byteOffset << std::endl;
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
    std::cout << "node: " << in << " name: " << node.name << std::endl;
    std::cout << "node: " << in << " iskin: " << node.skin << std::endl;
    std::cout << "node: " << in << " imesh: " << node.mesh << std::endl;
    std::cout << "node: " << in << " weightsize: " << node.weights.size() << std::endl;
    if( node.rotation.size() == 4 ){
      std::cout << "node: " << in << " rot:" << node.rotation[0] << " " << node.rotation[1] << " " << node.rotation[2] << " " << node.rotation[3] << std::endl;
    }
    if( node.translation.size() == 3 ){
      std::cout << "node: " << in << " trans:" << node.translation[0] << " " << node.translation[1] << " " << node.translation[2] << std::endl;
    }
    assert( node.scale.empty() || node.scale.size() == 3 );
    if( node.scale.size() == 3 ){
      std::cout << "node: " << in << " scale:" << node.scale[0] << " " << node.scale[1] << " " << node.scale[2] << std::endl;
    }
    assert( node.matrix.empty() || node.matrix.size() == 16 );
    if( node.matrix.size() == 16 ){
      for(int i=0;i<16;++i){
        std::cout << "    " << i << " " << node.matrix[i] << std::endl;
      }
    }
    std::cout << "node: " << in << " child: ";
    for(int ic : node.children){
      std::cout << ic << " ";
    }
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
    std::cout << "skin" << is << " inode_skeleton: " << skin.skeleton << std::endl;
    std::cout << "skin: " << is << " joints: ";
    for(int joint : skin.joints){
      std::cout << joint << " ";
    }
    std::cout << std::endl;
    std::cout << "skin" << is << " inverseBindMatrices: " << skin.inverseBindMatrices << std::endl;
    std::vector<double> M; GetArray_Double(M, model, skin.inverseBindMatrices);
    assert( M.size()%16 == 0 );
    assert( M.size()/16 == model.skins[is].joints.size() );
    const int nj = model.skins[is].joints.size();
    for(int ij=0;ij<nj;++ij){
      for(int i=0;i<16;++i){ std::cout << "   " << ij << " " << i << " " << M[ij*16+i] << std::endl; }
    }
  }
}

void GetMeshInfo
(std::vector<double>& aXYZ,
 std::vector<unsigned int>& aTri,
 std::vector<double>& aRigWeight,
 std::vector<unsigned int>& aRigJoint,
 const tinygltf::Model& model,
 int imsh, int iprimitive)
{
  aXYZ.clear();
  aTri.clear();
  aRigJoint.clear();
  aRigWeight.clear();
  //////
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
      std::cout << "has position: " << aXYZ.size()/3 << std::endl;
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
      std::cout << "has rig weight: " << aRigWeight.size()/4 << std::endl;
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
      std::cout << "has rig joint: " << aRigJoint.size()/4 << std::endl;
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

void GetBoneBinding
(std::vector<CRigBone>& aBone,
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
    for(int i=0;i<16;++i){ aBone[ij].Mat[i] = aBone[ij].invBindMat[i]; }
    int info; CalcInvMat(aBone[ij].Mat, 4, info);
  }
}

void SetBone
(std::vector<CRigBone>& aBone,
 const tinygltf::Model& model,
 int inode_cur, int ibone_p,
 const std::vector<int>& mapNode2Bone)
{
  assert( inode_cur >= 0 && inode_cur < (int)model.nodes.size() );
  const tinygltf::Node& node = model.nodes[inode_cur];
  const int ibone = mapNode2Bone[inode_cur];
  assert(ibone>=0&&ibone<(int)aBone.size());
  aBone[ibone].ibone_parent = ibone_p;
  aBone[ibone].trans[0] = node.translation[0];
  aBone[ibone].trans[1] = node.translation[1];
  aBone[ibone].trans[2] = node.translation[2];
  aBone[ibone].rot[0] = node.rotation[3];
  aBone[ibone].rot[1] = node.rotation[0];
  aBone[ibone].rot[2] = node.rotation[1];
  aBone[ibone].rot[3] = node.rotation[2];
  aBone[ibone].name = node.name;
  if( !node.scale.empty() ){ aBone[ibone].scale = node.scale[0]; }
  else{ aBone[ibone].scale = 1;  }
  // ----------------------------
  for(int inode_ch : node.children){
    SetBone(aBone, model, inode_ch, ibone, mapNode2Bone);
  }
}


bool CGLTF::Read(const std::string& fpath)
{
  std::string err;
  std::string warn;
  tinygltf::TinyGLTF loader;
  model = new tinygltf::Model;
  bool ret = loader.LoadBinaryFromFile(model, &err, &warn,
                                       fpath); // for binary glTF(.glb)
  if (!warn.empty()) { printf("Warn: %s\n", warn.c_str()); }
  if (!err.empty()) { printf("Err: %s\n", err.c_str()); }
  if (!ret) { printf("Failed to parse glTF\n"); return -1; }
  return true;
}



void CGLTF::Print() const
{
  ::Print(*model);
}


void CGLTF::GetMeshInfo
(std::vector<double>& aXYZ0,
 std::vector<unsigned int>& aTri,
 std::vector<double>& aRigWeight,
 std::vector<unsigned int>& aRigJoint,
 int imesh, int iprimitive) const
{
  ::GetMeshInfo(aXYZ0, aTri, aRigWeight, aRigJoint,
                *model, imesh, iprimitive);
}

void CGLTF::GetBone
(std::vector<CRigBone>& aBone,
 int iskin) const
{
  aBone.resize( model->skins[0].joints.size() );
  std::vector<int> mapNode2Bone( model->nodes.size(), -1);
  for(size_t ij=0;ij<model->skins[0].joints.size();++ij){
    int inode = model->skins[0].joints[ij];
    mapNode2Bone[inode] = ij;
  }
  SetBone(aBone,
          *model, model->skins[0].skeleton, -1, mapNode2Bone);
  GetBoneBinding(aBone,
                 *model);
}
