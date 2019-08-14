#include <iostream>
#include <sstream>
#include <vector>
#include <set>

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/vec3.h"
#include "delfem2/vec2.h"
#include "delfem2/mat3.h"
#include "delfem2/mshio.h"
#include "delfem2/funcs.h"
#include "delfem2/quat.h"

#include "delfem2/gl_funcs.h"
#include "delfem2/gl_v23q.h"
#include "delfem2/glut_funcs.h"

#include "delfem2/rigmesh.h"

// Define these only in *one* .cc file.
#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
// #define TINYGLTF_NOEXCEPTION // optional. disable exception handling.
#include "delfem2/../../external/tinygltf/tiny_gltf.h"


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
  for(int ib=0;ib<model.buffers.size();++ib){
    std::cout << "buffer: " << ib << " name: " << model.buffers[ib].name << std::endl;
    std::cout << "buffer: " << ib << " size: " << model.buffers[ib].data.size() << std::endl;
    //    const tinygltf::Buffer& buff = model.buffers[ib];
  }
  for(int ibv=0;ibv<model.bufferViews.size();++ibv){
    const tinygltf::BufferView& bv = model.bufferViews[ibv];
    std::cout << std::endl;
    std::cout << "buffer view " << ibv << " name: " << bv.name << std::endl;
    std::cout << "buffer view " << ibv << " offset: " << bv.byteOffset << std::endl;
    std::cout << "buffer view " << ibv << " length: " << bv.byteLength << std::endl;
    std::cout << "buffer view " << ibv << " buffer: " << bv.buffer << std::endl;
  }
  for(int iac=0;iac<model.accessors.size();++iac){
    const tinygltf::Accessor& ac = model.accessors[iac];
    std::cout << std::endl;
    std::cout << "accessor" << iac << " bufview: " << ac.bufferView << std::endl;
    std::cout << "accessor" << iac << " componentType: " << ac.componentType << std::endl;
    std::cout << "accessor" << iac << " type: " << ac.type << std::endl;
    std::cout << "accessor" << iac << " byteStride: " << ac.ByteStride(model.bufferViews[ac.bufferView]) << std::endl;
    std::cout << "accessor" << iac << " byteOffset: " << ac.byteOffset << std::endl;
  }
  for(int im=0;im<model.meshes.size();++im){
    const tinygltf::Mesh& mesh = model.meshes[im];
    for(int ipri=0;ipri<mesh.primitives.size();++ipri){
      std::cout << std::endl;
      const tinygltf::Primitive& primitive = mesh.primitives[ipri];
      for(auto itr = primitive.attributes.begin();itr!=primitive.attributes.end();++itr){
        std::cout << "mesh" << im << " primitive:" << ipri << " att: " << itr->first << " acc:" << itr->second << std::endl;
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
  for(int in=0;in<model.nodes.size();++in){
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
    assert( node.scale.size() == 0 || node.scale.size() == 3 );
    if( node.scale.size() == 3 ){
      std::cout << "node: " << in << " scale:" << node.scale[0] << " " << node.scale[1] << " " << node.scale[2] << std::endl;
    }
    assert( node.matrix.size() == 0 || node.matrix.size() == 16 );
    if( node.matrix.size() == 16 ){
      for(int i=0;i<16;++i){
        std::cout << "    " << i << " " << node.matrix[i] << std::endl;
      }
    }
    std::cout << "node: " << in << " child: ";
    for(int ic=0;ic<node.children.size();++ic){
      std::cout << node.children[ic] << " ";
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
  for(int is=0;is<model.skins.size();++is){
    const tinygltf::Skin& skin = model.skins[is];
    std::cout << std::endl;
    std::cout << "skin" << is << " inode_skeleton: " << skin.skeleton << std::endl;
    std::cout << "skin: " << is << " joints: ";
    for(int ij=0;ij<skin.joints.size();++ij){
      std::cout << skin.joints[ij] << " ";
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






//////////////////////////////////////////////////////////////////////

CGlutWindowManager window;
std::vector<double> aXYZ0;
std::vector<unsigned int> aTri;
std::vector<double> aRigWeight;
std::vector<unsigned int> aRigJoint;
////
std::vector<double> aXYZ;
std::vector<CRigBone> aBone;
int ibone_selected = -1;
int ielem_bone_selected = -1;
double rad_bone_sphere = 0.01;
double rad_rot_hndlr = 1.0;

bool is_animation = false;

void UpdateRigSkin()
{
  const int np = aXYZ0.size()/3;
  assert(aRigWeight.size()==np*4);
  assert(aRigJoint.size()==np*4);
  aXYZ.resize(aXYZ0.size());
  for(int ip=0;ip<np;++ip){
    double pos0[4] = {aXYZ0[ip*3+0],aXYZ0[ip*3+1],aXYZ0[ip*3+2],1.0};
    double pos1[3] = {0,0,0};
    double sum_w = 0.0;
    for(int iij=0;iij<4;++iij){
      double w = aRigWeight[ip*4+iij];
      if( w < 1.0e-30 ){ continue; }
      int ij = aRigJoint[ip*4+iij];
      sum_w += w;
      assert (ij>=0 && ij<aBone.size());
      double pos0a[4]; MatVec4(pos0a,aBone[ij].invBindMat,pos0);
      double pos0b[4]; MatVec4(pos0b,aBone[ij].Mat,pos0a);
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

void myGlutDisplay(void)
{
  //  ::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(0.5f, 0.8f, 1.0f ,1.0f);
  ::glClearStencil(0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  window.SetGL_Camera();
  
//  ::glDisable(GL_LIGHTING);
//  DrawMeshTri3D_Edge(aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
  
  ::glEnable(GL_LIGHTING);
  DrawMeshTri3D_FaceNorm(aXYZ.data(), aTri.data(), aTri.size()/3);
  
  DrawAxis(1);
  
  
  ::glDisable(GL_DEPTH_TEST);
  DrawBone(aBone,
           ibone_selected, ielem_bone_selected,
           rad_bone_sphere, rad_rot_hndlr);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glColor3d(0,0,0);
  ShowFPS();
  ::glutSwapBuffers();
}


void myGlutIdle()
{
  if( is_animation ){
  }
  ::glutPostRedisplay();
}


void myGlutResize(int w, int h)
{
  ::glViewport(0,0,w,h);
  ::glutPostRedisplay();
}

void myGlutSpecial(int Key, int x, int y)
{
  window.glutSpecial(Key, x, y);
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y )
{
  window.glutMotion(x, y);
  if( window.imodifier != 0 ) return;
  ////
  if( ibone_selected>=0 && ibone_selected<aBone.size() ){
    window.SetGL_Camera();
    float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
    float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
    CVector2 sp1(window.mouse_x, window.mouse_y);
    CVector2 sp0(window.mouse_x-window.dx, window.mouse_y-window.dy);
    CRigBone& bone = aBone[ibone_selected];
    DragHandlerRot_Mat4(bone.rot,
                        ielem_bone_selected, sp0, sp1, bone.Mat,
                        mMV, mPj);
    UpdateBoneRotTrans(aBone);
    UpdateRigSkin();
  }
  ////
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  window.glutMouse(button, state, x, y);
  /////
  window.SetGL_Camera();
  float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
  float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
  CVector3 src = screenUnProjection(CVector3(window.mouse_x,window.mouse_y,0.0), mMV, mPj);
  CVector3 dir = screenDepthDirection(src,mMV,mPj);
  if( state == GLUT_DOWN ){
    const double wh = 1.0/mPj[5];
    std::cout << wh << std::endl;
    PickBone(ibone_selected, ielem_bone_selected,
             aBone,
             src,dir,
             rad_rot_hndlr,
             wh*0.05);
    std::cout << ibone_selected << std::endl;
  }
  /////
  ::glutPostRedisplay();
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
  {
    case 'q':
    case 'Q':
    case '\033':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
    case 'a':
    {
      is_animation = !is_animation;
      break;
    }
  }
  ::glutPostRedisplay();
}



void SetBone
(std::vector<CRigBone>& aBone,
 const tinygltf::Model& model,
 int inode_cur, int ibone_p,
 const std::vector<int>& mapNode2Bone)
{
  assert( inode_cur >= 0 && inode_cur < model.nodes.size() );
  const tinygltf::Node& node = model.nodes[inode_cur];
  const int ibone = mapNode2Bone[inode_cur];
  assert(ibone>=0&&ibone<aBone.size());
  aBone[ibone].ibone_parent = ibone_p;
  aBone[ibone].trans[0] = node.translation[0];
  aBone[ibone].trans[1] = node.translation[1];
  aBone[ibone].trans[2] = node.translation[2];
  aBone[ibone].rot[0] = node.rotation[3];
  aBone[ibone].rot[1] = node.rotation[0];
  aBone[ibone].rot[2] = node.rotation[1];
  aBone[ibone].rot[3] = node.rotation[2];
  aBone[ibone].name = node.name;
  if( node.scale.size() > 0 ){ aBone[ibone].scale = node.scale[0]; }
  else{ aBone[ibone].scale = 1;  }
  ////
  for(int ich=0;ich<node.children.size();++ich){
    int inode_ch = node.children[ich];
    SetBone(aBone, model, inode_ch, ibone, mapNode2Bone);
  }
}

int main(int argc,char* argv[])
{
  glutInit(&argc, argv);
  
  // Initialize GLUT window 3D
  glutInitWindowPosition(200,200);
  glutInitWindowSize(400, 300);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_STENCIL);
  glutCreateWindow("3D View");
  glutDisplayFunc(myGlutDisplay);
  glutIdleFunc(myGlutIdle);
  glutReshapeFunc(myGlutResize);
  glutMotionFunc(myGlutMotion);
  glutMouseFunc(myGlutMouse);
  glutKeyboardFunc(myGlutKeyboard);
  glutSpecialFunc(myGlutSpecial);
  
  ////////////////////////
  
  
  {
    tinygltf::Model model;
    tinygltf::TinyGLTF loader;
    std::string err;
    std::string warn;
    
    //std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedSimple.gltf";
    //bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, path_gltf);
    
//    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Duck.glb";
//      std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedSimple.glb";
//    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedFigure.glb";
//    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/Monster.glb";
    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/CesiumMan.glb";
  
    bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, path_gltf); // for binary glTF(.glb)
    if (!warn.empty()) { printf("Warn: %s\n", warn.c_str()); }
    if (!err.empty()) { printf("Err: %s\n", err.c_str()); }
    if (!ret) { printf("Failed to parse glTF\n"); return -1; }
    
    Print(model);
    ////
    GetMeshInfo(aXYZ0, aTri, aRigWeight, aRigJoint,
                model, 0, 0);
    if( !model.skins.empty() ){
      aBone.resize( model.skins[0].joints.size() );
      std::vector<int> mapNode2Bone( model.nodes.size(), -1);
      for(int ij=0;ij<model.skins[0].joints.size();++ij){
        int inode = model.skins[0].joints[ij];
        mapNode2Bone[inode] = ij;
      }
      SetBone(aBone,
              model, model.skins[0].skeleton, -1, mapNode2Bone);
      const tinygltf::Skin& skin = model.skins[0];
      std::vector<double> M; GetArray_Double(M, model, skin.inverseBindMatrices);
      assert( M.size() == aBone.size()*16 );
      for(int ij=0;ij<M.size()/16;++ij){
        for(int i=0;i<4;++i){
          for(int j=0;j<4;++j){
            aBone[ij].invBindMat[i*4+j] = M[ij*16+j*4+i];
          }
        }
        for(int i=0;i<16;++i){ aBone[ij].Mat[i] = aBone[ij].invBindMat[i]; }
        int info; CalcInvMat(aBone[ij].Mat, 4, info);
      }
    }
    UpdateBoneRotTrans(aBone);
    UpdateRigSkin();
  }
  
  window.camera.view_height = 2.0;
  window.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  setSomeLighting();
  
  glutMainLoop();
  return 0;
}
