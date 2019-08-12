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
#include "delfem2/mshio.h"
#include "delfem2/funcs.h"

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


CGlutWindowManager window;
std::vector<double> aXYZ;
std::vector<unsigned int> aTri;

bool is_animation = false;

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
  
  DrawAxis(10);
  
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
  }
  else{
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

bool GetArray_UInt
(std::vector<unsigned int>& res,
 const tinygltf::Model& model, int iacc)
{
  const tinygltf::Accessor& acc = model.accessors[iacc];
  const int ibv = acc.bufferView;
  const tinygltf::BufferView& bv = model.bufferViews[ibv];
  const int ibuff = bv.buffer;
  const tinygltf::Buffer& buff = model.buffers[ibuff];
  int buffoffset = bv.byteOffset;
  int buffLength = bv.byteLength;
  if( acc.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT ){ // signed short
    const unsigned short* pdata = (unsigned short*)(buff.data.data() + buffoffset);
    res.assign(pdata,pdata+buffLength/sizeof(unsigned short));
    return true;
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
  const int buffLength = bv.byteLength;
  if( acc.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT ){ // signed short
    const float* pdata = (float*)(buff.data.data() + bv.byteOffset + acc.byteOffset);
    res.assign(pdata,pdata+buffLength/sizeof(float));
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
          GetArray_Double(aXYZ,
                          model, itr->second);
        }
      }
      {
        auto itr = primitive.attributes.find(std::string("NORMAL"));
        if( itr != primitive.attributes.end() ){
          std::cout << "mesh" << im << " primitive" << ipri << " normal acc: " << itr->second << std::endl;
        }
      }
      std::cout << "mesh" << im << " primitive" << ipri << " index: " << primitive.indices << std::endl;
      GetArray_UInt(aTri,
                    model, primitive.indices);
    }
  }
  for(int is=0;is<model.skins.size();++is){
    const tinygltf::Skin& skin = model.skins[is];
    std::cout << "skin" << is << " inode_skeleton: " << skin.skeleton << std::endl;
    std::cout << "skin: " << is << " joints: ";
    for(int ij=0;ij<skin.joints.size();++ij){
      std::cout << skin.joints[ij] << " ";
    }
    std::cout << std::endl;
  }
  for(int in=0;in<model.nodes.size();++in){
    const tinygltf::Node& node = model.nodes[in];
    std::cout << std::endl;
    std::cout << "node: " << in << " name: " << node.name << std::endl;
    std::cout << "node: " << in << " iskin: " << node.skin << std::endl;
    std::cout << "node: " << in << " imesh: " << node.mesh << std::endl;
    std::cout << "node: " << in << " weightsize: " << node.weights.size() << std::endl;
    std::cout << "node: " << in << " child: ";
    for(int ic=0;ic<node.children.size();++ic){
      std::cout << node.children[ic] << " ";
    }
    std::cout << std::endl;
  }
  for(int ia=0;ia<model.animations.size();++ia){
    const tinygltf::Animation& an = model.animations[ia];
    for(int ic=0;ic<an.channels.size();++ic){
      const tinygltf::AnimationChannel& ac = an.channels[ic];
      std::cout << "animation: " << ia << " " << "channel: " << ic << " sampler: " <<  ac.sampler << std::endl;
      std::cout << "animation: " << ia << " " << "channel: " << ic << " sampler: " <<  ac.target_path << std::endl;
    }
    for(int is=0;is<an.samplers.size();++is){
      const tinygltf::AnimationSampler& as = an.samplers[is];
      std::cout << "animation: " << ia << " " << "sampler: " << is << " " << as.input << " " << as.output << std::endl;
    }
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
    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedSimple.glb";
//    std::string path_gltf = std::string(PATH_INPUT_DIR)+"/RiggedFigure.glb";
    bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, path_gltf); // for binary glTF(.glb)
    
    if (!warn.empty()) {
      printf("Warn: %s\n", warn.c_str());
    }
    
    if (!err.empty()) {
      printf("Err: %s\n", err.c_str());
    }
    
    if (!ret) {
      printf("Failed to parse glTF\n");
      return -1;
    }
    Print(model);
    {
      const tinygltf::Primitive& primitive = model.meshes[0].primitives[0];
      GetArray_UInt(aTri,
                    model, primitive.indices);
      {
        auto itr = primitive.attributes.find(std::string("POSITION"));
        if( itr != primitive.attributes.end() ){
          GetArray_Double(aXYZ,
                          model, itr->second);
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
          std::cout << "has rig weight" << std::endl;
        }
      }
      {
        auto itr = primitive.attributes.find(std::string("JOINTS_0"));
        if( itr != primitive.attributes.end() ){
          std::cout << "has rig joint" << std::endl;
        }
      }

    }

  }
  
  window.camera.view_height = 150.0;
  window.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  setSomeLighting();
  
  glutMainLoop();
  return 0;
}
