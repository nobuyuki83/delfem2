/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <vector>
#include "delfem2/noise.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_funcs.h"

// -----------------------------
std::vector<double> aXYZ;
std::vector<unsigned int> aTri;
std::vector<int> aP;
std::vector<double> aGrad;
int nH, nW, nD;
std::vector<unsigned char> aV;
// -------------------------------

void myGlutDisplay()
{
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_TEXTURE_3D);
  glEnable(GL_TEXTURE_GEN_S);
  glEnable(GL_TEXTURE_GEN_T);
  glEnable(GL_TEXTURE_GEN_R);
  delfem2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aTri);
//  ::glutSolidTeapot(1.0);
  //    glutSolidSphere(1.0, 32, 16);
  //  glutSolidDodecahedron();
  ::glDisable(GL_TEXTURE_3D);
  glDisable(GL_TEXTURE_GEN_S);
  glDisable(GL_TEXTURE_GEN_T);
  glDisable(GL_TEXTURE_GEN_R);
}

void ComputePerlin(){
  
  aP.resize(256);
  for(int i=0;i<256;++i){ aP[i]=i; }
  delfem2::Shuffle(aP);
  aP.resize(512);
  for(int i=0;i<256;++i){ aP[256+i]=i; }
  
  aGrad.clear();
  
  /*
   for(int i=0;i<12;++i){
   double x = (rand()/(RAND_MAX+1.0))-0.5;
   double y = (rand()/(RAND_MAX+1.0))-0.5;
   double len_inv = 1.0/sqrt(x*x+y*y);
   x *= len_inv;
   y *= len_inv;
   aGrad.push_back(x);
   aGrad.push_back(y);
   }
   */
  
  aGrad.push_back(-1); aGrad.push_back(-1); aGrad.push_back(+0);
  aGrad.push_back(-1); aGrad.push_back(+1); aGrad.push_back(+0);
  aGrad.push_back(+1); aGrad.push_back(-1); aGrad.push_back(+0);
  aGrad.push_back(+1); aGrad.push_back(+1); aGrad.push_back(+0);
  aGrad.push_back(+0); aGrad.push_back(-1); aGrad.push_back(-1);
  aGrad.push_back(+0); aGrad.push_back(-1); aGrad.push_back(+1);
  aGrad.push_back(+0); aGrad.push_back(+1); aGrad.push_back(-1);
  aGrad.push_back(+0); aGrad.push_back(+1); aGrad.push_back(+1);
  aGrad.push_back(-1); aGrad.push_back(+0); aGrad.push_back(-1);
  aGrad.push_back(-1); aGrad.push_back(+0); aGrad.push_back(+1);
  aGrad.push_back(+1); aGrad.push_back(+0); aGrad.push_back(-1);
  aGrad.push_back(+1); aGrad.push_back(+0); aGrad.push_back(+1);
  
  
  nH = 128;
  nW = 128;
  nD = 128;
  aV.resize(nH*nW*nD*4);
  int nrep = 4;
  for(int id=0;id<nD;++id){
    for(int ih=0;ih<nH;++ih){
      for(int iw=0;iw<nW;++iw){
        double x = (double)iw/nH*nrep;
        double y = (double)ih/nW*nrep;
        double z = (double)id/nD*nrep;
        double v = delfem2::noise_perlin_3d_oct(x,y,z,nrep, 4,0.5, aGrad,aP);
        //        double v = noise_perlin_3d(x,y,z, aGrad,aP);
        double v0 = v*128+128;
        if( v0 < 0   ){ v0 =   0; }
        if( v0 > 255 ){ v0 = 255; }
        auto ucv = (unsigned char)v0;
        aV[(id*nW*nH+ih*nW+iw)*4+0] = ucv;
        aV[(id*nW*nH+ih*nW+iw)*4+1] = ucv;
        aV[(id*nW*nH+ih*nW+iw)*4+2] = ucv;
        aV[(id*nW*nH+ih*nW+iw)*4+3] = 255;
      }
    }
  }
  
  /*
   Noise3 noise(5, 5, 5);
   for (int id = 0; id < nD; ++id) {
   double p = (double)id / (double)nD;
   for (int ih = 0; ih < nH; ++ih) {
   double q = (double)ih / (double)nH;
   for (int iw = 0; iw < nW; ++iw) {
   double r = (double)iw / (double)nW;
   GLubyte ub = (GLubyte)(noise.noise(p, q, r) * 127.5 + 127.5);
   //        texture[k][j][i][0] = texture[k][j][i][1] = texture[k][j][i][2] = ub;
   //        texture[k][j][i][3] = 255;
   aV[(id*nW*nH+ih*nW+iw)*4+0] = ub;
   aV[(id*nW*nH+ih*nW+iw)*4+1] = ub;
   aV[(id*nW*nH+ih*nW+iw)*4+2] = ub;
   aV[(id*nW*nH+ih*nW+iw)*4+3] = 255;
   }
   }
   }
   */
}


int main(int argc,char* argv[])
{
  delfem2::Read_Ply(std::string(PATH_INPUT_DIR)+"/bunny_1k.ply",
           aXYZ,aTri);
  delfem2::Normalize_Points3(aXYZ);
  
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  
  ComputePerlin();
  
  /*
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
   */
  
  printf("OpenGL version supported by this platform (%s): \n",
         glGetString(GL_VERSION));
  
  glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
  glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
  glTexGeni(GL_R, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
  static double genfunc[][4] = {
    { 1.0, 0.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0, 0.0 },
    { 0.0, 0.0, 1.0, 0.0 },
  };
  glTexGendv(GL_S, GL_OBJECT_PLANE, genfunc[0]);
  glTexGendv(GL_T, GL_OBJECT_PLANE, genfunc[1]);
  glTexGendv(GL_R, GL_OBJECT_PLANE, genfunc[2]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, nW, nH, nD, 0,
               GL_RGBA, GL_UNSIGNED_BYTE, aV.data());
  
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  
  delfem2::opengl::setSomeLighting();
  

  
  while (!glfwWindowShouldClose(viewer.window))
  {
    viewer.DrawBegin_oldGL();
    
    myGlutDisplay();
    
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
