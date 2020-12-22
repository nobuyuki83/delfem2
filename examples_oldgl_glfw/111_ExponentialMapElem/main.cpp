/*
 * Copyright (c) 2019-2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/tex.h"
#include "delfem2/expmap_geo3dijk.h"
#include "delfem2/points.h"
#include "delfem2/mshio.h"
#include "delfem2/mshuni.h"
#include "delfem2/vec3.h"
#include "delfem2/imgio.h"
#include <GLFW/glfw3.h>
#include <vector>
#include <string>
#include <cstdlib>
#include <set>

namespace dfm2 = delfem2;

// ---------------------------

void Draw(
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& aTri,
    const std::vector<double>& aTexP,
    const std::vector<unsigned int>& elsup_ind,
    const std::vector<unsigned int>& elsup)
{
  ::glEnable(GL_TEXTURE_2D);
  ::glEnable(GL_LIGHTING);
  float gray[4] = {1,1,1,1};
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, gray);
  float shine[4] = {1,1,1,1};
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
  ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128.0);
  ::glBegin(GL_TRIANGLES);
  //
  for(unsigned int it0=0;it0<aTri.size()/3;++it0){
    {
      dfm2::CVec3d n0 = dfm2::Normal_Tri3(it0, aTri, aXYZ);
      ::glNormal3dv(n0.p);
    }
    for(unsigned int inotri=0;inotri<3;++inotri){
      const unsigned int ip0 = aTri[it0*3+inotri];
      ::glTexCoord2dv(aTexP.data()+ip0*2);
      ::glVertex3dv(aXYZ.data()+ip0*3);
    }
  }
  ::glEnd();
}

void Draw2(
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& aTri,
    const std::vector<double>& aTexP,
    const std::vector<unsigned int>& elsup_ind,
    const std::vector<unsigned int>& elsup,
    const dfm2::CVec3d lc[4])
{
  ::glEnable(GL_TEXTURE_2D);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,1,1);
  ::glBegin(GL_TRIANGLES);
  for(unsigned int it0=0;it0<aTri.size()/3;++it0){
    { // check the distortion
      const unsigned int i0 = aTri[it0*3+0];
      const unsigned int i1 = aTri[it0*3+1];
      const unsigned int i2 = aTri[it0*3+2];
      const double area2 = dfm2::Area_Tri2(aTexP.data()+i0*2,  aTexP.data()+i1*2, aTexP.data()+i2*2);
      const double area3 = dfm2::Area_Tri3(aXYZ.data()+i0*3,  aXYZ.data()+i1*3, aXYZ.data()+i2*3);
      const double scoreArea = 0.5*(area2/area3 + area3/area2);
      if( scoreArea < 0 || scoreArea > 2.0 ){ continue; }
      const double len12 = dfm2::Distance2(aTexP.data()+i1*2, aTexP.data()+i2*2);
      const double len20 = dfm2::Distance2(aTexP.data()+i2*2, aTexP.data()+i0*2);
      const double len01 = dfm2::Distance2(aTexP.data()+i0*2, aTexP.data()+i1*2);
      const double Len12 = dfm2::Distance3(aXYZ.data()+i1*3, aXYZ.data()+i2*3);
      const double Len20 = dfm2::Distance3(aXYZ.data()+i2*3, aXYZ.data()+i0*3);
      const double Len01 = dfm2::Distance3(aXYZ.data()+i0*3, aXYZ.data()+i1*3);
      if( 0.5*(len12/Len12+Len12/len12) > 1.5 ){ continue; }
      if( 0.5*(len20/Len20+Len20/len20) > 1.5 ){ continue; }
      if( 0.5*(len01/Len01+Len01/len01) > 1.5 ){ continue; }
    }
    for(unsigned int inotri=0;inotri<3;++inotri){
      const unsigned int ip0 = aTri[it0*3+inotri];
      const dfm2::CVec3d p0(aXYZ.data()+ip0*3);
      double tx = aTexP[ip0*2+0];
      double ty = aTexP[ip0*2+1];
      ::glTexCoord2d(tx,ty);
      dfm2::CVec3d p = lc[3] + lc[2]*0.05 + lc[0]*tx + lc[1]*ty;
      ::glVertex3dv(p.p);
    }
  }
  ::glEnd();
}

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;

  delfem2::Read_Ply(
      std::string(PATH_INPUT_DIR)+"/bunny_34k.ply",
//      std::string(PATH_INPUT_DIR)+"/arm_16k.ply",
      aXYZ,aTri);
  delfem2::Normalize_Points3(aXYZ);
  dfm2::Rotate_Points3(
      aXYZ,
      M_PI*1.0,M_PI*1.5,0.);

  std::vector<unsigned int> elsup_ind, elsup;
  dfm2::JArray_ElSuP_MeshElem(elsup_ind, elsup,
      aTri.data(), aTri.size()/3, 3, aXYZ.size()/3);

  std::vector<unsigned int> aTriSuTri;
  ElSuEl_MeshElem(aTriSuTri,
      aTri.data(), aTri.size()/3,
      delfem2::MESHELEM_TRI,
      aXYZ.size()/3);

  unsigned int ielm_ker = 10010;
  dfm2::CVec3d coordLocal[4]; // ex, ey, ez, org
  std::vector<double> aTexP;
  {
    std::vector<double> aTexE; // element-wise texture coordinate
    dfm2::CExpMap_DijkstraElem expmap(
        aTexE,
        ielm_ker, aXYZ, aTri, aTriSuTri);
    std::vector<double> aDist;
    std::vector<unsigned int> aOrder;
    dfm2::DijkstraElem_MeshElemGeo3(
        aDist, aOrder, expmap,
        ielm_ker,
        aTri, aTri.size() / 3,
        aXYZ,
        aTriSuTri);
    coordLocal[0] = expmap.aAxisX[ielm_ker];
    coordLocal[2] = dfm2::Normal_Tri3(ielm_ker,aTri,aXYZ).Normalize();
    assert( fabs(coordLocal[0] * coordLocal[2])<1.0e-10 );
    assert( fabs(coordLocal[0].Length()-1.0)<1.0e-10 );
    assert( fabs(coordLocal[2].Length()-1.0)<1.0e-10 );
    coordLocal[1] = coordLocal[2]^coordLocal[0];
    coordLocal[3] = dfm2::CG_Tri3(ielm_ker,aTri,aXYZ);
    //
    aTexP.resize(aXYZ.size()/3*2);
    for(unsigned int ip0=0;ip0<aXYZ.size()/3;++ip0) {
      const dfm2::CVec3d p0(aXYZ.data()+ip0*3);
      double tex[2] = {0, 0};
      double w0 = 0.0;
      for (unsigned int ielsup = elsup_ind[ip0]; ielsup < elsup_ind[ip0 + 1]; ++ielsup) {
        const unsigned int it1 = elsup[ielsup];
        dfm2::CVec3d p1 = dfm2::CG_Tri3(it1, aTri, aXYZ);
        double w1 = 1.0 / (p0 - p1).Length();
        tex[0] += w1 * aTexE[it1 * 2 + 0];
        tex[1] += w1 * aTexE[it1 * 2 + 1];
        w0 += w1;
      }
      tex[0] /= w0;
      tex[1] /= w0;
      aTexP[ip0*2+0] = tex[0];
      aTexP[ip0*2+1] = tex[1];
    }
  }

  // ------

//  std::random_device rd;
//  std::mt19937 rdeng(rd());
//  std::uniform_int_distribution<unsigned int> ncluster_gen(1,100);

  // above: data preparation
  // -----------------------
  // below: view

  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  delfem2::opengl::setSomeLighting();
  viewer.camera.view_height = 0.5;
  viewer.camera.camera_rot_mode  = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;

  int m_texName;
  {
    unsigned int w, h;
    std::vector<unsigned char> image;
    dfm2::LoadImage_PPMAscii(w, h, image,
        std::string(PATH_INPUT_DIR) + "/dep.ppm");
    assert(image.size() == w * h * 3);
    m_texName = dfm2::opengl::SetTexture_RGB(w,h,image);
    {
      glEnable(GL_TEXTURE_2D);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
      glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
      glShadeModel(GL_SMOOTH);
      glBindTexture(GL_TEXTURE_2D , m_texName);
    }
  }
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    //
    for(unsigned int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      ::glDisable(GL_LIGHTING);
      ::glColor3d(0, 0, 0);
      delfem2::opengl::DrawMeshTri3D_Edge(aXYZ, aTri);
      {
        ::glDisable(GL_TEXTURE_2D);
        ::glDisable(GL_LIGHTING);
        ::glColor3d(1,0,0);
        ::glBegin(GL_TRIANGLES);
        unsigned int i0 = aTri[ielm_ker*3+0];
        unsigned int i1 = aTri[ielm_ker*3+1];
        unsigned int i2 = aTri[ielm_ker*3+2];
        ::glVertex3dv(aXYZ.data()+i0*3);
        ::glVertex3dv(aXYZ.data()+i1*3);
        ::glVertex3dv(aXYZ.data()+i2*3);
        ::glEnd();
      }
      Draw(aXYZ,aTri,aTexP,elsup_ind,elsup);
      Draw2(aXYZ,aTri,aTexP,elsup_ind,elsup,coordLocal);
//      delfem2::opengl::DrawMeshTri3DFlag_FaceNorm(aXYZ, aTri, aFlgElm, aColor);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
    }
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
