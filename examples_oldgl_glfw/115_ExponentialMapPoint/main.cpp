/*
 * Copyright (c) 2019-2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/tex_gl.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/dijkstra.h"
#include "delfem2/points.h"
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/mshmisc.h"
#include "delfem2/color.h"
#include "delfem2/vec3.h"
#include "delfem2/imgio.h"
#include <GLFW/glfw3.h>
#include <vector>
#include <string>
#include <cassert>
#include <cstdlib>
#include <random>
#include <set>
#include <queue>
#include <climits>

namespace dfm2 = delfem2;

void myGlutDisplay(
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& aTri,
    const std::vector<double>& aNorm,
    const std::vector<double>& aTex)
{
  ::glBegin(GL_TRIANGLES);
  for(unsigned int itri=0;itri<aTri.size()/3;itri++){
//  for(unsigned int iit=0;iit<aIndTri.size();++iit){
//    unsigned int itri = aIndTri[iit];
    const unsigned int i0 = aTri[itri*3+0];
    const unsigned int i1 = aTri[itri*3+1];
    const unsigned int i2 = aTri[itri*3+2];
    double r = 2.0;
    //
    ::glNormal3dv(aNorm.data()+i0*3);
    ::glTexCoord2d(aTex[i0*2+0]*r,aTex[i0*2+1]*r);
    ::glVertex3dv(aXYZ.data()+i0*3);
    //
    ::glNormal3dv(aNorm.data()+i1*3);
    ::glTexCoord2d(aTex[i1*2+0]*r,aTex[i1*2+1]*r);
    ::glVertex3dv(aXYZ.data()+i1*3);
    //
    ::glNormal3dv(aNorm.data()+i2*3);
    ::glTexCoord2d(aTex[i2*2+0]*r,aTex[i2*2+1]*r);
    ::glVertex3dv(aXYZ.data()+i2*3);
  }
  ::glEnd();
}


class CExpMap{
public:
  const std::vector<double>& aXYZ;
  const std::vector<double>& aNorm;
  const std::vector<unsigned int>& aTri;
  std::vector<double>& aTex;
  const std::vector<unsigned int>& psup_ind;
  const std::vector<unsigned int>& psup;
  std::vector<double> aW;
  std::vector<dfm2::CVec3d> aAxisX;
public:
  CExpMap(
      unsigned int ip_ker,
      const std::vector<double>& aXYZ_,
      const std::vector<double>& aNorm_,
      const std::vector<unsigned int>& aTri_,
      std::vector<double>& aTex_,
      std::vector<unsigned int>& psup_ind_,
      std::vector<unsigned int>& psup_) :
      aXYZ(aXYZ_), aNorm(aNorm_), aTri(aTri_), aTex(aTex_),
      psup_ind(psup_ind_), psup(psup_)
  {
    const unsigned int np = aXYZ.size() / 3;
    aAxisX.resize(np, dfm2::CVec3d(0, 0, 0));
    aTex.resize(np * 2);
    aW.assign(np, 0.0);

    { // set kernel point information
      aTex[ip_ker * 2 + 0] = 0.0;
      aTex[ip_ker * 2 + 1] = 0.0;
      dfm2::CVec3d y0;
      dfm2::GetVertical2Vector(
          dfm2::CVec3d(aNorm.data() + ip_ker * 3),
          aAxisX[ip_ker], y0);
      aW[ip_ker] = 1.0;
    }
  }

  /**
   *
   * @param[in] ip0 point newly added
   * @param[in] io0 the order of newly added point. ( 0 if this is the first point)
   * @param[in] aOrder map from point index 2 the order of points added
   */
  void AddPoint(
      unsigned int ip0,
      unsigned int io0,
      std::vector<unsigned int>& aOrder)
  {
    assert( aOrder.size() == aXYZ.size()/3 );
    const dfm2::CVec3d n0 = dfm2::CVec3d(aNorm.data() + ip0 * 3).Normalize();
    aAxisX[ip0].SetNormalizedVector();
    const dfm2::CVec3d x0 = aAxisX[ip0];
    const dfm2::CVec3d y0 = n0 ^x0;
    aTex[ip0 * 2 + 0] /= aW[ip0];
    aTex[ip0 * 2 + 1] /= aW[ip0];
    for (unsigned int ipsup = psup_ind[ip0]; ipsup < psup_ind[ip0 + 1]; ++ipsup) {
      const unsigned int ip1 = psup[ipsup];
      if ( aOrder[ip1] < io0 ) { continue; } // effect propagate from fixed to unfixed
      const dfm2::CVec3d n1 = dfm2::CVec3d(aNorm.data() + ip1 * 3).Normalize();
      const dfm2::CVec3d d01 = dfm2::CVec3d(aXYZ.data() + ip1 * 3) - dfm2::CVec3d(aXYZ.data() + ip0 * 3);
      const double len01 = d01.Length();
      const dfm2::CVec3d e01 = (d01 - (d01 * n0) * n0).Normalize() * len01; // projected edge and same length
      const double w01 = 1.0 / len01;
      const dfm2::CVec3d x1 = dfm2::Mat3_MinimumRotation(n0, n1) * x0;
      aAxisX[ip1] += w01 * x1;
      aW[ip1] += w01;
      aTex[ip1 * 2 + 0] += w01 * (aTex[ip0 * 2 + 0] + (e01 * x0));
      aTex[ip1 * 2 + 1] += w01 * (aTex[ip0 * 2 + 1] + (e01 * y0));
    }
  }
};

// ---------------------------

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;

  delfem2::Read_Ply(
//      std::string(PATH_INPUT_DIR)+"/bunny_34k.ply",
      std::string(PATH_INPUT_DIR)+"/arm_16k.ply",
      aXYZ,aTri);
  delfem2::Normalize_Points3(aXYZ);

  std::vector<unsigned int> psup_ind, psup;
  {
    std::vector<unsigned int> elsup_ind, elsup;
    dfm2::JArray_ElSuP_MeshElem(elsup_ind, elsup,
        aTri.data(), aTri.size()/3, 3, aXYZ.size()/3);
    dfm2::JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
        aTri.data(), elsup_ind, elsup, 3, aXYZ.size()/3);
  }
  std::vector<double> aNorm(aXYZ.size());
  delfem2::Normal_MeshTri3D(
      aNorm.data(),
      aXYZ.data(), aXYZ.size() / 3,
      aTri.data(), aTri.size() / 3);

  unsigned int ip_ker = 0;
  std::vector<double> aTex;
  {
    CExpMap expmap(ip_ker,
                   aXYZ,aNorm,aTri,aTex,
                   psup_ind,psup);
    std::vector<unsigned int> mapIp2Io;
    std::vector<double> aDist;
    dfm2::DijkstraPoint_MeshTri3D(
        aDist, mapIp2Io, expmap,
        ip_ker, aXYZ, psup_ind, psup);
  }

  // above: data preparation
  // -----------------------
  // below: view

  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  delfem2::opengl::setSomeLighting();
  viewer.nav.camera.view_height = 0.5;
  viewer.nav.camera.camera_rot_mode  = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;

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
    std::vector<std::pair<double, dfm2::CColor> > colorMap;
    dfm2::ColorMap_BlueCyanGreenYellowRed(colorMap, 0,1);
    for(unsigned int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      {
        ::glDisable(GL_LIGHTING);
        ::glDisable(GL_TEXTURE_2D);
        ::glColor3d(0, 0, 0);
        delfem2::opengl::DrawMeshTri3D_Edge(aXYZ, aTri);
      }
      {
        ::glEnable(GL_TEXTURE_2D);
        ::glEnable(GL_LIGHTING);
        float gray[4] = {1,1,1,1};
        ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, gray);
        float shine[4] = {1,1,1,1};
        ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
        ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128.0);
      }
      myGlutDisplay( aXYZ,aTri,aNorm,aTex);
      {
      }
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
    }
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
