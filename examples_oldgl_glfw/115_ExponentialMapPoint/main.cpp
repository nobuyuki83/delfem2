/*
 * Copyright (c) 2019-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <string>
#include <cassert>
#include <cstdlib>
#include <set>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/expmap_geo3dijk.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/mshmisc.h"
#include "delfem2/color.h"
#include "delfem2/img_ioppm.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/tex.h"

namespace dfm2 = delfem2;

void myGlutDisplay(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aNorm,
    const std::vector<double> &aTex) {
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < aTri.size() / 3; itri++) {
//  for(unsigned int iit=0;iit<aIndTri.size();++iit){
//    unsigned int itri = aIndTri[iit];
    const unsigned int i0 = aTri[itri * 3 + 0];
    const unsigned int i1 = aTri[itri * 3 + 1];
    const unsigned int i2 = aTri[itri * 3 + 2];
    double r = 2.0;
    //
    ::glNormal3dv(aNorm.data() + i0 * 3);
    ::glTexCoord2d(aTex[i0 * 2 + 0] * r, aTex[i0 * 2 + 1] * r);
    ::glVertex3dv(aXYZ.data() + i0 * 3);
    //
    ::glNormal3dv(aNorm.data() + i1 * 3);
    ::glTexCoord2d(aTex[i1 * 2 + 0] * r, aTex[i1 * 2 + 1] * r);
    ::glVertex3dv(aXYZ.data() + i1 * 3);
    //
    ::glNormal3dv(aNorm.data() + i2 * 3);
    ::glTexCoord2d(aTex[i2 * 2 + 0] * r, aTex[i2 * 2 + 1] * r);
    ::glVertex3dv(aXYZ.data() + i2 * 3);
  }
  ::glEnd();
}



// ---------------------------

int main() {
  std::vector<double> vtx_xyz;
  std::vector<unsigned int> tri_vtx;

  delfem2::Read_Ply(
//      std::string(PATH_INPUT_DIR)+"/bunny_34k.ply",
      vtx_xyz, tri_vtx,
      std::filesystem::path(PATH_INPUT_DIR) / "arm_16k.ply");
  delfem2::Normalize_Points3(vtx_xyz);

  std::vector<unsigned int> psup_ind, psup;
  {
    std::vector<unsigned int> elsup_ind, elsup;
    dfm2::JArray_ElSuP_MeshElem(
        elsup_ind, elsup,
        tri_vtx.data(), tri_vtx.size() / 3, 3,
        vtx_xyz.size() / 3);
    dfm2::JArrayPointSurPoint_MeshOneRingNeighborhood(
        psup_ind, psup,
        tri_vtx.data(), elsup_ind, elsup, 3,
        vtx_xyz.size() / 3);
  }
  std::vector<double> vtx_normal(vtx_xyz.size());
  delfem2::Normal_MeshTri3D(
      vtx_normal.data(),
      vtx_xyz.data(), vtx_xyz.size() / 3,
      tri_vtx.data(), tri_vtx.size() / 3 );

  unsigned int ip_ker = 0;
  std::vector<double> aTex;
  {
    dfm2::CExpMap_DijkstraPoint expmap(
        ip_ker,
        vtx_xyz, vtx_normal, tri_vtx, aTex,
        psup_ind, psup);
    std::vector<unsigned int> mapIp2Io;
    std::vector<double> aDist;
    dfm2::DijkstraPoint_MeshTri3D(
        aDist, mapIp2Io, expmap,
        ip_ker, vtx_xyz, psup_ind, psup);
  }

  // above: data preparation
  // -----------------------
  // below: view

  delfem2::glfw::CViewer3 viewer(0.5);
  //
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();

  int m_texName;
  {
    unsigned int image_width, image_height;
    std::vector<unsigned char> image_data_rgb;
    dfm2::LoadImage_PPMAscii(
        image_width, image_height, image_data_rgb,
        std::string(PATH_INPUT_DIR) + "/dep.ppm");
    assert(image_data_rgb.size() == image_width * image_height * 3);
    m_texName = dfm2::opengl::SetTexture_RGB(image_width, image_height, image_data_rgb);
    {
      glEnable(GL_TEXTURE_2D);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
      glShadeModel(GL_SMOOTH);
      glBindTexture(GL_TEXTURE_2D, m_texName);
    }
  }

  while (!glfwWindowShouldClose(viewer.window)) {
    std::vector<std::pair<double, dfm2::CColor> > colorMap;
    dfm2::ColorMap_BlueCyanGreenYellowRed(colorMap, 0, 1);
    for (unsigned int iframe = 0; iframe < 30; ++iframe) {
      viewer.DrawBegin_oldGL();
      {
        ::glDisable(GL_LIGHTING);
        ::glDisable(GL_TEXTURE_2D);
        ::glColor3d(0, 0, 0);
        delfem2::opengl::DrawMeshTri3D_Edge(vtx_xyz, tri_vtx);
      }
      {
        ::glEnable(GL_TEXTURE_2D);
        ::glEnable(GL_LIGHTING);
        float gray[4] = {1, 1, 1, 1};
        ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, gray);
        float shine[4] = {1, 1, 1, 1};
        ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
        ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128.0);
      }
      myGlutDisplay(vtx_xyz, tri_vtx, vtx_normal, aTex);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
    }
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
