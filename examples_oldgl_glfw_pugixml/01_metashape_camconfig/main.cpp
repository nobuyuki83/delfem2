/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include <sstream>
#include <vector>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#define STB_IMAGE_IMPLEMENTATION
#include "delfem2/msh_io_obj.h"
#include "delfem2/msh_normal.h"
#include "delfem2/metashape_camera.h"
#include "drawer_camera_config.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/stb_opengl/img2tex.h"
#include "delfem2/pugixml/metashape_camera.h"

void View(
    const std::filesystem::path &path_xml,
    const std::filesystem::path &path_obj,
    const std::filesystem::path &path_tex,
    float depth_scale) {
  delfem2::MetashapeCamConfig cam_config;
  delfem2::ReadMetashapeCameraXml(cam_config, path_xml.string().c_str());

  std::vector<double> vtx_xyz;
  std::vector<double> vtx_tex;
  std::vector<double> vtx_nrm;
  std::vector<unsigned int> tri_vtx;
  std::vector<unsigned int> tri_vtx_tex;
  {
    std::string fname_mtl;
    std::vector<std::string> group_names;
    std::vector<unsigned int> group_elem_index;
    std::vector<unsigned int> elem_vtx_index;
    std::vector<unsigned int> elem_vtx_nrm;
    delfem2::Read_WavefrontObjWithMaterialMixedElem(
      fname_mtl,
      vtx_xyz, vtx_tex, vtx_nrm,
      elem_vtx_index,
      tri_vtx, tri_vtx_tex, elem_vtx_nrm,
      group_names, group_elem_index,
      path_obj);
    if( vtx_nrm.empty() ) {
      vtx_nrm.resize( vtx_xyz.size() );
      delfem2::Normal_MeshTri3D(
          vtx_nrm.data(),
          vtx_xyz.data(), vtx_xyz.size() / 3,
          tri_vtx.data(), tri_vtx.size() / 3);
    }
  }
  delfem2::opengl::CTexRGB tex0;
  delfem2::openglstb::SetRgbToTex(
      tex0,
      path_tex.string(), true);
  //
  delfem2::glfw::CViewer3 viewer;
  bool is_view_texture = true;
  {  // viewer setting
    // setting anchor
    auto *pmv = dynamic_cast<delfem2::ModelView_Trackball *>(viewer.view_rotation.get());
    pmv->anchor[0] = cam_config.region_center.x;
    pmv->anchor[1] = cam_config.region_center.y;
    pmv->anchor[2] = cam_config.region_center.z;
    // setting keyboard callback
    std::function<void(int,int)> rr = [&is_view_texture](int key,int mods) {
      if( key == GLFW_KEY_T && mods == GLFW_MOD_SHIFT ){
        is_view_texture = !is_view_texture;
      }
    };
    viewer.keypress_callbacks.push_back(rr);
  }
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  tex0.InitGL();
  delfem2::opengl::setSomeLighting();

  while (!glfwWindowShouldClose(viewer.window)) {
    viewer.DrawBegin_oldGL();

    ::glDisable(GL_TEXTURE_2D);
    ::glDisable(GL_LIGHTING);
    for (const auto &cam: cam_config.aCamera) {
      DrawCamera(cam, cam_config, depth_scale);
    }
    ::glDisable(GL_LIGHTING);
    delfem2::opengl::DrawAxis(1);

    DrawRegion(cam_config);

    ::glDisable(GL_LIGHTING);
    ::glColor3d(1, 1, 1);

    if( is_view_texture ) {
      ::glEnable(GL_TEXTURE_2D);
      ::glBindTexture(GL_TEXTURE_2D, tex0.id_tex);
      delfem2::opengl::DrawMeshTri3D_FaceNorm_TexVtx(
          vtx_xyz, tri_vtx, vtx_tex, tri_vtx_tex);
    } else{
      ::glEnable(GL_LIGHTING);
      delfem2::opengl::DrawMeshTri3D_FaceNorm(
          vtx_xyz,
          tri_vtx,
          vtx_nrm);
    }

    // ------------
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
}

int main() {
  const auto dir_path = std::filesystem::path(PATH_SOURCE_DIR) / "mvs";
  const auto path_xml = dir_path / "camera.xml";
  const auto path_obj = dir_path / "foo.obj";
  const auto path_tex = dir_path / "foo.jpg";
  View(path_xml, path_obj, path_tex, 1);
}


