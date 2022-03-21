/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#define STB_IMAGE_IMPLEMENTATION
#include "delfem2/cvcamera.h"
#include "delfem2/cam_projection.h"
#include "delfem2/cam_modelview.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/msh_normal.h"
#include "delfem2/glfw/util.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/tex.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/stb_opengl/img2tex.h"
#include "delfem2/pugixml/metashape_camera.h"

namespace delfem2 {

class Projection_FrameCamera : public delfem2::Projection {
 public:
  Projection_FrameCamera(
      int width0, int height0,
      float f, float cx, float cy,
      float depth) {
    const auto wf = static_cast<float>(width0);
    const auto hf = static_cast<float>(height0);
    const delfem2::CMat4f mT = delfem2::CMat4f::Translation({0.f, 0.f, depth});
    const delfem2::CMat4f glb2img =
        delfem2::Mat4_CameraInternal_MetashapePinhole(f, cx, cy, wf, hf);
    const delfem2::CMat4f img2scr =
        delfem2::Mat4_Image2Screen(wf, hf, -1);
    const delfem2::CMat4f zyinv = delfem2::CMat4f::ScaleXYZ(1, -1, -1);
    const delfem2::CMat4f zinv = delfem2::CMat4f::ScaleXYZ(1, 1, -1);
    mM = zinv * img2scr * glb2img * zyinv * mT;
  }

  [[nodiscard]] std::array<float, 16> GetMatrix([[maybe_unused]] float asp) const override {
    std::array<float, 16> m{};
    mM.CopyTo(m.data());
    return m;
  }
 public:
  delfem2::CMat4f mM;
};

}

void SetCameraToViewer(
    delfem2::glfw::CViewer3 &viewer,
    unsigned int idx_camera,
    const delfem2::MetashapeCamConfig &cam_config) {
  namespace dfm2 = delfem2;
  assert(idx_camera < cam_config.aCamera.size());
  const dfm2::MetashapeCamConfig::Camera &cam = cam_config.aCamera[idx_camera];
  assert(cam.id_sensor >= 0 && cam.id_sensor < static_cast<int>(cam_config.aSensor.size()));
  const dfm2::MetashapeCamConfig::Sensor &snsr = cam_config.aSensor[cam.id_sensor];
  static auto camera_width = viewer.width = snsr.width / 8;
  static auto camera_height = viewer.height = snsr.height / 8;

  ::glfwSetWindowSizeCallback(
      viewer.window, [](GLFWwindow *window, int width, int height) {
        camera_width = width;
        camera_height = height;
      });
  ::glfwSetWindowSize(
      viewer.window,
      static_cast<int>(camera_width),
      static_cast<int>(camera_height));
  ::glfwSetWindowTitle(
      viewer.window,
      (std::string("camera:") + std::to_string(idx_camera)).c_str());

  const auto[mat_modelview, mat_projection]
  = CameraTransformationMVP(cam_config, idx_camera);

  viewer.scale = 1.0;

  float view_height = cam_config.region_size.norm();
  float depth = view_height * snsr.f / static_cast<float>(snsr.height);

  {  // set modelview matrix
    auto *pmv = dynamic_cast<delfem2::ModelView_Trackball *>(viewer.view_rotation.get());
    const delfem2::CMat3f m3 = mat_modelview.GetMat3();
    const delfem2::CVec3f cnt = cam_config.region_center;
    const delfem2::CVec3f tmp1 = m3.MatVec(cnt.data());
    const auto q1 = m3.GetQuaternion();
    delfem2::Copy_Quat(pmv->quaternion, q1.data());
    pmv->anchor[0] = cnt.x;
    pmv->anchor[1] = cnt.y;
    pmv->anchor[2] = cnt.z;
    viewer.trans[0] = tmp1[0] + mat_modelview(0, 3);
    viewer.trans[1] = tmp1[1] + mat_modelview(1, 3);
    viewer.trans[2] = tmp1[2] + mat_modelview(2, 3) + depth;
  }

  { // set projection matrix
    const dfm2::MetashapeCamConfig::Camera &cam0 = cam_config.aCamera[idx_camera];
    const dfm2::MetashapeCamConfig::Sensor &snsr0 = cam_config.aSensor[cam0.id_sensor];
    const auto width0 = static_cast<float>(snsr0.width);
    const auto height0 = static_cast<float>(snsr0.height);
    viewer.projection = std::make_unique<delfem2::Projection_FrameCamera>(
        width0, height0,
        snsr0.f, snsr0.cx, snsr0.cy,
        -depth);
  }
}

class MyViewer :
    public delfem2::glfw::CViewer3 {
 public:
  explicit MyViewer(
      std::filesystem::path path_xml,
      std::string image_folder_name) :
      cam_config_path_xml(std::move(path_xml)),
      image_folder_name(std::move(image_folder_name)) {
    delfem2::ReadMetashapeCameraXml(cam_config, cam_config_path_xml.string().c_str());
    for (idx_camera = 0; idx_camera < cam_config.aCamera.size(); ++idx_camera) {
      if (cam_config.aCamera[idx_camera].id_sensor == -1) { continue; }
      break;
    }
    bgcolor[0] = 0.8;
  }

  void OpenWindow() override {
    CViewer3::OpenWindow();
  }

  void key_press(int key, [[maybe_unused]] int mods) override {
  }

  void mouse_press(
      [[maybe_unused]] const float src[3],
      [[maybe_unused]] const float dir[3]) override {
    namespace dfm2 = delfem2;
#ifdef IMGUI_VERSION
    if (ImGui::GetIO().WantCaptureMouse) { return; }
#endif
    if (idx_camera >= cam_config.aCamera.size()) { return; }
  }

  void mouse_drag(
      [[maybe_unused]] const float src0[3],
      [[maybe_unused]] const float src1[3],
      [[maybe_unused]] const float dir[3]) override {
    if (idx_camera >= cam_config.aCamera.size()) { return; }
  }

  void mouse_release() override {
#ifdef IMGUI_VERSION
    if (ImGui::GetIO().WantCaptureMouse) { return; }
#endif
    if (nav.imodifier == GLFW_MOD_ALT || nav.imodifier == GLFW_MOD_SHIFT) {
      const unsigned int icam_best = SnapCamera(
          this->GetModelViewMatrix(),
          this->cam_config);
      SetCamera(icam_best);
    }
    if (idx_camera >= cam_config.aCamera.size()) { return; }
  }

  void mouse_wheel(double) override {
    idx_camera = UINT_MAX;
    ::glfwSetWindowTitle(this->window, "free view");
  }

  void CursorPosition(double xpos, double ypos) override {
#ifdef IMGUI_VERSION
    if (ImGui::GetIO().WantCaptureMouse) { return; }
#endif
    ::delfem2::glfw::CViewer3::CursorPosition(xpos, ypos);
    if (this->nav.ibutton == GLFW_MOUSE_BUTTON_LEFT &&
        (nav.imodifier == GLFW_MOD_ALT || nav.imodifier == GLFW_MOD_SHIFT)) {
      idx_camera = UINT_MAX;
      ::glfwSetWindowTitle(this->window, "free view");
    }
  }

  void SetCamera(unsigned int idx_camera_) {
    this->idx_camera = idx_camera_ % cam_config.aCamera.size();
    SetCameraToViewer(
        *this,
        idx_camera, cam_config);
    // -----------------------
    if (idx_camera != idx_camera_texture) {  // set background image
      std::filesystem::path path_img = cam_config_path_xml.parent_path()
          / image_folder_name
          / (cam_config.aCamera[idx_camera].label + ".jpg");
      std::cout << path_img << std::endl;
      if (std::filesystem::exists(path_img)) {
        idx_camera_texture = idx_camera;
        if (!glIsTexture(id_texture)) { glGenTextures(1, &id_texture); }
        glBindTexture(GL_TEXTURE_2D, id_texture);
        auto imgshape = delfem2::openglstb::LoadImageFileSetToTexture(
            path_img.string().c_str());
        /*
        std::cout << cam_config.aCamera[idx_camera].label << " ";
        std::cout << std::get<0>(imgshape) << " ";
        std::cout << std::get<1>(imgshape) << " ";
        std::cout << std::get<2>(imgshape) << std::endl;
         */
      }
    }
  }

 public:
  // camera
  std::filesystem::path cam_config_path_xml;
  delfem2::MetashapeCamConfig cam_config{};  // camera configuration
  unsigned int idx_camera = 0;  // index of current camera

  // background image
  unsigned int id_texture = 0;  // id of texture for background image
  std::string image_folder_name;  // path for the background image
  unsigned int idx_camera_texture = UINT_MAX;  // index of camera which current texture has the image
};

void View(
    const std::filesystem::path &path_xml,
    const std::filesystem::path &path_obj,
    const std::filesystem::path &path_tex,
    const std::string &image_folder_name) {

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
    std::vector<unsigned int> elem_vtx_xyz;
    std::vector<unsigned int> elem_vtx_nrm;
    delfem2::Read_WavefrontObjWithMaterialMixedElem(
        fname_mtl,
        vtx_xyz, vtx_tex, vtx_nrm,
        elem_vtx_index,
        tri_vtx, tri_vtx_tex, elem_vtx_nrm,
        group_names, group_elem_index,
        path_obj);
    if (vtx_nrm.empty()) {
      vtx_nrm.resize(vtx_xyz.size());
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

  //head.Init(path_obj, path_tex);
  MyViewer viewer(
      path_xml,
      image_folder_name);
  // const float length_edge_ini = 0.05f;
  // -------------
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  viewer.SetCamera(viewer.idx_camera);
  tex0.InitGL();
  delfem2::opengl::setSomeLighting();
  while (!glfwWindowShouldClose(viewer.window)) {
    viewer.DrawBegin_oldGL();

    ::glEnable(GL_TEXTURE_2D);
    ::glBindTexture(GL_TEXTURE_2D, tex0.id_tex);
    delfem2::opengl::DrawMeshTri3D_FaceNorm_TexVtx(
        vtx_xyz, tri_vtx, vtx_tex, tri_vtx_tex);

    if (viewer.idx_camera_texture == viewer.idx_camera && viewer.idx_camera != UINT_MAX) {
      ::glBindTexture(GL_TEXTURE_2D, viewer.id_texture);
      delfem2::opengl::DrawTextureBackground();
    }
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
  std::cout << path_xml << std::endl;
  View(path_xml, path_obj, path_tex, "");
}


