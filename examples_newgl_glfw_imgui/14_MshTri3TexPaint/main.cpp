

#include <cstdio>
#include <iostream>
#include <vector>
#include <set>
#if defined(_MSC_VER)
#  include <windows.h>
#endif
//
#ifdef EMSCRIPTEN
#  include <emscripten/emscripten.h>
#  define GLFW_INCLUDE_ES3
#  define GL_GLEXT_PROTOTYPES
#  define EGL_EGLEXT_PROTOTYPES
#else
#include <glad/glad.h>
#endif
#include <GLFW/glfw3.h>

//#include "ImGuiFileDialog/ImGuiFileDialog.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/srch_bruteforce.h"
#include "delfem2/opengl/new/drawer_mshtex.h"
#include "delfem2/opengl/new/drawer_sphere.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "delfem2/openglstb/img2tex.h"
#include "stb/stb_image_write.h"

namespace dfm2 = delfem2;
// ---------------------------------------------------------------

class MyViewer : public delfem2::glfw::CViewer3 {
 public:
  void OpenFile(const std::string &obj_file_path) {
    this->obj_file_path_ = obj_file_path;
    std::string fname_mtl;
    {
      std::vector<unsigned int> elem_vtx_index;
      std::vector<std::string> group_names;
      std::vector<unsigned int> group_elem_index;
      delfem2::Read_WavefrontObjWithMaterialMixedElem(
          fname_mtl,
          vtx_xyz, vtx_tex, vtx_nrm,
          elem_vtx_index,
          elem_vtx_xyz, elem_vtx_tex, elem_vtx_nrm,
          group_names, group_elem_index, obj_file_path);
      delfem2::Normalize_Points3(
          vtx_xyz,
          1.);
      std::vector<double> uni_xyz0, uni_tex0;
      std::vector<unsigned int> tri_uni, uni_xyz, uni_tex;
      delfem2::UnifySeparateIndexing_PosTex(
          uni_xyz0, uni_tex0,
          tri_uni, uni_xyz, uni_tex,
          vtx_xyz, vtx_tex,
          elem_vtx_xyz, elem_vtx_tex);
      drawer_tritex.SetElement(tri_uni, GL_TRIANGLES);
      drawer_tritex.SetCoords(uni_xyz0, 3);
      drawer_tritex.SetTexUV(uni_tex0);
    }
    // ------------
    auto mtlFilePathName = std::filesystem::path(obj_file_path).parent_path() / fname_mtl;
    std::vector<delfem2::MaterialWavefrontObj> materials;
    delfem2::Read_WavefrontMaterial(
        mtlFilePathName,
        materials);
    std::cout << materials.size() << std::endl;
    for (const auto &mat: materials) {
      tex_file_path_ = std::filesystem::path(obj_file_path).parent_path() / mat.map_Kd;
      if (!std::filesystem::exists(tex_file_path_)) {
        tex_file_path_ = "";
        continue;
      }
      ::glBindTexture(GL_TEXTURE_2D, id_tex_obj);
      delfem2::openglstb::LoadImageFileSetToTexture(
          tex_file_path_.c_str());
      break;
    }
    assert(::glIsTexture(id_tex_ano));
    ::glBindTexture(GL_TEXTURE_2D, id_tex_ano);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    {
      int N = 1024;
      texInfo = {N,N,3};
      std::vector<unsigned char> pixels(N*N*3, 255);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
                   N, N,
                   0, GL_RGB, GL_UNSIGNED_BYTE,
                   pixels.data());
    }
  }

  void mouse_press(
      const float src[3], const float dir[3]) override {
    std::map<double, delfem2::PointOnSurfaceMesh<double>> mapDepthPES;
    delfem2::IntersectionRay_MeshTri3(
        mapDepthPES,
        dfm2::CVec3d(src), dfm2::CVec3d(dir), elem_vtx_xyz, vtx_xyz,
        1.e-3);
    if (mapDepthPES.empty()) { return; }
    delfem2::PointOnSurfaceMesh<double> pom = mapDepthPES.begin()->second;
    const unsigned int ip0 = elem_vtx_tex[pom.itri * 3 + 0];
    const unsigned int ip1 = elem_vtx_tex[pom.itri * 3 + 1];
    const unsigned int ip2 = elem_vtx_tex[pom.itri * 3 + 2];
    pick_uv[0] =
        pom.r0 * vtx_tex[ip0 * 2 + 0] +
        pom.r1 * vtx_tex[ip1 * 2 + 0] +
        (1 - pom.r0 - pom.r1) * vtx_tex[ip2 * 2 + 0];
    pick_uv[1] =
        pom.r0 * vtx_tex[ip0 * 2 + 1] +
        pom.r1 * vtx_tex[ip1 * 2 + 1] +
        (1 - pom.r0 - pom.r1) * vtx_tex[ip2 * 2 + 1];
    std::cout << pick_uv[0] << " " << pick_uv[1] << std::endl;
    //
    std::cout << frame_buffer_idx << std::endl;
    std::cout << std::get<0>(texInfo) << " " << std::get<1>(texInfo) << std::endl;
    glBindFramebuffer(GL_FRAMEBUFFER, frame_buffer_idx);
    ::glActiveTexture(GL_TEXTURE0);
    ::glIsTexture(id_tex_ano);
    ::glBindTexture(GL_TEXTURE_2D, id_tex_ano);
    glViewport(0, 0, std::get<0>(texInfo), std::get<1>(texInfo));
    drawer_sphere.color = {1, 0, 0, 1};
    drawer_sphere.Draw(
        0.1,
        std::array<float, 3>{
            2 * pick_uv[0] - 1,
            2 * pick_uv[1] - 1, 0.0}.data(),
        dfm2::CMat4f::Identity().data(),
        dfm2::CMat4f::Identity().data());
    ::glfwSwapBuffers(window);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    {
      int width, height;
      ::glfwGetFramebufferSize(window, &width, &height);
      glViewport(0, 0, width, height);
    }
  }

  void InitGL() {
    drawer_sphere.InitGL();
    drawer_tritex.InitGL();
    if (!::glIsTexture(id_tex_ano)) { ::glGenTextures(1, &id_tex_ano); }
    if (!::glIsTexture(id_tex_obj)) { ::glGenTextures(1, &id_tex_obj); }
    ::glBindTexture(GL_TEXTURE_2D, id_tex_obj);
    ::glBindTexture(GL_TEXTURE_2D, id_tex_ano);
    assert( glIsTexture(id_tex_obj) );
    assert( glIsTexture(id_tex_ano) );
    ::glGenFramebuffers(1, &frame_buffer_idx);
    ::glBindFramebuffer(GL_FRAMEBUFFER, frame_buffer_idx);
    ::glBindTexture(GL_TEXTURE_2D, id_tex_ano);
    ::glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, id_tex_ano, 0);
    ::glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

  void SaveAnnotation() {
    assert(std::get<2>(texInfo) == 3);
    const int w0 = std::get<0>(texInfo);
    const int h0 = std::get<1>(texInfo);
    std::vector<char> pixels(w0 * h0 * 3, 0);
    assert(::glIsTexture(id_tex_ano));
    ::glBindTexture(GL_TEXTURE_2D, id_tex_ano);
    glGetTexImage(
        GL_TEXTURE_2D, 0,
        GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
    stbi_flip_vertically_on_write(1);
    stbi_write_png(
        ano_file_path_.c_str(),
        w0, h0, 3, pixels.data(), w0 * 3);
  }

 public:
  std::string obj_file_path_;
  std::string tex_file_path_;
  std::vector<double> vtx_xyz, vtx_tex, vtx_nrm;
  std::vector<unsigned int> elem_vtx_xyz, elem_vtx_tex, elem_vtx_nrm;
  std::array<float, 2> pick_uv = {0.f, 0.f};
  //
  dfm2::opengl::Drawer_MeshMixTwoTextures drawer_tritex;
  float slider_mix_ratio = 0.5f;
  unsigned int id_tex_obj = 0;  // original texture
  unsigned int id_tex_ano = 0;  // annotation
  //
  std::string ano_file_path_;
  std::tuple<int, int, int> texInfo = {0, 0, 0};
  unsigned int frame_buffer_idx = 0;
  dfm2::opengl::Drawer_Sphere drawer_sphere;
} viewer;



void draw(GLFWwindow *window) {
  glfwPollEvents();
  ::glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);

  // feed inputs to dear imgui, start new frame
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  { // render your GUI
    ImGui::Begin("Triangle Mesh");
    /*
    ImGui::Text(
        "%s", std::filesystem::path(viewer.obj_file_path_).filename().c_str());
    // open Dialog Simple
    if (ImGui::Button("Load Wavefront Obj file")) {
      ImGuiFileDialog::Instance()->OpenDialog(
          "ChooseMeshFileDlgKey",
          "Choose File",
          ".obj",
          ".");
    }
    ImGui::Text(
        "%s", std::filesystem::path(viewer.tex_file_path_).filename().c_str());
    if (ImGui::Button("Load texture image")) {
      std::string path_dir = ".";
      if( std::filesystem::exists( std::filesystem::path(viewer.obj_file_path_).remove_filename() ) ){
        path_dir = std::filesystem::path(viewer.obj_file_path_).remove_filename().string();
      }
      ImGuiFileDialog::Instance()->OpenDialog(
          "ChooseTexFileDlgKey",
          "Choose File",
          "{.png,.jpg}",
          path_dir.c_str());
    }
    ImGui::Separator();
    //
    ImGui::SliderFloat("slider 1", &viewer.slider_mix_ratio, 0.0f, 1.0f);
    ImGui::Separator();
    if (ImGui::Button("Load Annotation")) {
      ImGuiFileDialog::Instance()->OpenDialog(
          "ChooseAnnoFileDlgKey",
          "Choose Annotation File",
          ".png",
          ".");
    }
    if (ImGui::Button("Save Annotation As")) {
      ImGuiFileDialog::Instance()->OpenDialog(
          "SaveAnnoFileDlgKey",
          "Choose Annotation File",
          "{.png}",
          ".");
    }
    if( std::filesystem::exists(viewer.ano_file_path_) ) {
      ImGui::Text("%s", std::filesystem::path(viewer.ano_file_path_).filename().c_str());
      if (ImGui::Button("Save Annotation")) {
        viewer.SaveAnnotation();
      }
    }
    //
    if (ImGuiFileDialog::Instance()->Display("ChooseMeshFileDlgKey")) {
      if (ImGuiFileDialog::Instance()->IsOk()) {
        std::string filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
        viewer.OpenFile(filePathName);
      }
      ImGuiFileDialog::Instance()->Close();
    }
    if (ImGuiFileDialog::Instance()->Display("ChooseTexFileDlgKey")) {
      if (ImGuiFileDialog::Instance()->IsOk()) {
        std::string filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
        std::cout << "texPath: " << filePathName << std::endl;
        if (std::filesystem::exists(filePathName)) {
          assert(::glIsTexture(viewer.id_tex_obj));
          ::glBindTexture(GL_TEXTURE_2D, viewer.id_tex_obj);
          delfem2::openglstb::LoadImageFileSetToTexture(filePathName.c_str());
        }
      }
      ImGuiFileDialog::Instance()->Close();
    }
    if (ImGuiFileDialog::Instance()->Display("ChooseAnnoFileDlgKey")) {
      if (ImGuiFileDialog::Instance()->IsOk()) {
        viewer.ano_file_path_ = ImGuiFileDialog::Instance()->GetFilePathName();
        if (std::filesystem::exists(viewer.ano_file_path_)) {
          assert(::glIsTexture(viewer.id_tex_ano));
          ::glBindTexture(GL_TEXTURE_2D, viewer.id_tex_ano);
          viewer.texInfo = delfem2::openglstb::LoadImageFileSetToTexture(
              viewer.ano_file_path_.c_str());
        }
      }
      ImGuiFileDialog::Instance()->Close();
    }
    if (ImGuiFileDialog::Instance()->Display("SaveAnnoFileDlgKey")) {
      if (ImGuiFileDialog::Instance()->IsOk()) {
        viewer.ano_file_path_ = ImGuiFileDialog::Instance()->GetFilePathName();
        viewer.SaveAnnotation();
      }
      ImGuiFileDialog::Instance()->Close();
    }
     */
    ImGui::End();
  }

  {
    ::glEnable(GL_DEPTH_TEST);
    ::glEnable(GL_TEXTURE_2D);
    assert(::glIsTexture(viewer.id_tex_obj));
    assert(::glIsTexture(viewer.id_tex_ano));
    ::glActiveTexture(GL_TEXTURE0);
    ::glBindTexture(GL_TEXTURE_2D, viewer.id_tex_obj);
    ::glActiveTexture(GL_TEXTURE1);
    ::glBindTexture(GL_TEXTURE_2D, viewer.id_tex_ano);
    viewer.drawer_tritex.Draw(
        viewer.GetProjectionMatrix().data(),
        viewer.GetModelViewMatrix().data(),
        viewer.slider_mix_ratio);
  }

  // Render dear imgui into screen
  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);

  glfwSwapBuffers(window);
}

int main(int, char **) {

  delfem2::glfw::InitGLNew();
  viewer.OpenWindow();
#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif
  viewer.InitGL();

  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGui_ImplGlfw_InitForOpenGL(viewer.window, true);// Setup Platform/Renderer bindings
  ImGui::StyleColorsClassic(); // Setup Dear ImGui style

  viewer.OpenFile(
      std::filesystem::path(PATH_SRC_DIR)
          / ".." / ".." / "test_inputs" / "Babi" / "Babi.obj");

#ifdef EMSCRIPTEN
  ImGui_ImplOpenGL3_Init("#version 100");
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, window, 60, 1);
#else
  ImGui_ImplOpenGL3_Init("#version 150");
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif

  // Cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  glfwDestroyWindow(viewer.window);
  glfwTerminate();

  return 0;
}
