

#include <cstdio>
#include <iostream>
#include <vector>
#include <set>
#if defined(_MSC_VER)
#include <windows.h>
#endif
//
#ifdef EMSCRIPTEN
#include <emscripten/emscripten.h>
#define GLFW_INCLUDE_ES3
#define GL_GLEXT_PROTOTYPES
#define EGL_EGLEXT_PROTOTYPES
#else
#include <glad/glad.h>
#endif
#include <GLFW/glfw3.h>

#include "ImGuiFileDialog/ImGuiFileDialog.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mshuni.h"
#include "delfem2/srchuni_v3.h"
#include "delfem2/points.h"
#include "delfem2/opengl/new/shdr_mshtex.h"
#include "delfem2/opengl/new/drawer_sphere.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "delfem2/openglstb/img2tex.h"
#include "stb/stb_image_write.h"

namespace dfm2 = delfem2;
// ---------------------------------------------------------------

class MyViewer : public delfem2::glfw::CViewer3{
 public:
  void OpenFile(const std::string& filePathName){

    std::cout << filePathName << std::endl;
    std::string fname_mtl;
    std::vector<unsigned int> elem_vtx_index;
    std::vector<std::string> group_names;
    std::vector<unsigned int> group_elem_index;
    const std::filesystem::path file_path;
    delfem2::Read_WavefrontObjWithMaterialMixedElem(
        fname_mtl,
        vtx_xyz, vtx_tex, vtx_nrm,
        elem_vtx_index,
        elem_vtx_xyz, elem_vtx_tex, elem_vtx_nrm,
        group_names, group_elem_index, filePathName);
    delfem2::Normalize_Points3(
        vtx_xyz,
        1.);
    drawer_tritex.SetMesh(vtx_xyz, vtx_tex, elem_vtx_xyz, elem_vtx_tex);
    // ------------
    auto mtlFilePathName = std::filesystem::path(filePathName).parent_path() / fname_mtl;
    std::cout << "material name: " << mtlFilePathName << std::endl;
    std::vector<delfem2::MaterialWavefrontObj> materials;
    delfem2::Read_WavefrontMaterial(
        mtlFilePathName,
        materials);
    std::cout << materials.size() << std::endl;
    for(const auto& mat : materials ) {
      auto texPath = std::filesystem::path(filePathName).parent_path() / mat.map_Kd;
      std::cout << "texPath: " << texPath << std::endl;
      if( !std::filesystem::exists(texPath) ){ continue; }
      if( !::glIsTexture(id_tex_org) ) { ::glGenTextures(1, &id_tex_org); }
      ::glBindTexture(GL_TEXTURE_2D, id_tex_org);
      delfem2::openglstb::LoadImageFileSetToTexture(
          texPath.string().c_str());
      break;
    }
    if( !::glIsTexture(id_tex_ano) ) { ::glGenTextures(1, &id_tex_ano); }
    ::glBindTexture(GL_TEXTURE_2D, id_tex_ano);
  }

  void mouse_press(
      const float src[3], const float dir[3]) override {
    std::map<double, delfem2::PointOnSurfaceMesh<double>> mapDepthPES;
    delfem2::IntersectionRay_MeshTri3(
        mapDepthPES,
        dfm2::CVec3d(src), dfm2::CVec3d(dir), elem_vtx_xyz, vtx_xyz,
        1.e-3);
    if( mapDepthPES.empty() ){ return; }
    delfem2::PointOnSurfaceMesh<double> pom = mapDepthPES.begin()->second;
    const unsigned int ip0 = elem_vtx_tex[pom.itri*3+0];
    const unsigned int ip1 = elem_vtx_tex[pom.itri*3+1];
    const unsigned int ip2 = elem_vtx_tex[pom.itri*3+2];
    pick_uv[0] = pom.r0*vtx_tex[ip0*2+0] + pom.r1*vtx_tex[ip1*2+0] + (1-pom.r0-pom.r1)*vtx_tex[ip2*2+0];
    pick_uv[1] = pom.r0*vtx_tex[ip0*2+1] + pom.r1*vtx_tex[ip1*2+1] + (1-pom.r0-pom.r1)*vtx_tex[ip2*2+1];
    std::cout << pick_uv[0] << " " << pick_uv[1] << std::endl;
    //
    std::cout << FramebufferName << std::endl;
    std::cout << std::get<0>(texInfo) << " " << std::get<1>(texInfo) << std::endl;
    glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);
    ::glActiveTexture(GL_TEXTURE0);
    if( glIsTexture(id_tex_ano) ) {
      ::glBindTexture(GL_TEXTURE_2D, id_tex_ano);
    }
    glViewport(0,0, std::get<0>(texInfo),std::get<1>(texInfo));
    drawer_sphere.color = {1,0,0,1};
    drawer_sphere.Draw(
        0.1,
        std::array<float,3>{
          2*pick_uv[0]-1,
          2*pick_uv[1]-1,0.0}.data(),
        dfm2::CMat4f::Identity().data(),
        dfm2::CMat4f::Identity().data() );
    ::glfwSwapBuffers(window);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    {
      int width, height;
      ::glfwGetFramebufferSize(window, &width, &height);
      glViewport(0,0,width,height);
    }
  }

  void CompileShader() {
    drawer_sphere.InitGL();
    drawer_tritex.CompileShader();
    if (!::glIsTexture(id_tex_ano)) { ::glGenTextures(1, &id_tex_ano); }
    ::glGenFramebuffers(1, &FramebufferName);
    ::glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);
    ::glBindTexture(GL_TEXTURE_2D, id_tex_ano);
    ::glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, id_tex_ano, 0);
    ::glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

 public:
  std::vector<double> vtx_xyz, vtx_tex, vtx_nrm;
  std::vector<unsigned int> elem_vtx_xyz, elem_vtx_tex, elem_vtx_nrm;
  std::array<float,2> pick_uv;
  //
  dfm2::opengl::Drawer_MeshTexSeparateIndexingMixTwoTexture drawer_tritex;
  float slider_mix_ratio;
  unsigned int id_tex_org = 0;  // original texture
  unsigned int id_tex_ano = 0;  // annotation
  //
  std::tuple<int,int,int> texInfo;
  unsigned int FramebufferName;
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
    // open Dialog Simple
    if (ImGui::Button("Load Mesh")) {
      ImGuiFileDialog::Instance()->OpenDialog(
          "ChooseMeshFileDlgKey",
          "Choose File",
          ".obj",
          ".");
    }
    if (ImGui::Button("Load Annotation")) {
      ImGuiFileDialog::Instance()->OpenDialog(
          "ChooseTexFileDlgKey",
          "Choose File",
          "{.png,.jpg}",
          ".");
    }
    if (ImGui::Button("Save Annotation")) {
      ImGuiFileDialog::Instance()->OpenDialog(
          "SaveTexFileDlgKey",
          "Choose File",
          "{.png,.jpg}",
          ".");
    }
    ImGui::SliderFloat("slider 1", &viewer.slider_mix_ratio, 0.0f, 1.0f);
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
        if( std::filesystem::exists(filePathName) ) {
          assert( ::glIsTexture(viewer.id_tex_ano) );
          ::glBindTexture(GL_TEXTURE_2D, viewer.id_tex_ano);
          viewer.texInfo = delfem2::openglstb::LoadImageFileSetToTexture(
              filePathName.c_str());
        }
      }
      ImGuiFileDialog::Instance()->Close();
    }
    if (ImGuiFileDialog::Instance()->Display("SaveTexFileDlgKey")) {
      if (ImGuiFileDialog::Instance()->IsOk()) {
        std::string filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
        std::cout << "texPath: " << filePathName << std::endl;
        if( std::get<2>(viewer.texInfo) == 3 ) {
          std::vector<char> pixels(
              std::get<0>(viewer.texInfo) * std::get<1>(viewer.texInfo) * 3, 0);
          assert(::glIsTexture(viewer.id_tex_ano));
          ::glBindTexture(GL_TEXTURE_2D, viewer.id_tex_ano);
          glGetTexImage(
              GL_TEXTURE_2D,
              0,
              GL_RGB,
              GL_UNSIGNED_BYTE,
              pixels.data());
          stbi_flip_vertically_on_write(1);
          stbi_write_png(
              filePathName.c_str(),
              std::get<0>(viewer.texInfo),
              std::get<1>(viewer.texInfo),
              3,
              pixels.data(),
              std::get<0>(viewer.texInfo) * 3 );
        }
      }
      ImGuiFileDialog::Instance()->Close();
    }
    ImGui::End();
  }

  {
    ::glEnable(GL_DEPTH_TEST);
    ::glEnable(GL_TEXTURE_2D);
    ::glActiveTexture(GL_TEXTURE0);
    if( glIsTexture(viewer.id_tex_org) ) {
      ::glBindTexture(GL_TEXTURE_2D, viewer.id_tex_org);
    }
    ::glActiveTexture(GL_TEXTURE1);
    if( glIsTexture(viewer.id_tex_ano) ) {
      ::glBindTexture(GL_TEXTURE_2D, viewer.id_tex_ano);
    }
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
  viewer.InitGL();
#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif
  viewer.CompileShader();

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
