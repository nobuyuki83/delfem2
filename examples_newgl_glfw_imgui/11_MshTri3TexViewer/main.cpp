

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
#  include <glad/glad.h>
#endif
#include <GLFW/glfw3.h>

#include "ImGuiFileDialog/ImGuiFileDialog.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/msh_points.h"
#include "delfem2/opengl/new/drawer_mshtex.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#define  STB_IMAGE_IMPLEMENTATION
#include "delfem2/openglstb/img2tex.h"

namespace dfm2 = delfem2;
// ---------------------------------------------------------------

delfem2::glfw::CViewer3 viewer;
dfm2::opengl::Drawer_MeshTex drawer_tritex;
unsigned int id_tex;

void OpenFile(const std::string& filePathName){
  std::cout << filePathName << std::endl;
  std::string fname_mtl;
  {
    std::vector<double> vtx_xyz, vtx_tex, vtx_nrm;
    std::vector<unsigned int> elem_vtx_index;
    std::vector<unsigned int> elem_vtx_xyz, elem_vtx_tex, elem_vtx_nrm;
    std::vector<std::string> group_names;
    std::vector<unsigned int> group_elem_index;
    delfem2::Read_WavefrontObjWithMaterialMixedElem(
        fname_mtl,
        vtx_xyz, vtx_tex, vtx_nrm,
        elem_vtx_index,
        elem_vtx_xyz, elem_vtx_tex, elem_vtx_nrm,
        group_names, group_elem_index, filePathName);
    delfem2::Normalize_Points3(
        vtx_xyz,
        1.);
    std::vector<double> uni_xyz0,uni_tex0;
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
    if( !::glIsTexture(id_tex) ) { ::glGenTextures(1, &id_tex); }
    ::glBindTexture(GL_TEXTURE_2D, id_tex);
    delfem2::openglstb::LoadImageFileSetToTexture(
        texPath.string().c_str());
  }
  // -------------------
}

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
    if (ImGui::Button("Load Texture")) {
      ImGuiFileDialog::Instance()->OpenDialog(
          "ChooseTexFileDlgKey",
          "Choose File",
          "{.png,.jpg}",
          ".");
    }
    if (ImGuiFileDialog::Instance()->Display("ChooseMeshFileDlgKey")) {
      if (ImGuiFileDialog::Instance()->IsOk()) {
        std::string filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
        OpenFile(filePathName);
      }
      ImGuiFileDialog::Instance()->Close();
    }
    if (ImGuiFileDialog::Instance()->Display("ChooseTexFileDlgKey")) {
      if (ImGuiFileDialog::Instance()->IsOk()) {
        std::string filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
        std::cout << "texPath: " << filePathName << std::endl;
        if( std::filesystem::exists(filePathName) ) {
          if (!::glIsTexture(id_tex)) { ::glGenTextures(1, &id_tex); }
          ::glBindTexture(GL_TEXTURE_2D, id_tex);
          delfem2::openglstb::LoadImageFileSetToTexture(
              filePathName.c_str());
        }
      }
      ImGuiFileDialog::Instance()->Close();
    }
    ImGui::End();
  }

  {
    ::glEnable(GL_DEPTH_TEST);
    ::glEnable(GL_TEXTURE_2D);
    if( glIsTexture(id_tex) ) {
      ::glad_glBindTexture(GL_TEXTURE_2D, id_tex);
    }
    drawer_tritex.Draw(
        viewer.GetProjectionMatrix().data(),
        viewer.GetModelViewMatrix().data());
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
  drawer_tritex.InitGL();

  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGui_ImplGlfw_InitForOpenGL(viewer.window, true);// Setup Platform/Renderer bindings
  ImGui::StyleColorsClassic(); // Setup Dear ImGui style

  OpenFile(
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
