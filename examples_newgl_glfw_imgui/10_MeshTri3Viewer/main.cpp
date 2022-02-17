

#include <cstdio>
#include <iostream>
#include <vector>
#if defined(_MSC_VER)
#  include <windows.h>
#endif
//
#define GL_SILENCE_DEPRECATION
#ifdef EMSCRIPTEN
#  include <emscripten/emscripten.h>
#  define GLFW_INCLUDE_ES3
#  define GL_GLEXT_PROTOTYPES
#  define EGL_EGLEXT_PROTOTYPES
#else
#  include <glad/glad.h>
#endif
#include <GLFW/glfw3.h>

//#include "ImGuiFileDialog/ImGuiFileDialog.h"
//#include "ImGui-Addons/FileBrowser/ImGuiFileBrowser.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/opengl/new/drawer_mshtri.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

namespace dfm2 = delfem2;
// ---------------------------------------------------------------

delfem2::glfw::CViewer3 viewer;
dfm2::opengl::CShader_TriMesh drawer_mesh_edge;

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
    /*
    if (ImGui::Button("Open File Dialog")) {
      ImGuiFileDialog::Instance()->OpenDialog(
          "ChooseFileDlgKey",
          "Choose File",
          "{.ply,.obj}",
          ".");
    }
    if (ImGuiFileDialog::Instance()->Display("ChooseFileDlgKey")) {
      if (ImGuiFileDialog::Instance()->IsOk()) {
        std::string filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
        std::cout << filePathName << std::endl;
        std::vector<double> vtx_xyz;
        std::vector<unsigned int> tri_vtx;
        if( std::filesystem::path(filePathName).extension() == ".obj" ){
          delfem2::Read_Obj3(
              filePathName,
              vtx_xyz, tri_vtx);
        }
        if( std::filesystem::path(filePathName).extension() == ".ply" ) {
          delfem2::Read_Ply(
              vtx_xyz,
              tri_vtx,filePathName);
        }
        delfem2::Normalize_Points3(
            vtx_xyz,
            1.);
        drawer_mesh_edge.Initialize(vtx_xyz, 3, tri_vtx);
      }
      ImGuiFileDialog::Instance()->Close();
    }
     */
    ImGui::End();
  }

  {
    ::glEnable(GL_DEPTH_TEST);
    drawer_mesh_edge.Draw(
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
  drawer_mesh_edge.InitGL();

  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGui_ImplGlfw_InitForOpenGL(viewer.window, true);// Setup Platform/Renderer bindings
  ImGui::StyleColorsClassic(); // Setup Dear ImGui style

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
