

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

//#include "ImGuiFileDialog/ImGuiFileDialog.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/opengl/new/drawer_mshtex.h"
#include "delfem2/opengl/new/drawer_sphere.h"
#include "delfem2/glfw/viewer2.h"
#include "delfem2/glfw/util.h"
#define  STB_IMAGE_IMPLEMENTATION
#include "delfem2/openglstb/img2tex.h"

namespace dfm2 = delfem2;
// ---------------------------------------------------------------

class MyViewer : public delfem2::glfw::CViewer2 {
 public:
  void mouse_press(const float src[2]) override {
    glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);
    glViewport(0,0, std::get<0>(texInfo),std::get<1>(texInfo));
    drawer_sphere.color = {1,0,0,1};
    drawer_sphere.Draw(
        0.1,
        std::array<float,3>{src[0],src[1],0.0}.data(),
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

  void InitGL(){
    drawer_sphere.InitGL();
    drawer_rect.InitGL();
    if (!::glIsTexture(id_tex)) { ::glGenTextures(1, &id_tex); }
    ::glGenFramebuffers(1, &FramebufferName);
    ::glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);
    ::glBindTexture(GL_TEXTURE_2D, id_tex);
    ::glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, id_tex, 0);
    ::glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

  void Draw()  {
    ::glEnable(GL_DEPTH_TEST);
    ::glEnable(GL_TEXTURE_2D);
    if( glIsTexture(id_tex) ) {
      ::glad_glBindTexture(GL_TEXTURE_2D, id_tex);
    }
    drawer_rect.Draw(
        GetProjectionMatrix().data(),
        GetModelViewMatrix().data());
  }

  void SetTextureImageFile(const std::filesystem::path& filePathName){
    if( std::filesystem::exists(filePathName) ) {
      if ( !::glIsTexture(id_tex)) { ::glGenTextures(1, &id_tex); }
      ::glBindTexture(GL_TEXTURE_2D, id_tex);
      texInfo = delfem2::openglstb::LoadImageFileSetToTexture(
          filePathName.c_str());
    }
  }

 public:
  unsigned int FramebufferName = 0;
  dfm2::opengl::Drawer_Sphere drawer_sphere;
  unsigned int id_tex;
  std::tuple<int,int,int> texInfo;
  dfm2::opengl::Drawer_RectangleTex drawer_rect;
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
    /*
    if (ImGui::Button("Load Texture")) {
      ImGuiFileDialog::Instance()->OpenDialog(
          "ChooseTexFileDlgKey",
          "Choose File",
          "{.png,.jpg}",
          ".");
    }
    if (ImGuiFileDialog::Instance()->Display("ChooseTexFileDlgKey")) {
      if (ImGuiFileDialog::Instance()->IsOk()) {
        std::string filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
        viewer.SetTextureImageFile(filePathName);
        std::cout << "texPath: " << filePathName << std::endl;
      }
      ImGuiFileDialog::Instance()->Close();
    }
     */
    ImGui::End();
  }

  viewer.Draw();

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
