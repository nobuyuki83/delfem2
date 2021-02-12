
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
//
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "delfem2/opengl/funcs.h"
#include <cstdio>
#include <iostream>
#include <vector>
#include <GLFW/glfw3.h>

namespace dfm2 = delfem2;

class CShader
{
public:
  void Compile()
  {
    const std::string glsl33vert_projection =
    "layout (location = 0) in vec3 posIn;\n"
    "layout (location = 1) in vec3 colorIn;\n"
    "out vec3 color;\n"
    "void main()\n"
    "{\n"
    "  gl_Position = vec4(posIn.x, posIn.y, posIn.z, 1.0);\n"
    "  color = colorIn;\n"
    "}\0";
    
    const std::string glsl33frag =
    "out vec4 FragColor;\n"
    "in vec3 color;\n"
    "void main()\n"
    "{\n"
    "  FragColor = vec4(color.x, color.y, color.z, 1.0f);\n"
    "}\n\0";
    
#ifdef EMSCRIPTEN
    shaderProgram = dfm2::opengl::GL24_CompileShader((std::string("#version 300 es\n")+
                                                      glsl33vert_projection).c_str(),
                                                     (std::string("#version 300 es\n")+
                                                      std::string("precision highp float;\n")+
                                                      glsl33frag).c_str());
#else
    shaderProgram = dfm2::opengl::GL24_CompileShader((std::string("#version 330 core\n")+
                                                      glsl33vert_projection).c_str(),
                                                     (std::string("#version 330 core\n")+
                                                      glsl33frag).c_str());
#endif

    if( !glIsProgram(shaderProgram) ){
      std::cout << "shader doesnot exist" << std::endl;
    }
    glUseProgram(shaderProgram);
  }
  void create_triangle()
  {
      // create the triangle
    float triangle_vertices[] = {
      0.0f, 0.25f, 0.0f,  // position vertex 1
      1.0f, 0.0f, 0.0f,   // color vertex 1
      0.25f, -0.25f, 0.0f,  // position vertex 1
      0.0f, 1.0f, 0.0f,   // color vertex 1
      -0.25f, -0.25f, 0.0f, // position vertex 1
      0.0f, 0.0f, 1.0f,   // color vertex 1
    };
    unsigned int triangle_indices[] = {
      0, 1, 2};
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glGenBuffers(1, &ebo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(triangle_vertices), triangle_vertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(triangle_indices), triangle_indices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
  }
  void Draw(){ // rendering our geometries
    glUseProgram( shaderProgram );
    glBindVertexArray( vao );
    glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
  }
public:
//  CGL4_VAO_Mesh vao;
  int shaderProgram;
  unsigned int vbo;
  unsigned int vao;
  unsigned int ebo;
};

static void callback_error(int error, const char* description)
{
    fputs(description, stderr);
}

static GLFWwindow* myGLFW_OpenWindow
        (const unsigned int SCR_WIDTH,
         const unsigned int SCR_HEIGHT)
{
    glfwSetErrorCallback(callback_error);
    if (!glfwInit()){
        exit(EXIT_FAILURE);
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this statement to fix compilation on OS X
#endif
/*
    // Decide GL+GLSL versions
#if __APPLE__
    // GL 3.2 + GLSL 150
  const char *glsl_version = "#version 150";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // 3.2+ only
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);       // Required on Mac
#else
    // GL 3.0 + GLSL 130
  const char *glsl_version = "#version 130";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
#endif
*/
// glfw window creation
// --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH,
                                          SCR_HEIGHT,
                                          "LearnOpenGL",
                                          NULL,
                                          NULL);
    if (window == NULL)
    {
//    std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return 0;
    }
    return window;
}

// ---------------------------------------------------------------

CShader shdr0;

void draw(GLFWwindow* window)
{
  glfwPollEvents();
  glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
  glClear(GL_COLOR_BUFFER_BIT);
  
    // feed inputs to dear imgui, start new frame
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();
  
  { // render your GUI
    ImGui::Begin("Triangle Position/Color");
    static float rotation = 0.0;
    ImGui::SliderFloat("rotation", &rotation, 0, 2 * 3.1415);
    static float translation[] = {0.0, 0.0};
    ImGui::SliderFloat2("position", translation, -1.0, 1.0);
    static float color[4] = { 1.0f,1.0f,1.0f,1.0f };
      // pass the parameters to the shader
      //      triangle_shader.setUniform("rotation", rotation);
      //      triangle_shader.setUniform("translation", translation[0], translation[1]);
      // color picker
    ImGui::ColorEdit3("color", color);
      // multiply triangle's color with this color
      //      triangle_shader.setUniform("color", color[0], color[1], color[2]);
    ImGui::End();
  }
  
  shdr0.Draw();
  
    // Render dear imgui into screen
  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  
  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  glfwSwapBuffers(window);
}

int main(int, char **)
{
  GLFWwindow* window = myGLFW_OpenWindow(800,600);
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1); // Enable vsync
  
#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif

	int screen_width, screen_height;
	glfwGetFramebufferSize(window, &screen_width, &screen_height);
	glViewport(0, 0, screen_width, screen_height);
  
  shdr0.Compile();
  shdr0.create_triangle();
	
  // Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGui_ImplGlfw_InitForOpenGL(window, true);// Setup Platform/Renderer bindings
  ImGui::StyleColorsDark(); // Setup Dear ImGui style
    
#ifdef EMSCRIPTEN
  ImGui_ImplOpenGL3_Init("#version 100");
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, window, 60, 1);
#else
  ImGui_ImplOpenGL3_Init("#version 150");
  while (!glfwWindowShouldClose(window)) { draw(window); }
#endif

	// Cleanup
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}
