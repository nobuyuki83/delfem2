#include <iostream>
#include <math.h>

#if defined(_MSC_VER)
  #include <windows.h>
#endif

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#endif

#include "delfem2/cad2d.h"

#include "delfem2/gl24_funcs.h"
#include "delfem2/gl4_funcs.h"
#include "../glfw_funcs.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

CNav3D_GLFW nav;
CCad2D cad;
int shaderProgram;
int Loc_MatrixProjection;
int Loc_MatrixModelView;
int Loc_Color;
CGL4_VAO_Mesh vao_face;


void draw(GLFWwindow* window)
{
  float asp;
  {
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    glViewport(0, 0, width, height);
    asp = width / (float) height;
//    std::cout << width << " " << height << " " << asp << std::endl;
  }
  
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );

  
  glUseProgram(shaderProgram);
  float mP[16]; nav.camera.Affine4f_Projection(mP, asp, 10);
  glUniformMatrix4fv(Loc_MatrixProjection, 1, GL_FALSE, mP);
  float mMV[16]; nav.camera.Affine4f_ModelView(mMV);
  glUniformMatrix4fv(Loc_MatrixModelView, 1, GL_FALSE, mMV);
  glUniform3f(Loc_Color, 1,1,1);
  vao_face.Draw(0);
  glUniform3f(Loc_Color, 0,0,0);
  vao_face.Draw(1);
  
  
  glfwSwapBuffers(window);
  glfwPollEvents();
}

void callback_key(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS){
    glfwSetWindowShouldClose(window, GL_TRUE);
  }
}

void callback_resize(GLFWwindow* window, int width, int height)
{
  glViewport(0, 0, width, height);
}

void callback_mouse_button(GLFWwindow* window, int button, int action, int mods)
{
  nav.Mouse(window,button,action,mods);

}

void callback_cursor_position(GLFWwindow* window, double xpos, double ypos)
{
  nav.Motion(window,xpos,ypos);
}

void callback_scroll(GLFWwindow* window, double xoffset, double yoffset)
{
  nav.camera.scale *= pow(1.01,yoffset);
}


int main(void)
{
  GLFWwindow* window = myGLFW_OpenWindow(800,600);
  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, callback_resize);
  glfwSetKeyCallback(            window, callback_key);
  glfwSetMouseButtonCallback(    window, callback_mouse_button);
  glfwSetCursorPosCallback(      window, callback_cursor_position);
  glfwSetScrollCallback(         window, callback_scroll);
  
  // glad: load all OpenGL function pointers
  // ---------------------------------------
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

  {
    std::vector<double> aXY = {-1,-1, +1,-1, +1,+1, -1,+1};
    cad.AddPolygon(aXY);
    std::vector<float> aXYf;
    std::vector<unsigned int> aLine;
    std::vector<unsigned int> aTri;
    {
      const std::vector<CVector2>& aVec2 = cad.aVec2_Tessellation;
      aXYf.reserve(aVec2.size()*2);
      for(int iv=0;iv<aVec2.size();++iv){
        aXYf.push_back(aVec2[iv].x);
        aXYf.push_back(aVec2[iv].y);
      }
      for(int ie=0;ie<cad.aEdge.size();++ie){
        const CCadTopo& topo = cad.topo;
        int iv0 = topo.aEdge[ie].iv0;
        int iv1 = topo.aEdge[ie].iv1;
        int ip0 = cad.aEdge[ie].ip0;
        int nseg = cad.aEdge[ie].aP.size()+1;
        for(int iseg=0;iseg<nseg;++iseg){
          int i0 = (iseg==0)?iv0:ip0+iseg-1;
          int i1 = (iseg==nseg-1)?iv1:iseg;
          aLine.push_back(i0);
          aLine.push_back(i1);
          std::cout << i0 << " " << i1 << std::endl;
        }
      }
      for(int ifc=0;ifc<cad.aFace.size();++ifc){
        const CCad2D_FaceGeo& fc = cad.aFace[ifc];
        aTri.insert(aTri.end(),fc.aTri.begin(),fc.aTri.end());
      }
    }
    {
      unsigned int VAO;
      unsigned int VBO_pos;
      GL4_VAO_MeshTri3D(VAO,VBO_pos,
                        aXYf.data(),aXYf.size()/2,2);
      unsigned int EBO_Tri;
      {
        glBindVertexArray(VAO); // opengl4
        glGenBuffers(1, &EBO_Tri);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_Tri);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*aTri.size(), aTri.data(), GL_STATIC_DRAW);
      }
      unsigned int EBO_Line;
      {
        glBindVertexArray(VAO); // opengl4
        glGenBuffers(1, &EBO_Line);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_Line);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*aLine.size(), aLine.data(), GL_STATIC_DRAW);
      }
      std::cout << "VAO: " << VAO << std::endl;
      vao_face.VAO = VAO;
      {
        CGL4_VAO_Mesh::CElem e0;
        e0.size = aTri.size();
        e0.GL_MODE = GL_TRIANGLES;
        e0.EBO = EBO_Tri;
        vao_face.aElem.push_back(e0);
      }
      {
        CGL4_VAO_Mesh::CElem e0;
        e0.size = aLine.size();
        e0.GL_MODE = GL_LINES;
        e0.EBO = EBO_Line;
        vao_face.aElem.push_back(e0);
      }
    }
    /*
    {
      int VAO = GL4_VAO_MeshTri3D(aXYf.data(),aXYf.size()/2,2,
                                  aLine.data(),aLine.size()/2);
      std::cout << "VAO: " << VAO << std::endl;
      vao_edge.VAO = VAO;
      vao_edge.nElem = aLine.size()/2;
      vao_edge.nNoel = 2;
      vao_edge.GL_ELEM_TYPE = GL_LINES;
    }
     */
  }

#ifdef EMSCRIPTEN
  shaderProgram = GL24_CompileShader((std::string("#version 300 es\n")+
                                      glsl33vert_projection).c_str(),
                                     (std::string("#version 300 es\n")+
                                      std::string("precision highp float;\n")+
                                      glsl33frag).c_str());
#else
  shaderProgram = GL24_CompileShader((std::string("#version 330 core\n")+
                                      glsl33vert_projection).c_str(),
                                     (std::string("#version 330 core\n")+
                                      glsl33frag).c_str());
#endif
  
  
  {
    if( !glIsProgram(shaderProgram) ){
      std::cout << "shader doesnot exist" << std::endl;
    }
    glUseProgram(shaderProgram);
    Loc_MatrixProjection = glGetUniformLocation(shaderProgram,  "matrixProjection");
    Loc_MatrixModelView  = glGetUniformLocation(shaderProgram,  "matrixModelView");
    Loc_Color            = glGetUniformLocation(shaderProgram,  "color");
    std::cout << "projectionMatrixLoc: " << Loc_MatrixProjection << "   shaderProgram: " << shaderProgram << "  LocColor: " << Loc_Color << std::endl;
  }
  
  nav.camera.view_height = 1.5;
  nav.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  
  
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, window, 60, 1);
#else
  while (!glfwWindowShouldClose(window)) { draw(window); }
#endif
  
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

