//
//  glfw_viewer.hpp
//  00_openwin
//
//  Created by Nobuyuki Umetani on 2019-10-31.
//

#ifndef glfw_viewer_hpp
#define glfw_viewer_hpp

#include <stdio.h>

#include "delfem2/opengl/glfw_cam.h"

// -----------

class CViewer_GLFW{
public:
  void Init_GLold();
  void Init_GLnew();
  
  void DrawBegin_Glold();
  void DrawEnd_oldGL();
  
public:
  GLFWwindow* window;
  CNav3D_GLFW nav;
};


#endif /* glfw_viewer_hpp */
