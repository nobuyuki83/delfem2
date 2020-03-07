/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief SMPL model
 * @details hogehoge
 */

#include <cstdlib>

#include "cnpy/cnpy.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_funcs.h"

namespace dfm2 = delfem2;

int main()
{
  std::vector<double> aXYZ0;
  std::vector<double> aW;
  std::vector<unsigned int> aTri;
  {
    cnpy::npz_t my_npz = cnpy::npz_load(std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    {
      cnpy::NpyArray& npT = my_npz["face_indices"];
      aTri = npT.as_vec<unsigned>();
      for(auto &i: aTri){ i -= 1; }
    }
    aXYZ0 = my_npz["vertices_template"].as_vec<double>();
    aW = my_npz["weights"].as_vec<double>();
  }
        
  // -----------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  dfm2::opengl::setSomeLighting();

  while (true)
  {
    viewer.DrawBegin_oldGL();
    
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ0.data(), aTri.data(), aTri.size()/3);
    
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
