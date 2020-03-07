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
#include "delfem2/vecxitrsol.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_funcs.h"

namespace dfm2 = delfem2;

int main()
{
  std::vector<double> aXYZ0;
  std::vector<double> aW;
  std::vector<unsigned int> aTri;
  std::vector<int> aIndBoneParent;
  std::vector<double> aJntRgrs;
  {
    cnpy::npz_t my_npz = cnpy::npz_load(std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    const unsigned int nP = my_npz["vertices_template"].shape[0];
    assert( my_npz["vertices_template"].shape[1] == 3 );
    aXYZ0 = my_npz["vertices_template"].as_vec<double>();
    {
      cnpy::NpyArray& npT = my_npz["face_indices"];
      assert( my_npz["face_indices"].shape[1] == 3 );
      aTri = npT.as_vec<unsigned>();
      for(auto &i: aTri){ i -= 1; }
    }
    const unsigned int nBone = my_npz["weights"].shape[1];
    assert( my_npz["weights"].shape[0] == nP );
    aW = my_npz["weights"].as_vec<double>();
    {
      const cnpy::NpyArray& npT = my_npz["kinematic_tree"];
      const int* tree = npT.data<int>();
      aIndBoneParent.assign(tree, tree+nBone);
    }
    assert( my_npz["joint_regressor"].shape[0] == nBone );
    assert( my_npz["joint_regressor"].shape[1] == nP );
    std::cout << my_npz["joint_regressor"].fortran_order << std::endl;
    aJntRgrs = my_npz["joint_regressor"].as_vec<double>();
    assert( aJntRgrs.size() == nBone*nP );
  }
    
  std::vector<double> aJntPos;
  {
    const unsigned int nP = aXYZ0.size()/3;
    const unsigned int nBone = aIndBoneParent.size();
    aJntPos.assign(nBone*3, 0.0);
    for(int ib=0;ib<nBone;++ib){
      aJntPos[ib*3+0] = 0;
      aJntPos[ib*3+1] = 0;
      aJntPos[ib*3+2] = 0;
      for(int ip=0;ip<nP;++ip){
        aJntPos[ib*3+0] += aJntRgrs[ip*nBone+ib]*aXYZ0[ip*3+0];
        aJntPos[ib*3+1] += aJntRgrs[ip*nBone+ib]*aXYZ0[ip*3+1];
        aJntPos[ib*3+2] += aJntRgrs[ip*nBone+ib]*aXYZ0[ip*3+2];
      }
      std::cout <<  ib << " ";
      std::cout << aJntPos[ib*3+0] << " ";
      std::cout << aJntPos[ib*3+1] << " ";
      std::cout << aJntPos[ib*3+2] << std::endl;
    }
  }
  
  
  // -----------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  dfm2::opengl::setSomeLighting();

  while (true)
  {
    viewer.DrawBegin_oldGL();
    
    ::glEnable(GL_LIGHTING);
    ::glEnable(GL_DEPTH_TEST);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ0.data(), aTri.data(), aTri.size()/3);
    for(int ib=0;ib<aJntPos.size()/3;++ib){
      const double* p = aJntPos.data()+ib*3;
      ::glColor3d(1,0,0);
      ::glDisable(GL_LIGHTING);
      ::glDisable(GL_DEPTH_TEST);
      dfm2::opengl::DrawSphereAt(8, 8, 0.01, p[0], p[1], p[2]);
      int ibp = aIndBoneParent[ib];
      if( ibp == -1 ){ continue; }
      const double* pp = aJntPos.data()+ibp*3;
      ::glColor3d(0,0,0);
      ::glBegin(GL_LINES);
      ::glVertex3dv(p);
      ::glVertex3dv(pp);
      ::glEnd();
    }
    
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
