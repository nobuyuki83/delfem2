/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <stack>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/cagedef.h"
#include "delfem2/cad2_dtri2.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/v2.h"
#include "delfem2/opengl/old/cad2dtriv2.h"

#ifndef M_PI
#  define M_PI 3.141592653589793
#endif

namespace dfm2 = delfem2;

// ----------------------------------------------

int main()
{
  class CCadMesh2DVeiwer : public delfem2::glfw::CViewer3
  {
  public:
    CCadMesh2DVeiwer(){
      const double poly[8] = {-1,-1, +1,-1, +1,+1, -1,+1};
      cad.AddPolygon(std::vector<double>(poly,poly+8));
      cad.is_draw_face = false;
      {
        std::vector<int> aFlgPnt, aFlgTri;
        delfem2::CMeshDynTri2D dmsh;
        delfem2::CMesher_Cad2D mesher;
        mesher.edge_length = 0.5;
        mesher.Meshing(dmsh, cad);
        dmsh.Export_StlVectors(aXY, aTri);
      }
      
      const int nv = 4;
      std::vector<double> aXY_bound = cad.XY_VtxCtrl_Face(0);
      const auto nXY = static_cast<unsigned int>(aXY.size()/2);
      aW.resize(nXY*nv);
      for(unsigned int ip=0;ip<nXY;++ip){
        dfm2::MeanValueCoordinate_Polygon2<dfm2::CVec2d>(
            aW.data()+nv*ip,
            aXY[ip*2+0], aXY[ip*2+1],
            aXY_bound.data(), static_cast<unsigned int>(aXY_bound.size()/2));
        double sum = 0.0;
        for(unsigned int ipb=0;ipb<aXY_bound.size()/2;++ipb){
          sum += aW[nv*ip+ipb];
        }
        assert( fabs(sum-1)<1.0e-10 );
      }
    }
    void mouse_press(
		const float src[3], 
		[[maybe_unused]] const float dir[3]) override{
      cad.Pick(src[0], src[1], this->camera.view_height);
    }
    void mouse_drag(const float src0[3], const float src1[3], 
		[[maybe_unused]] const float dir[3]) override{
      cad.DragPicked(src1[0], src1[1], src0[0], src0[1]);
      std::vector<double> aXY_bound = cad.XY_VtxCtrl_Face(0);
      auto npb = static_cast<unsigned int>(aXY_bound.size()/2);
      auto np = static_cast<unsigned int>(aXY.size()/2);
      for(unsigned int ip=0;ip<np;++ip){
        aXY[ip*2+0] = 0.0;
        aXY[ip*2+1] = 0.0;
        for(unsigned int ipb=0;ipb<npb;++ipb){
          aXY[ip*+2+0] += aW[ip*npb+ipb]*aXY_bound[ipb*2+0];
          aXY[ip*+2+1] += aW[ip*npb+ipb]*aXY_bound[ipb*2+1];
        }
      }
    }
    void Draw(){
      this->DrawBegin_oldGL();
      delfem2::opengl::Draw_CCad2D(cad);
      ::glDisable(GL_LIGHTING);
      ::glColor3d(0.8, 0.8, 0.8);
      delfem2::opengl::DrawMeshTri2D_Face(aTri,aXY);
      ::glLineWidth(1);
      delfem2::opengl::DrawMeshTri2D_Edge(aTri,aXY);
      this->SwapBuffers();
    }
  public:
    std::vector<double> aXY;
    std::vector<unsigned int> aTri;
    std::vector<double> aW;
    delfem2::CCad2D cad;
  };
  
  // ----------------------------------
  CCadMesh2DVeiwer viewer;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.view_height = 1.5;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::YTOP;
  delfem2::opengl::setSomeLighting();
  while(!glfwWindowShouldClose(viewer.window)){
    viewer.Draw();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


