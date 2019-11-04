#include <iostream>
#include <sstream>
#include <math.h>
#include <complex>
#include <set>
#include <stack>
#include "delfem2/dtri_v2.h"
#include "delfem2/cad2d.h"

// --------------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/gl2_funcs.h"
#include "delfem2/opengl/gl2_v23.h"
#include "delfem2/opengl/gl2_v23dtricad.h"

#ifndef M_PI
  #define M_PI 3.141592653589793
#endif

// ----------------------------------------------

int main(int argc,char* argv[])
{
  class CCadMesh2DVeiwer : public CViewer_GLFW
  {
  public:
    CCadMesh2DVeiwer(){
      const double poly[8] = {-1,-1, +1,-1, +1,+1, -1,+1};
      cad.AddPolygon(std::vector<double>(poly,poly+8));
      cad.is_draw_face = false;
      {
        std::vector<int> aFlgPnt, aFlgTri;
        CMeshDynTri2D dmsh;
        CMesher_Cad2D mesher;
        mesher.edge_length = 0.5;
        mesher.Meshing(dmsh, cad);
        dmsh.Export_StlVectors(aXY, aTri);
      }
      
      const int nv = 4;
      std::vector<double> aXY_bound = cad.XY_VtxCtrl_Face(0);
      const int nXY = aXY.size()/2;
      aW.resize(nXY*nv);
      for(int ip=0;ip<nXY;++ip){
        MeanValueCoordinate2D(aW.data()+nv*ip,
                              aXY[ip*2+0], aXY[ip*2+1],
                              aXY_bound.data(), aXY_bound.size()/2);
        double sum = 0.0;
        for(int ipb=0;ipb<aXY_bound.size()/2;++ipb){
          sum += aW[nv*ip+ipb];
        }
        assert( fabs(sum-1)<1.0e-10 );
      }
    }
    virtual void mouse_press(const float src[3], const float dir[3]){
      cad.Pick(src[0], src[1], this->nav.camera.view_height);
    }
    virtual void mouse_drag(const float src0[3], const float src1[3], const float dir[3]){
      cad.DragPicked(src1[0], src1[1], src0[0], src0[1]);
      std::vector<double> aXY_bound = cad.XY_VtxCtrl_Face(0);
      int npb = aXY_bound.size()/2;
      int np = aXY.size()/2;
      for(int ip=0;ip<np;++ip){
        aXY[ip*2+0] = 0.0;
        aXY[ip*2+1] = 0.0;
        for(int ipb=0;ipb<npb;++ipb){
          aXY[ip*+2+0] += aW[ip*npb+ipb]*aXY_bound[ipb*2+0];
          aXY[ip*+2+1] += aW[ip*npb+ipb]*aXY_bound[ipb*2+1];
        }
      }
    }
    void Draw(){
      this->DrawBegin_Glold();
      Draw_CCad2D(cad);
      ::glDisable(GL_LIGHTING);
      ::glColor3d(0.8, 0.8, 0.8);
      opengl::DrawMeshTri2D_Face(aTri,aXY);
      ::glLineWidth(1);
      opengl::DrawMeshTri2D_Edge(aTri,aXY);
      this->DrawEnd_oldGL();
    }
  public:
    std::vector<double> aXY;
    std::vector<unsigned int> aTri;
    std::vector<double> aW;
    CCad2D cad;
  };
  
  // ----------------------------------
  CCadMesh2DVeiwer viewer;
  viewer.Init_GLold();
  viewer.nav.camera.view_height = 1.5;
  viewer.nav.camera.camera_rot_mode = CAMERA_ROT_YTOP;
  opengl::setSomeLighting();
  while(!glfwWindowShouldClose(viewer.window)){
    viewer.Draw();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


