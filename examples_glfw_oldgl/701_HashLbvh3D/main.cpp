#include <vector>
#include <algorithm>
#include <random>
#include "delfem2/bv.h"
#include "delfem2/bvh.h"

// ---------------------------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/glold_funcs.h"

namespace dfm2 = delfem2;

// ------------------------------------
// input parameter for simulation
std::vector<double> aXYZ; // 3d points
std::vector<dfm2::CNodeBVH2> aNodeBVH;

// ----------------------------------------

void myGlutDisplay()
{
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  ::glPointSize(3);
  ::glBegin(GL_POINTS);
  for(size_t ip=0;ip<aXYZ.size()/3;++ip){
    ::glVertex3d(aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]);
  }
  ::glEnd();
}

int main(int argc,char* argv[])
{
  {
    const double min_xyz[3] = {-1,-1,-1};
    const double max_xyz[3] = {+1,+1,+1};
    dfm2::CBV3D_AABB bb(min_xyz, max_xyz);
    {
      const unsigned int N = 10000;
      aXYZ.resize(N*3);
      std::random_device dev;
      std::mt19937 rng(dev());
      std::uniform_real_distribution<> udist(0.0, 1.0);
      for(unsigned int i=0;i<N;++i) {
        aXYZ[i * 3 + 0] = (bb.x_max - bb.x_min) * udist(rng) + bb.x_min;
        aXYZ[i * 3 + 1] = (bb.y_max - bb.y_min) * udist(rng) + bb.y_min;
        aXYZ[i * 3 + 2] = (bb.z_max - bb.z_min) * udist(rng) + bb.z_min;
      }
      srand(3);
      for(int iip=0;iip<10;++iip){ // hash collision
        const unsigned int ip = N*(rand()/(RAND_MAX+1.0));
        assert( N >= 0 && ip < N);
        double x0 = aXYZ[ip*3+0];
        double y0 = aXYZ[ip*3+1];
        double z0 = aXYZ[ip*3+2];
        for(int itr=0;itr<2;itr++){
          aXYZ.insert(aXYZ.begin(), z0);
          aXYZ.insert(aXYZ.begin(), y0);
          aXYZ.insert(aXYZ.begin(), x0);
        }
      }
    }
    std::vector<unsigned int> aSortedId;
    std::vector<std::uint32_t> aSortedMc;
    dfm2::GetSortedMortenCode(aSortedId,aSortedMc,
                              aXYZ,min_xyz,max_xyz);
    {
      const double bbmin[3] = {bb.x_min, bb.y_min, bb.z_min};
      const double bbmax[3] = {bb.x_max, bb.y_max, bb.z_max};
      dfm2::Check_MortonCode_Sort(aSortedId, aSortedMc, aXYZ, bbmin, bbmax);
    }
    dfm2::Check_MortonCode_RangeSplit(aSortedMc);
    dfm2::BVH_TreeTopology_Morton(aNodeBVH,
                                  aSortedId,aSortedMc);
    dfm2::Check_BVH(aNodeBVH,aXYZ.size()/3);
  }
  
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.5;
  delfem2::opengl::setSomeLighting();
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


