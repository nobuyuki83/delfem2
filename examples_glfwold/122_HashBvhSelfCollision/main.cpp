/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <vector>
#include "delfem2/bv.h"
#include "delfem2/bvh.h"
#include "delfem2/primitive.h"
#include "delfem2/points.h"
#include "delfem2/mshtopo.h"
#include "delfem2/srchbi_v3bvh.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// ------------------------------------

int main(int argc,char* argv[])
{
  // shape definition
  std::vector<double> aXYZ0; // undeformed vertex positions
  std::vector<double> aXYZ; // deformed vertex positions
  std::vector<double> aUVW; // vertex deformation modes
  std::vector<unsigned int> aTri;  // index of triangles
  
  // variables for self-collision
  int iroot_bvh; // index BVH root node
  std::vector<dfm2::CNodeBVH2> aNodeBVH; // array of BVH node
  std::vector<dfm2::CBV3d_Sphere> aBB_BVH; // array of AABB same size as aNodeBVH
  std::vector<dfm2::CIntersectTriPair<double>> aITP;
  
  // data for camera
  double cur_time = 0;
  {
    delfem2::MeshTri3D_Sphere(aXYZ0, aTri, 1.0, 16, 16);
    delfem2::Rotate_Points3(aXYZ0,
                            0.2, 0.3, 0.4);
    aXYZ = aXYZ0;
    {
      const size_t ntri = aTri.size()/3;
      std::vector<double> aElemCenter(ntri*3);
      for(unsigned int itri=0;itri<ntri;++itri){
        dfm2::CVec3d p0 = dfm2::CG_Tri3(itri, aTri, aXYZ);
        aElemCenter[itri*3+0] = p0.x();
        aElemCenter[itri*3+1] = p0.y();
        aElemCenter[itri*3+2] = p0.z();
      }
      std::vector<unsigned int> aTriSuTri;
      dfm2::ElSuEl_MeshElem(
          aTriSuTri,
          aTri.data(), aTri.size()/3,
          delfem2::MESHELEM_TRI, aXYZ.size()/3);
      iroot_bvh = dfm2::BVHTopology_TopDown_MeshElem(
          aNodeBVH,
          3,aTriSuTri,
          aElemCenter);
      std::cout << "aNodeBVH.size(): " << aNodeBVH.size() << std::endl;
    }
  }
  {
    aUVW.assign(aXYZ.size(),0.0);
    for(size_t ixyz=0;ixyz<aXYZ.size()/3;++ixyz){
      double x0 = aXYZ[ixyz*3+0];
      aUVW[ixyz*3+0] = -3*x0*x0*x0*x0*x0;
    }
  }
  
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  delfem2::opengl::setSomeLighting();
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    {
      cur_time += 0.02;
      double d = sin(cur_time);
      for(int ip=0;ip<(int)aXYZ.size()/3;ip++){
        aXYZ[ip*3+0] =  aXYZ0[ip*3+0] + aUVW[ip*3+0]*d;
        aXYZ[ip*3+1] =  aXYZ0[ip*3+1] + aUVW[ip*3+1]*d;
        aXYZ[ip*3+2] =  aXYZ0[ip*3+2] + aUVW[ip*3+2]*d;
      }
      dfm2::BVH_BuildBVHGeometry_Mesh(
                                      aBB_BVH,
                                      iroot_bvh,aNodeBVH,
                                      1.0e-5,
                                      aXYZ.data(),aXYZ.size()/3,
                                      aTri.data(),3,aTri.size()/3);
      aITP.clear();
      dfm2::GetIntersectTriPairs(aITP,
                                 aXYZ,aTri,
                                 iroot_bvh,
                                 aNodeBVH,aBB_BVH); // output
      std::cout << aITP.size() << std::endl;
    }
    viewer.DrawBegin_oldGL();
    //  Draw_SurfaceMeshNorm(aXYZ, aTri, aNormal);
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0,0,0);
    delfem2::opengl::DrawMeshTri3D_Edge(aXYZ,aTri);
    ::glDisable(GL_LIGHTING);
    ::glLineWidth(2);
    ::glColor3d(1,0,0);
    ::glBegin(GL_LINES);
    for(const auto & tp : aITP){
      glVertex3d(tp.P[0].x(), tp.P[0].y(), tp.P[0].z());
      glVertex3d(tp.P[1].x(), tp.P[1].y(), tp.P[1].z());
    }
    ::glEnd();
    viewer.DrawEnd_oldGL();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


