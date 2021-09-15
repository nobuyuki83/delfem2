/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <vector>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/srchbi_v3bvh.h"
#include "delfem2/srchbv3sphere.h"
#include "delfem2/srchbvh.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/points.h"
#include "delfem2/mshuni.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// ------------------------------------

int main() {
  // shape definition
  std::vector<double> vtx_xyz_ini; // undeformed vertex positions
  std::vector<double> vtx_xyz; // deformed vertex positions
  std::vector<double> vtx_uvw; // vertex deformation modes
  std::vector<unsigned int> tri_vtx;  // index of triangles

  // variables for self-collision
  int root_bvh_idx; // index BVH root node
  std::vector<dfm2::CNodeBVH2> bvh_nodes; // array of BVH node
  std::vector<dfm2::CBV3d_Sphere> bounding_volumes; // array of AABB same size as aNodeBVH
  std::vector<dfm2::CIntersectTriPair<double>> intersected_triangle_pairs;

  // data for camera
  double cur_time = 0;
  {
    delfem2::MeshTri3D_Sphere(
        vtx_xyz_ini, tri_vtx,
        1.0, 16, 16);
    delfem2::Rotate_Points3(
        vtx_xyz_ini,
        0.2, 0.3, 0.4);
    vtx_xyz = vtx_xyz_ini;
    {
      const size_t ntri = tri_vtx.size() / 3;
      std::vector<double> aElemCenter(ntri * 3);
      for (unsigned int itri = 0; itri < ntri; ++itri) {
        dfm2::CVec3d p0 = dfm2::CG_Tri3(itri, tri_vtx, vtx_xyz);
        aElemCenter[itri * 3 + 0] = p0.x;
        aElemCenter[itri * 3 + 1] = p0.y;
        aElemCenter[itri * 3 + 2] = p0.z;
      }
      std::vector<unsigned int> tri_adjtri;
      dfm2::ElSuEl_MeshElem(
          tri_adjtri,
          tri_vtx.data(), tri_vtx.size() / 3,
          delfem2::MESHELEM_TRI, vtx_xyz.size() / 3);
      root_bvh_idx = dfm2::BVHTopology_TopDown_MeshElem(
          bvh_nodes,
          3, tri_adjtri,
          aElemCenter);
      std::cout << "aNodeBVH.size(): " << bvh_nodes.size() << std::endl;
    }
  }
  {
    vtx_uvw.assign(vtx_xyz.size(), 0.0);
    for (size_t ixyz = 0; ixyz < vtx_xyz.size() / 3; ++ixyz) {
      double x0 = vtx_xyz[ixyz * 3 + 0];
      vtx_uvw[ixyz * 3 + 0] = -3 * x0 * x0 * x0 * x0 * x0;
    }
  }

  dfm2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  delfem2::opengl::setSomeLighting();

  while (!glfwWindowShouldClose(viewer.window)) {
    {
      cur_time += 0.02;
      double d = sin(cur_time);
      for (unsigned int ip = 0; ip < vtx_xyz.size() / 3; ip++) {
        vtx_xyz[ip * 3 + 0] = vtx_xyz_ini[ip * 3 + 0] + vtx_uvw[ip * 3 + 0] * d;
        vtx_xyz[ip * 3 + 1] = vtx_xyz_ini[ip * 3 + 1] + vtx_uvw[ip * 3 + 1] * d;
        vtx_xyz[ip * 3 + 2] = vtx_xyz_ini[ip * 3 + 2] + vtx_uvw[ip * 3 + 2] * d;
      }
      dfm2::CLeafVolumeMaker_Mesh<dfm2::CBV3d_Sphere, double> lvm(
          1.0e-5,
          vtx_xyz.data(), vtx_xyz.size() / 3,
          tri_vtx.data(), tri_vtx.size() / 3, 3);
      dfm2::BVH_BuildBVHGeometry(
          bounding_volumes,
          root_bvh_idx, bvh_nodes,
          lvm);
      intersected_triangle_pairs.clear();
      dfm2::GetIntersectTriPairs(
          intersected_triangle_pairs,
          vtx_xyz, tri_vtx,
          root_bvh_idx,
          bvh_nodes, bounding_volumes); // output
      std::cout << intersected_triangle_pairs.size() << std::endl;
    }
    viewer.DrawBegin_oldGL();
    //  Draw_SurfaceMeshNorm(aXYZ, aTri, aNormal);
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0, 0, 0);
    delfem2::opengl::DrawMeshTri3D_Edge(vtx_xyz, tri_vtx);
    ::glDisable(GL_LIGHTING);
    ::glLineWidth(2);
    ::glColor3d(1, 0, 0);
    ::glBegin(GL_LINES);
    for (const auto &tp : intersected_triangle_pairs) {
      glVertex3d(tp.P[0].x, tp.P[0].y, tp.P[0].z);
      glVertex3d(tp.P[1].x, tp.P[1].y, tp.P[1].z);
    }
    ::glEnd();
    viewer.SwapBuffers();
    glfwPollEvents();
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


