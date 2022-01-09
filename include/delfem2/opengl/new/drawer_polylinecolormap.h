//
// Created by Nobuyuki Umetani on 2021-10-08.
//

#ifndef DRAWER_POLYLINECOLORMAP_H_
#define DRAWER_POLYLINECOLORMAP_H_

#include <array>

#include "delfem2/opengl/new/drawer_mshcolormap.h"
#include "delfem2/vec3.h"
#include "delfem2/geo3_v23m34q.h"

namespace delfem2::opengl {

class Drawer_PolylineColormap {
 public:
  void InitGL() {
    sphere.InitGL();
    {
      std::vector<float> vtx_xyz;
      std::vector<unsigned int> tri_vtx;
      delfem2::MeshTri3D_Sphere(
          vtx_xyz, tri_vtx,
          1.f, 32, 32);
      sphere.AddConnectivity(tri_vtx, GL_TRIANGLES);
      sphere.SetCoordinates(vtx_xyz, 3);
      std::vector<float> vtx_val(vtx_xyz.size() / 3, 1.f);
      sphere.SetValues(vtx_val);
    }
    cylinder.InitGL();
    {
      std::vector<float> vtx_xyz;
      std::vector<unsigned int> tri_vtx;
      delfem2::MeshTri3D_CylinderOpen(
          vtx_xyz, tri_vtx,
          1.f, 1.f, 32, 1);
      cylinder.AddConnectivity(tri_vtx, GL_TRIANGLES);
      cylinder.SetCoordinates(vtx_xyz, 3);
      std::vector<float> vtx_val(vtx_xyz.size() / 3);
      for(unsigned int iv=0;iv<vtx_xyz.size()/3;++iv){
        vtx_val[iv] = (vtx_xyz[iv*3+1]+0.5);
      }
      cylinder.SetValues(vtx_val);
    }
  }

  /**
   * @detail use delfem2::ViewAsArray2D in "delfem2/view_array2d.h" for row pointer array
   * @tparam ARRAY2D has size() and [][] e.g., std::vector<dfm2::CVecX>, std::vector<Eigen::VectorX>
   * @param vtx_pos array of coordinate
   * @param ndim dimension of the points
   * @param mat_projection
   * @param mat_modelview
   */
  template<typename ARRAY2D>
  void Draw(
      const ARRAY2D &vtx_pos,
      unsigned int ndim,
      const float mat_projection[16],
      const float mat_modelview[16],
      const std::vector<float> &vtx_val,
      double val_min,
      double val_max) const {
    namespace dfm2 = ::delfem2;
    const unsigned int nvtx = vtx_pos.size();
    assert(vtx_val.size() == nvtx);
    for (unsigned int ivtx = 0; ivtx < nvtx; ++ivtx) {
      float p[3] = {static_cast<float>(vtx_pos[ivtx][0]), static_cast<float>(vtx_pos[ivtx][1]), 0.f};
      if (ndim == 3) { p[2] = vtx_pos[ivtx][2]; }
      const dfm2::CMat4f mmv = dfm2::CMat4f(mat_modelview)
          * dfm2::CMat4f::Translation(p)
          * dfm2::CMat4f::AffineScale(radius_sphere);
      double v0 = (vtx_val[ivtx] - val_min)/(val_max-val_min);
      const double umin = 1-v0;
      const double umax = 2-v0;
      sphere.Draw(mat_projection, mmv.data(),umin,umax);
    }
    if (nvtx < 2) { return; }
    for (unsigned int ivtx = 0; ivtx < nvtx - 1; ++ivtx) {
      unsigned int jvtx = ivtx + 1;
      dfm2::CVec3f p0_(vtx_pos[ivtx][0], vtx_pos[ivtx][1], 0.f);
      dfm2::CVec3f p1_(vtx_pos[jvtx][0], vtx_pos[jvtx][1], 0.f);
      if (ndim == 3) {
        p0_.z = vtx_pos[ivtx][2];
        p1_.z = vtx_pos[jvtx][2];
      }
      dfm2::CMat3f R3 = Mat3_MinimumRotation(
          dfm2::CVec3f{0, 1, 0},
          p1_ - p0_);
      dfm2::CMat4f R4 = dfm2::CMat4f::Mat3(R3.data());
      const dfm2::CMat4f mmv = dfm2::CMat4f(mat_modelview)
          * dfm2::CMat4f::Translation(((p0_ + p1_) * 0.5).data())
          * R4
          * dfm2::CMat4f::ScaleXYZ(radius_cylinder, (p1_ - p0_).norm(), radius_cylinder);
      const double va = (vtx_val[ivtx]-val_min)/(val_max-val_min);
      const double vb = (vtx_val[jvtx]-val_min)/(val_max-val_min);
      const double umin = va/(va-vb);
      const double umax = (va-1)/(va-vb);
      cylinder.Draw(mat_projection, mmv.data(),umin,umax);
    }
  }

 public:
  Drawer_MeshColormap sphere;
  Drawer_MeshColormap cylinder;
  float radius_sphere = 1.f;
  float radius_cylinder = 1.f;
};

}

#endif //DRAWER_POLYLINECOLORMAP_H_
