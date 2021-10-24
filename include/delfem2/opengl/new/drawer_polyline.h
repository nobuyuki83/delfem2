//
// Created by Nobuyuki Umetani on 2021-10-08.
//

#ifndef DRAWER_POLYLINE_H_
#define DRAWER_POLYLINE_H_

#include "delfem2/opengl/new/drawer_sphere.h"
#include "delfem2/vec3.h"
#include "delfem2/geo3_v23m34q.h"

namespace delfem2::opengl {

class Drawer_Polyline {
 public:
  void InitGL() {
    {
      sphere.Compile();
      std::vector<float> vtx_xyz;
      std::vector<unsigned int> tri_vtx;
      delfem2::MeshTri3D_Sphere(
          vtx_xyz, tri_vtx,
          1.f, 32, 32);
      sphere.Initialize(vtx_xyz, 3, tri_vtx, GL_TRIANGLES);
    }
    {
      cylinder.Compile();
      std::vector<float> vtx_xyz;
      std::vector<unsigned int> tri_vtx;
      delfem2::MeshTri3D_CylinderOpen(
          vtx_xyz, tri_vtx,
          1.f, 1.f, 32, 1);
      cylinder.Initialize(vtx_xyz, 3, tri_vtx, GL_TRIANGLES);
    }
  }
  template <typename T0>
  void Draw3(
      const std::vector<T0>& vtx_xyz,
      float radius,
      const float mat_projection[16],
      const float mat_modelview[16]) const {
    namespace dfm2 = ::delfem2;
    unsigned int nvtx = vtx_xyz.size()/3;
    for(unsigned int ivtx=0;ivtx<nvtx;++ivtx) {
      const dfm2::CMat4f mmv = dfm2::CMat4f(mat_modelview)
          * dfm2::CMat4f::Translation({vtx_xyz[ivtx*3+0], vtx_xyz[ivtx*3+1], vtx_xyz[ivtx*3+2]})
          * dfm2::CMat4f::AffineScale(radius);
      sphere.Draw(mat_projection, mmv.data());
    }
    if( nvtx < 2 ){ return; }
    for(unsigned int ivtx=0;ivtx<nvtx-1;++ivtx) {
      unsigned int jvtx = ivtx+1;
      dfm2::CVec3f p0(vtx_xyz.data()+ivtx*3);
      dfm2::CVec3f p1(vtx_xyz.data()+jvtx*3);
      dfm2::CMat3f R3 = Mat3_MinimumRotation(
          dfm2::CVec3f{0,1,0},
          p1-p0);
      dfm2::CMat4f R4 = dfm2::CMat4f::Mat3(R3.data());
      const dfm2::CMat4f mmv = dfm2::CMat4f(mat_modelview)
          * dfm2::CMat4f::Translation(((p0+p1)*0.5).data())
          * R4
          * dfm2::CMat4f::ScaleXYZ(radius,(p1-p0).norm(), radius);
      cylinder.Draw(mat_projection, mmv.data());
    }
  }
 public:
  std::array<float,4> color{1.f, 1.f, 1.f, 1.f};
  CShader_Mesh sphere;
  CShader_Mesh cylinder;
};

}


#endif //DRAWER_SPHERE_H_
