//
// Created by Nobuyuki Umetani on 2021-10-08.
//

#ifndef DRAWER_SPHERE_H_
#define DRAWER_SPHERE_H_

#include "delfem2/mshprimitive.h"
#include "delfem2/mat4.h"
#include "delfem2/opengl/new/shdr_msh.h"

namespace delfem2::opengl {

class Drawer_Sphere {
 public:
  void InitGL() {
    drawer_mesh.Compile();
    std::vector<float> vtx_xyz;
    std::vector<unsigned int> tri_vtx;
    delfem2::MeshTri3D_Sphere(
        vtx_xyz, tri_vtx,
        1.f, 32, 32);
    drawer_mesh.Initialize(vtx_xyz, 3, tri_vtx, GL_TRIANGLES);
  }
  void Draw(
      float radius,
      const float trans[3],
      const float mat_projection[16],
      const float mat_modelview[16]) const {
    namespace dfm2 = ::delfem2;
    dfm2::CMat4f mmv = dfm2::CMat4f(mat_modelview) * dfm2::CMat4f::Translate(trans) * dfm2::CMat4f::Scale(radius);
    drawer_mesh.Draw(mat_projection, mmv.data());
  }
 public:
  delfem2::opengl::CShader_Mesh drawer_mesh;
};

}


#endif //DRAWER_SPHERE_H_
