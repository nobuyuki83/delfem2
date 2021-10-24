//
// Created by Nobuyuki Umetani on 2021-10-08.
//

#ifndef DRAWER_SPHERE_H_
#define DRAWER_SPHERE_H_

#include "delfem2/mshprimitive.h"
#include "delfem2/mat4.h"
#include "delfem2/opengl/new/shdr_msh.h"

#include "delfem2/vec3.h"
#include "delfem2/geo3_v23m34q.h"

namespace delfem2::opengl {

class Drawer_Sphere : public CShader_Mesh {
 public:
  void InitGL() {
    CShader_Mesh::Compile();
    std::vector<float> vtx_xyz;
    std::vector<unsigned int> tri_vtx;
    delfem2::MeshTri3D_Sphere(
        vtx_xyz, tri_vtx,
        1.f, 32, 32);
    CShader_Mesh::Initialize(vtx_xyz, 3, tri_vtx, GL_TRIANGLES);
  }
  void Draw(
      float radius,
      const float trans[3],
      const float mat_projection[16],
      const float mat_modelview[16]) const {
    namespace dfm2 = ::delfem2;
    const dfm2::CMat4f mmv = dfm2::CMat4f(mat_modelview)
        * dfm2::CMat4f::Translation({trans[0],trans[1],trans[2]})
        * dfm2::CMat4f::AffineScale(radius);
    CShader_Mesh::Draw(mat_projection, mmv.data());
  }
};

}


#endif //DRAWER_SPHERE_H_
