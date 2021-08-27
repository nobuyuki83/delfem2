#ifndef DFM2_OPENGLEIGEN_FUNCS_H
#define DFM2_OPENGLEIGEN_FUNCS_H

#include <Eigen/Core>

#include "delfem2/dfm2_inline.h"
#include "delfem2/opengl/old/mshuni.h"

namespace delfem2 {
namespace opengleigen {

template<class VEC>
void myGlVertex3(const VEC &v) {
  ::glVertex3f(v(0), v(1), v(2));
}

template<>
void myGlVertex3(const Eigen::Vector3d &v) {
  ::glVertex3d(v(0), v(1), v(2));
}

template<typename VEC>
void DrawPoints(std::vector<VEC, Eigen::aligned_allocator<VEC> > &aXYZ) {
  ::glBegin(GL_POINTS);
  for (const auto &xyz : aXYZ) {
    myGlVertex3(xyz);
  }
  ::glEnd();
}

template<typename VEC>
void DrawMeshTri3(
    std::vector<unsigned int> &aTri,
    std::vector<VEC, Eigen::aligned_allocator<VEC> > &aXYZ) {
  const unsigned int nTri = aTri.size() / 3;
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < nTri; itri++) {
    const unsigned int i1 = aTri[itri * 3 + 0];
    const unsigned int i2 = aTri[itri * 3 + 1];
    const unsigned int i3 = aTri[itri * 3 + 2];
    myGlVertex3(aXYZ[i1]);
    myGlVertex3(aXYZ[i2]);
    myGlVertex3(aXYZ[i3]);
  }
  ::glEnd();
}

template<typename VEC>
void DrawMeshTri3_Edge(
    std::vector<unsigned int> &aTri,
    std::vector<VEC, Eigen::aligned_allocator<VEC> > &aXYZ) {
  const unsigned int nTri = aTri.size() / 3;
  ::glBegin(GL_LINES);
  for (unsigned int itri = 0; itri < nTri; itri++) {
    const unsigned int i1 = aTri[itri * 3 + 0];
    const unsigned int i2 = aTri[itri * 3 + 1];
    const unsigned int i3 = aTri[itri * 3 + 2];
    myGlVertex3(aXYZ[i1]);
    myGlVertex3(aXYZ[i2]);
    myGlVertex3(aXYZ[i2]);
    myGlVertex3(aXYZ[i3]);
    myGlVertex3(aXYZ[i3]);
    myGlVertex3(aXYZ[i1]);
  }
  ::glEnd();
}

void DrawMeshTri3_Edge_EigenMats(
    const Eigen::Matrix<double, -1, 3, Eigen::RowMajor> &V,
    const Eigen::Matrix<unsigned int, -1, 3, Eigen::RowMajor> &F) {
  ::delfem2::opengl::DrawMeshTri3D_Edge(
      V.data(), V.rows(),
      F.data(), F.rows());
}

} // opengleigen
} // delfem2


#endif // #define DFM2_OPENGLEIGEN_FUNCS_H