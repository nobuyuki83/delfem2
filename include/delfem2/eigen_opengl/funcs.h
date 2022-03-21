#ifndef DFM2_EIGEN_OPENGL_FUNCS_H
#define DFM2_EIGEN_OPENGL_FUNCS_H

#include <Eigen/Core>

#include "delfem2/dfm2_inline.h"
#include "delfem2/opengl/old/mshuni.h"

namespace delfem2::eigen_opengl {

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

template <typename REAL>
void DrawMeshTri3_Edge(
    const Eigen::Matrix<REAL, -1, 3, Eigen::RowMajor>& V0,
    const Eigen::Matrix<unsigned int, -1, 3, Eigen::RowMajor>& F0) {
  ::glBegin(GL_LINES);
  for(unsigned int ifc=0;ifc<F0.rows();++ifc){
    unsigned int i0 = F0(ifc,0);
    unsigned int i1 = F0(ifc,1);
    unsigned int i2 = F0(ifc,2);
    ::glVertex3dv(V0.row(i0).data());
    ::glVertex3dv(V0.row(i1).data());
    ::glVertex3dv(V0.row(i1).data());
    ::glVertex3dv(V0.row(i2).data());
    ::glVertex3dv(V0.row(i2).data());
    ::glVertex3dv(V0.row(i0).data());
  }
  ::glEnd();
}

template <typename REAL>
void DrawMeshTri3_FaceFlatNorm(
    const Eigen::Matrix<REAL, -1, 3, Eigen::RowMajor>& V0,
    const Eigen::Matrix<unsigned int, -1, 3, Eigen::RowMajor>& F0) {
  ::glBegin(GL_TRIANGLES);
  for(unsigned int ifc=0;ifc<F0.rows();++ifc){
    unsigned int i0 = F0(ifc,0);
    unsigned int i1 = F0(ifc,1);
    unsigned int i2 = F0(ifc,2);
    const Eigen::Vector3<REAL> v0 = V0.row(i0);
    const Eigen::Vector3<REAL> v1 = V0.row(i1);
    const Eigen::Vector3<REAL> v2 = V0.row(i2);
    const Eigen::Vector3<REAL> n = (v1-v0).cross(v2-v0).normalized();
    ::glNormal3dv(n.data());
    ::glVertex3dv(v0.data());
    ::glVertex3dv(v1.data());
    ::glVertex3dv(v2.data());

  }
  ::glEnd();
}

} // delfem2::eigen_opengl

#endif // #define DFM2_EIGEN_OPENGL_FUNCS_H
