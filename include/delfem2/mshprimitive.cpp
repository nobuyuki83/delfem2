/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/mshprimitive.h"

#include <vector>
#include <cmath>
#include <cassert>

DFM2_INLINE void delfem2::MeshQuad2D_Grid(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aQuad,
    unsigned int nx,
    unsigned int ny) {
  unsigned int np = (nx + 1) * (ny + 1);
  aXYZ.resize(np * 2);
  for (unsigned int iy = 0; iy < ny + 1; ++iy) {
    for (unsigned int ix = 0; ix < nx + 1; ++ix) {
      int ip = iy * (nx + 1) + ix;
      aXYZ[ip * 2 + 0] = ix;
      aXYZ[ip * 2 + 1] = iy;
    }
  }
  aQuad.resize(nx * ny * 4);
  for (unsigned int iy = 0; iy < ny; ++iy) {
    for (unsigned int ix = 0; ix < nx; ++ix) {
      unsigned int iq = iy * nx + ix;
      aQuad[iq * 4 + 0] = (iy + 0) * (nx + 1) + (ix + 0);
      aQuad[iq * 4 + 1] = (iy + 0) * (nx + 1) + (ix + 1);
      aQuad[iq * 4 + 2] = (iy + 1) * (nx + 1) + (ix + 1);
      aQuad[iq * 4 + 3] = (iy + 1) * (nx + 1) + (ix + 0);
    }
  }
}

DFM2_INLINE void delfem2::MeshTri3D_Disk(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    double r,
    int nr,
    int nth) {
  aXYZ.clear();
  aTri.clear();
  const double pi = 3.1415926535;
  { // make coordinates
    const int npo = 1 + nr * nth;
    double dr = r / nr;
    double dth = 2.0 * pi / nth;
    aXYZ.reserve(npo * 3);
    aXYZ.push_back(0.0);
    aXYZ.push_back(0.0);
    aXYZ.push_back(0.0);
    for (int ir = 1; ir <= nr; ir++) {
      double ri = dr * ir;
      for (int ith = 0; ith < nth; ith++) {
        aXYZ.push_back(ri * cos(ith * dth));
        aXYZ.push_back(0);
        aXYZ.push_back(ri * sin(ith * dth));
      }
    }
  }
  int ntri = nth * (nr - 1) * 2 + nth;
  aTri.reserve(ntri * 3);
  for (int ith = 0; ith < nth; ith++) {
    aTri.push_back(0);
    aTri.push_back((ith + 1) % nth + 1);
    aTri.push_back((ith + 0) % nth + 1);
  }
  for (int ir = 0; ir < nr - 1; ir++) {
    for (int ith = 0; ith < nth; ith++) {
      int i1 = (ir + 0) * nth + 1 + (ith + 0) % nth;
      int i2 = (ir + 0) * nth + 1 + (ith + 1) % nth;
      int i3 = (ir + 1) * nth + 1 + (ith + 1) % nth;
      int i4 = (ir + 1) * nth + 1 + (ith + 0) % nth;
      aTri.push_back(i3);
      aTri.push_back(i1);
      aTri.push_back(i2);
      aTri.push_back(i4);
      aTri.push_back(i1);
      aTri.push_back(i3);
    }
  }
}

DFM2_INLINE void delfem2::MeshTri3D_CylinderOpen(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    double r, double l,
    int nr,
    int nl) {
  aXYZ.clear();
  aTri.clear();
  const double pi = 3.1415926535;
  double dl = l / nl;
  double dr = 2.0 * pi / nr;
  const int npo = (nl + 1) * nr;
  aXYZ.reserve(npo * 3);
  for (int il = 0; il < nl + 1; il++) {
    double y0 = -0.5 * l + il * dl;
    for (int ir = 0; ir < nr; ir++) {
      double x0 = r * cos(dr * ir);
      double z0 = r * sin(dr * ir);
      aXYZ.push_back(x0);
      aXYZ.push_back(y0);
      aXYZ.push_back(z0);
    }
  }
  // -------
  const int ntri = nl * nr * 2;
  aTri.reserve(ntri * 3);
  for (int il = 0; il < nl; il++) {
    for (int ir = 0; ir < nr; ir++) {
      const int i1 = (il + 0) * nr + (ir + 0) % nr;
      const int i2 = (il + 0) * nr + (ir + 1) % nr;
      const int i3 = (il + 1) * nr + (ir + 1) % nr;
      const int i4 = (il + 1) * nr + (ir + 0) % nr;
      //      std::cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<" "<<npo<<std::endl;
      aTri.push_back(i3);
      aTri.push_back(i2);
      aTri.push_back(i1);
      aTri.push_back(i4);
      aTri.push_back(i3);
      aTri.push_back(i1);
    }
  }
}

template<typename T>
DFM2_INLINE void delfem2::MeshTri3D_CylinderClosed(
    std::vector<T> &aXYZ,
    std::vector<unsigned int> &aTri,
    T r,
    T l,
    unsigned int nr,
    unsigned int nl) {
  aXYZ.clear();
  aTri.clear();
  if (nl < 1 || nr <= 2) { return; }
  const T pi = static_cast<T>(3.1415926535);
  T dl = l / nl;
  T dr = 2 * pi / nr;
  aXYZ.reserve((nr * (nl + 1) + 2) * 3);
  {
    aXYZ.push_back(0);
    aXYZ.push_back(-l / 2);
    aXYZ.push_back(0);
  }
  for (unsigned int il = 0; il < nl + 1; il++) {
    T y0 = -l / 2 + dl * il;
    for (unsigned int ilo = 0; ilo < nr; ilo++) {
      T x0 = r * std::cos(dr * ilo);
      T z0 = r * std::sin(dr * ilo);
      aXYZ.push_back(x0);
      aXYZ.push_back(y0);
      aXYZ.push_back(z0);
    }
  }
  {
    aXYZ.push_back(0);
    aXYZ.push_back(+l / 2);
    aXYZ.push_back(0);
  }
  // ------------------------------------
  unsigned int nla = nl + 2;
  unsigned int ntri = nr * (nla - 1) * 2 + nr * 2;
  aTri.reserve(ntri * 3);
  for (unsigned int ilo = 0; ilo < nr; ilo++) {
    aTri.push_back(0);
    aTri.push_back((ilo + 0) % nr + 1);
    aTri.push_back((ilo + 1) % nr + 1);
  }
  for (unsigned int ila = 0; ila < nla - 2; ila++) {
    for (unsigned int ilo = 0; ilo < nr; ilo++) {
      const unsigned int i1 = (ila + 0) * nr + 1 + (ilo + 0) % nr;
      const unsigned int i2 = (ila + 0) * nr + 1 + (ilo + 1) % nr;
      const unsigned int i3 = (ila + 1) * nr + 1 + (ilo + 1) % nr;
      const unsigned int i4 = (ila + 1) * nr + 1 + (ilo + 0) % nr;
      aTri.push_back(i3);
      aTri.push_back(i2);
      aTri.push_back(i1);
      aTri.push_back(i4);
      aTri.push_back(i3);
      aTri.push_back(i1);
    }
  }
  for (unsigned int ilo = 0; ilo < nr; ilo++) {
    aTri.push_back(nr * (nla - 1) + 1);
    aTri.push_back((nla - 2) * nr + 1 + (ilo + 1) % nr);
    aTri.push_back((nla - 2) * nr + 1 + (ilo + 0) % nr);
  }
  /*
   for(int itri=0;itri<aTri.size()/3;itri++){
   for(int inotri=0;inotri<3;++inotri){
   const int i0 = aTri[itri*3+inotri];
   assert( i0 >=0 && i0 < aXYZ.size()/3 );
   }
   }
   */
}
#ifdef DFM2_STATIC_LIBRARY
template DFM2_INLINE void delfem2::MeshTri3D_CylinderClosed(
    std::vector<float> &aXYZ,
    std::vector<unsigned int> &aTri,
    float r, float l,
    unsigned int nr, unsigned int nl);
template DFM2_INLINE void delfem2::MeshTri3D_CylinderClosed(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    double r, double l,
    unsigned int nr, unsigned int nl);
#endif

// ------------------------

DFM2_INLINE void delfem2::MeshTri3D_Sphere(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    double radius,
    int nlong,
    int nlat) {
  aXYZ.clear();
  aTri.clear();
  if (nlong <= 1 || nlat <= 2) { return; }
  const double pi = 3.1415926535;
  double dl = pi / nlong;
  double dr = 2.0 * pi / nlat;
  aXYZ.reserve((nlat * (nlong - 1) + 2) * 3);
  for (int ila = 0; ila < nlong + 1; ila++) {
    double y0 = cos(dl * ila);
    double r0 = sin(dl * ila);
    for (int ilo = 0; ilo < nlat; ilo++) {
      double x0 = r0 * sin(dr * ilo);
      double z0 = r0 * cos(dr * ilo);
      aXYZ.push_back(radius * x0);
      aXYZ.push_back(radius * y0);
      aXYZ.push_back(radius * z0);
      if (ila == 0 || ila == nlong) { break; }
    }
  }
  //
  int ntri = nlat * (nlong - 1) * 2 + nlat * 2;
  aTri.reserve(ntri * 3);
  for (int ilo = 0; ilo < nlat; ilo++) {
    aTri.push_back(0);
    aTri.push_back((ilo + 0) % nlat + 1);
    aTri.push_back((ilo + 1) % nlat + 1);
  }
  for (int ila = 0; ila < nlong - 2; ila++) {
    for (int ilo = 0; ilo < nlat; ilo++) {
      int i1 = (ila + 0) * nlat + 1 + (ilo + 0) % nlat;
      int i2 = (ila + 0) * nlat + 1 + (ilo + 1) % nlat;
      int i3 = (ila + 1) * nlat + 1 + (ilo + 1) % nlat;
      int i4 = (ila + 1) * nlat + 1 + (ilo + 0) % nlat;
      aTri.push_back(i3);
      aTri.push_back(i2);
      aTri.push_back(i1);
      aTri.push_back(i4);
      aTri.push_back(i3);
      aTri.push_back(i1);
    }
  }
  for (int ilo = 0; ilo < nlat; ilo++) {
    aTri.push_back(nlat * (nlong - 1) + 1);
    aTri.push_back((nlong - 2) * nlat + 1 + (ilo + 1) % nlat);
    aTri.push_back((nlong - 2) * nlat + 1 + (ilo + 0) % nlat);
  }
}

template<typename T>
DFM2_INLINE void delfem2::MeshTri3D_Cube(
    std::vector<T> &aXYZ,
    std::vector<unsigned int> &aTri,
    unsigned int n) {
  aXYZ.clear();
  aTri.clear();
  if (n < 1) { return; }
  const T r = static_cast<T>(1.0) / n;
  const unsigned int np = 4 * n * (n + 1) + (n - 1) * (n - 1) * 2;
  aXYZ.reserve(np * 3);
  constexpr T half(0.5);
  for (unsigned int iz = 0; iz < n + 1; ++iz) {  // height
    for (unsigned int ix = 0; ix < n; ++ix) {
      aXYZ.push_back(-half + r * ix);
      aXYZ.push_back(-half);
      aXYZ.push_back(-half + r * iz);
    }
    for (unsigned int iy = 0; iy < n; ++iy) {
      aXYZ.push_back(+half);
      aXYZ.push_back(-half + r * iy);
      aXYZ.push_back(-half + r * iz);
    }
    for (int ix = n; ix > 0; --ix) {
      aXYZ.push_back(-half + r * ix);
      aXYZ.push_back(+half);
      aXYZ.push_back(-half + r * iz);
    }
    for (int iy = n; iy > 0; --iy) {
      aXYZ.push_back(-half);
      aXYZ.push_back(-half + r * iy);
      aXYZ.push_back(-half + r * iz);
    }
  }
  for (unsigned int iy = 1; iy < n; ++iy) {
    for (unsigned int ix = 1; ix < n; ++ix) {
      aXYZ.push_back(-half + r * ix);
      aXYZ.push_back(-half + r * iy);
      aXYZ.push_back(-half);
    }
  }
  for (unsigned int iy = 1; iy < n; ++iy) {
    for (unsigned int ix = 1; ix < n; ++ix) {
      aXYZ.push_back(-half + r * ix);
      aXYZ.push_back(-half + r * iy);
      aXYZ.push_back(+half);
    }
  }
  // -------------------------------------------------
  const unsigned int ntri = n * n * 6 * 2;
  aTri.reserve(ntri * 3);
  for (unsigned int iz = 0; iz < n; ++iz) {
    for (unsigned int ixy = 0; ixy < 4 * n; ++ixy) {
      unsigned int i0 = ixy + 4 * n * iz;
      unsigned int i1 = (ixy + 1) % (4 * n) + 4 * n * iz;
      unsigned int i2 = (ixy + 1) % (4 * n) + 4 * n * (iz + 1);
      unsigned int i3 = ixy + 4 * n * (iz + 1);
      aTri.push_back(i0);
      aTri.push_back(i1);
      aTri.push_back(i2);
      // ----------------------------
      aTri.push_back(i2);
      aTri.push_back(i3);
      aTri.push_back(i0);
    }
  }
  // bottom
  for (unsigned int ix = 0; ix < n; ++ix) {
    for (unsigned int iy = 0; iy < n; ++iy) {
      unsigned int i0, i1, i2, i3;
      i0 = 4 * n * (n + 1) + (iy - 1) * (n - 1) + (ix - 1);
      i1 = 4 * n * (n + 1) + (iy - 1) * (n - 1) + (ix + 0);
      i2 = 4 * n * (n + 1) + (iy + 0) * (n - 1) + (ix + 0);
      i3 = 4 * n * (n + 1) + (iy + 0) * (n - 1) + (ix - 1);
      if (ix == 0) {
        i0 = (iy == 0) ? 0 : 4 * n - iy;
        i3 = 4 * n - iy - 1;
      }
      if (ix == n - 1) {
        i1 = n + iy;
        i2 = n + iy + 1;
      }
      if (iy == 0) {
        i0 = ix;
        i1 = ix + 1;
      }
      if (iy == n - 1) {
        i2 = 3 * n - ix - 1;
        i3 = 3 * n - ix + 0;
      }
      aTri.push_back(i1);
      aTri.push_back(i0);
      aTri.push_back(i2);
      //
      aTri.push_back(i3);
      aTri.push_back(i2);
      aTri.push_back(i0);
    }
  }
  // top
  unsigned int nps = 4 * n * (n + 1); // side vertex
  unsigned int nps0 = 4 * n * n; // side vertex
  for (unsigned int ix = 0; ix < n; ++ix) {
    for (unsigned int iy = 0; iy < n; ++iy) {
      unsigned int i0, i1, i2, i3;
      i0 = nps + (n - 1) * (n - 1) + (iy - 1) * (n - 1) + (ix - 1);
      i1 = nps + (n - 1) * (n - 1) + (iy - 1) * (n - 1) + (ix + 0);
      i2 = nps + (n - 1) * (n - 1) + (iy + 0) * (n - 1) + (ix + 0);
      i3 = nps + (n - 1) * (n - 1) + (iy + 0) * (n - 1) + (ix - 1);
      if (ix == 0) {
        i0 = (iy == 0) ? nps0 : nps0 + 4 * n - iy;
        i3 = nps0 + 4 * n - iy - 1;
      }
      if (ix == n - 1) {
        i1 = nps0 + n + iy;
        i2 = nps0 + n + iy + 1;
      }
      if (iy == 0) {
        i0 = nps0 + ix;
        i1 = nps0 + ix + 1;
      }
      if (iy == n - 1) {
        i2 = nps0 + 3 * n - ix - 1;
        i3 = nps0 + 3 * n - ix + 0;
      }
      aTri.push_back(i0);
      aTri.push_back(i1);
      aTri.push_back(i2);
      //
      aTri.push_back(i2);
      aTri.push_back(i3);
      aTri.push_back(i0);
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MeshTri3D_Cube(
    std::vector<float> &aXYZ,
    std::vector<unsigned int> &aTri,
    unsigned int n);
template void delfem2::MeshTri3D_Cube(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    unsigned int n);
#endif

// ---------------------------------------

// p0: -x, -y, -z
// p1: +x, -y, -z
// p2: -x, +y, -z
// p3: +x, +y, -z
// p4: -x, -y, +z
// p5: +x, -y, +z
// p6: -x, +y, +z
// p7: +x, +y, +z
// f0: -x
// f1: +x
// f2: -y
// f3: +y
// f4: -z
// f5: +z
std::vector<unsigned int>
DFM2_INLINE delfem2::MeshQuadTopo_CubeVox() {
  return {0, 4, 6, 2,
          1, 3, 7, 5,
          0, 1, 5, 4,
          2, 6, 7, 3,
          0, 2, 3, 1,
          4, 5, 7, 6};
}

template<typename REAL>
void delfem2::MeshQuad3_CubeVox(
    std::vector<REAL> &aXYZ,
    std::vector<unsigned int> &aQuad,
    const REAL bbmin[3],
    const REAL bbmax[3]) {
  aXYZ = std::vector<REAL>{
      bbmin[0], bbmin[1], bbmin[2],
      bbmax[0], bbmin[1], bbmin[2],
      bbmin[0], bbmax[1], bbmin[2],
      bbmax[0], bbmax[1], bbmin[2],
      bbmin[0], bbmin[1], bbmax[2],
      bbmax[0], bbmin[1], bbmax[2],
      bbmin[0], bbmax[1], bbmax[2],
      bbmax[0], bbmax[1], bbmax[2]};
  aQuad = MeshQuadTopo_CubeVox();
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MeshQuad3_CubeVox(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aQuad,
    const double bbmin[3], const double bbmax[3]);
template void delfem2::MeshQuad3_CubeVox(
    std::vector<float> &aXYZ,
    std::vector<unsigned int> &aQuad,
    const float bbmin[3], const float bbmax[3]);
#endif

// --------------------------------------

DFM2_INLINE void delfem2::MeshTri3D_Icosahedron(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri) {
  const double p = (1 + std::sqrt(5.)) * 0.5;
  aXYZ = std::vector<double>{
      +0, -1, -p,
      +0, -1, +p,
      +0, +1, -p,
      +0, +1, +p,
      -p, +0, -1,
      +p, +0, -1,
      -p, +0, +1,
      +p, +0, +1,
      -1, -p, +0,
      -1, +p, +0,
      +1, -p, +0,
      +1, +p, +0};
  //
  aTri = std::vector<unsigned int>{
      7, 11, 3,
      11, 9, 3,
      9, 6, 3,
      6, 1, 3,
      1, 7, 3,
      2, 5, 0,
      4, 2, 0,
      8, 4, 0,
      10, 8, 0,
      5, 10, 0,
      11, 7, 5,
      9, 11, 2,
      6, 9, 4,
      1, 6, 8,
      7, 1, 10,
      5, 2, 11,
      2, 4, 9,
      4, 8, 6,
      8, 10, 1,
      10, 5, 7};
}

// ------------------------------------------------------

template<typename T>
void delfem2::MeshTri3_Capsule(
    std::vector<T> &aXYZ,
    std::vector<unsigned int> &aTri,
    T r,
    T l,
    unsigned int nc,
    unsigned int nr,
    unsigned int nl) {
  constexpr T half(0.5);
  constexpr T pi = static_cast<T>(M_PI);
  MeshTri3D_CylinderClosed(aXYZ, aTri,
                           (T) 1., (T) 1.,
                           nc, 2 * nr + nl - 2);
  assert(aXYZ.size() / 3 == (2 * nr + nl - 1) * nc + 2);
  {
    aXYZ[0 * 3 + 0] = 0;
    aXYZ[0 * 3 + 1] = -l * half - r;
    aXYZ[0 * 3 + 2] = 0;
  }
  for (unsigned int ir = 0; ir < nr; ++ir) {
    const T t0 = pi * half * (nr - 1 - ir) / nr;
    const T y0 = -l * half - r * sin(t0);
    const T c0 = r * std::cos(t0);
    for (unsigned int ic = 0; ic < nc; ++ic) {
      aXYZ[(1 + ir * nc + ic) * 3 + 0] = c0 * std::cos(2 * pi * ic / nc);
      aXYZ[(1 + ir * nc + ic) * 3 + 1] = y0;
      aXYZ[(1 + ir * nc + ic) * 3 + 2] = c0 * std::sin(2 * pi * ic / nc);
    }
  }
  for (unsigned int il = 0; il < nl - 1; ++il) {
    const T y0 = -l * half + (il + 1) * l / nl;
    for (unsigned int ic = 0; ic < nc; ++ic) {
      aXYZ[(1 + (il + nr) * nc + ic) * 3 + 0] = r * cos(2 * pi * ic / nc);
      aXYZ[(1 + (il + nr) * nc + ic) * 3 + 1] = y0;
      aXYZ[(1 + (il + nr) * nc + ic) * 3 + 2] = r * sin(2 * pi * ic / nc);
    }
  }
  for (unsigned int ir = 0; ir < nr; ++ir) {
    const T t0 = pi * half * ir / nr;
    const T y0 = +l * half + r * std::sin(t0);
    const T c0 = r * std::cos(t0);
    for (unsigned int ic = 0; ic < nc; ++ic) {
      aXYZ[(1 + (ir + nl + nr - 1) * nc + ic) * 3 + 0] = c0 * std::cos(2 * pi * ic / nc);
      aXYZ[(1 + (ir + nl + nr - 1) * nc + ic) * 3 + 1] = y0;
      aXYZ[(1 + (ir + nl + nr - 1) * nc + ic) * 3 + 2] = c0 * std::sin(2 * pi * ic / nc);
    }
  }
  {
    const size_t np = aXYZ.size() / 3;
    aXYZ[(np - 1) * 3 + 0] = 0;
    aXYZ[(np - 1) * 3 + 1] = +l * half + r;
    aXYZ[(np - 1) * 3 + 2] = 0;
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MeshTri3_Capsule(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    double r, double l,
    unsigned int nc, unsigned int nr, unsigned int nl);
template void delfem2::MeshTri3_Capsule(
    std::vector<float> &aXYZ,
    std::vector<unsigned int> &aTri,
    float r, float l,
    unsigned int nc, unsigned int nr, unsigned int nl);
#endif

// ------------------------------------------------------------

template<typename T>
DFM2_INLINE void delfem2::MeshTri3_Torus(
    std::vector<T> &aXYZ,
    std::vector<unsigned int> &aTri,
    T radius_, // latitude
    T radius_tube_, // meridian
    unsigned int nlg, // latitude
    unsigned int nlt) // meridian
{
  const T rlg = static_cast<T>(M_PI * 2) / nlg;  // latitude
  const T rlt = static_cast<T>(M_PI * 2) / nlt;  // meridian
  aXYZ.resize(nlg * nlt * 3);
  for (unsigned int ilg = 0; ilg < nlg; ilg++) {
    for (unsigned int ilt = 0; ilt < nlt; ilt++) {
      aXYZ[(ilg * nlt + ilt) * 3 + 0] = (radius_ + radius_tube_ * std::cos(ilt * rlt)) * std::sin(ilg * rlg);
      aXYZ[(ilg * nlt + ilt) * 3 + 1] = (radius_ + radius_tube_ * std::cos(ilt * rlt)) * std::cos(ilg * rlg);
      aXYZ[(ilg * nlt + ilt) * 3 + 2] = radius_tube_ * sin(ilt * rlt);
    }
  }
  aTri.resize(nlg * nlt * 2 * 3);
  for (unsigned int ilg = 0; ilg < nlg; ilg++) {
    for (unsigned int ilt = 0; ilt < nlt; ilt++) {
      unsigned int iug = (ilg == nlg - 1) ? 0 : ilg + 1;
      unsigned int iut = (ilt == nlt - 1) ? 0 : ilt + 1;
      aTri[(ilg * nlt + ilt) * 6 + 0] = ilg * nlt + ilt;
      aTri[(ilg * nlt + ilt) * 6 + 2] = iug * nlt + ilt;
      aTri[(ilg * nlt + ilt) * 6 + 1] = iug * nlt + iut;
      //
      aTri[(ilg * nlt + ilt) * 6 + 3] = ilg * nlt + ilt;
      aTri[(ilg * nlt + ilt) * 6 + 5] = iug * nlt + iut;
      aTri[(ilg * nlt + ilt) * 6 + 4] = ilg * nlt + iut;
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MeshTri3_Torus(
    std::vector<float> &aXYZ,
    std::vector<unsigned int> &aTri,
    float radius_,      // latitude
    float radius_tube_,  // meridian
    unsigned int nr,     // latitude
    unsigned int nl);    // meridian

template void delfem2::MeshTri3_Torus(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    double radius_, double radius_tube_,
    unsigned int nr, unsigned int nl);
#endif

template<typename T>
DFM2_INLINE void delfem2::MeshHex3_Grid(
    std::vector<T> &aXYZ,
    std::vector<unsigned int> &aHex,
    unsigned int nx,
    unsigned int ny,
    unsigned int nz,
    T elen) {
  aXYZ.resize((nz + 1) * (ny + 1) * (nx + 1) * 3, 0);
  for (unsigned int iz = 0; iz < nz + 1; ++iz) {
    for (unsigned int iy = 0; iy < ny + 1; ++iy) {
      for (unsigned int ix = 0; ix < nx + 1; ++ix) {
        aXYZ[(iz * (ny + 1) * (nx + 1) + iy * (nx + 1) + ix) * 3 + 0] = ix * elen;
        aXYZ[(iz * (ny + 1) * (nx + 1) + iy * (nx + 1) + ix) * 3 + 1] = iy * elen;
        aXYZ[(iz * (ny + 1) * (nx + 1) + iy * (nx + 1) + ix) * 3 + 2] = iz * elen;
      }
    }
  }
  aHex.resize(nx * ny * nz * 8);
  for (unsigned int iz = 0; iz < nz; ++iz) {
    for (unsigned int iy = 0; iy < ny; ++iy) {
      for (unsigned int ix = 0; ix < nx; ++ix) {
        unsigned int ih0 = iz * ny * nx + iy * nx + ix;
        aHex[ih0 * 8 + 0] = (iz + 0) * (ny + 1) * (nx + 1) + (iy + 0) * (nx + 1) + (ix + 0);
        aHex[ih0 * 8 + 1] = (iz + 0) * (ny + 1) * (nx + 1) + (iy + 0) * (nx + 1) + (ix + 1);
        aHex[ih0 * 8 + 2] = (iz + 0) * (ny + 1) * (nx + 1) + (iy + 1) * (nx + 1) + (ix + 1);
        aHex[ih0 * 8 + 3] = (iz + 0) * (ny + 1) * (nx + 1) + (iy + 1) * (nx + 1) + (ix + 0);
        aHex[ih0 * 8 + 4] = (iz + 1) * (ny + 1) * (nx + 1) + (iy + 0) * (nx + 1) + (ix + 0);
        aHex[ih0 * 8 + 5] = (iz + 1) * (ny + 1) * (nx + 1) + (iy + 0) * (nx + 1) + (ix + 1);
        aHex[ih0 * 8 + 6] = (iz + 1) * (ny + 1) * (nx + 1) + (iy + 1) * (nx + 1) + (ix + 1);
        aHex[ih0 * 8 + 7] = (iz + 1) * (ny + 1) * (nx + 1) + (iy + 1) * (nx + 1) + (ix + 0);
      }
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MeshHex3_Grid(
    std::vector<float> &aXYZ,
    std::vector<unsigned int> &aHex,
    unsigned int nx,
    unsigned int ny,
    unsigned int nz,
    float elen);
template void delfem2::MeshHex3_Grid(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aHex,
    unsigned int nx,
    unsigned int ny,
    unsigned int nz,
    double elen);
#endif


// above: function
// =========================================================================
// below: class

template<typename REAL>
delfem2::CPlane<REAL>::CPlane(
    const double n[3], const double o[3]) {
  normal_[0] = n[0];
  normal_[1] = n[1];
  normal_[2] = n[2];
  //
  origin_[0] = o[0];
  origin_[1] = o[1];
  origin_[2] = o[2];
}

template<typename REAL>
double delfem2::CPlane<REAL>::Projection(
    double n[3],
    double px, double py, double pz) const // normal
{
  n[0] = normal_[0];
  n[1] = normal_[1];
  n[2] = normal_[2];
  return -(normal_[0] * (px - origin_[0]) + normal_[1] * (py - origin_[1]) + normal_[2] * (pz - origin_[2]));
}

// -------------------------------------------------------

template<typename REAL>
delfem2::CSphere<REAL>::CSphere(
    double r, const std::vector<double> &c, bool is_out) {
  cent_.resize(3);
  cent_[0] = c[0];
  cent_[1] = c[1];
  cent_[2] = c[2];
  radius_ = r;
  this->is_out_ = is_out;
}
template delfem2::CSphere<double>::CSphere
    (double r, const std::vector<double> &c, bool is_out);

// ----------------

// return penetration depth (inside is positive)
template<typename REAL>
double delfem2::CSphere<REAL>::Projection(
    double n[3],
    double px, double py, double pz) const // normal outward
{
  double dir[3] = {px - cent_[0], py - cent_[1], pz - cent_[2]};
  const double len = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  const double invlen = 1.0 / len;
  if (!is_out_) {
    n[0] = -dir[0] * invlen;
    n[1] = -dir[1] * invlen;
    n[2] = -dir[2] * invlen;
    return +len - radius_;
  }
  n[0] = dir[0] * invlen;
  n[1] = dir[1] * invlen;
  n[2] = dir[2] * invlen;
  return radius_ - len;
}
template double delfem2::CSphere<double>::Projection
    (double n[3],
     double px, double py, double pz) const; // normal outward

// -----------------------------------------

template<typename REAL>
unsigned int delfem2::CSphere<REAL>::FindInOut(
    double px, double py, double pz) const {
  double n[3];
  double pd = this->Projection(n, px, py, pz);
  if (!is_out_) pd *= -1.0;
  if (pd > 0) { return 0; }
  return 1;
}

template<typename REAL>
bool delfem2::CSphere<REAL>::IntersectionPoint(
    double p[3],
    const double o[3],
    const double d[3]) const {
  const double q[3] = {o[0] - cent_[0], o[1] - cent_[1], o[2] - cent_[2]};
  const double a = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  const double b = q[0] * d[0] + q[1] * d[1] + q[2] * d[2];
  const double c = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] - radius_ * radius_;
  const double det = b * b - a * c;
  if (det < 0) return false;
  const double t = (-b + sqrt(det)) / a;
  p[0] = o[0] + t * d[0];
  p[1] = o[1] + t * d[1];
  p[2] = o[2] + t * d[2];
  return true;
}


// --------------------------------------------------

template<typename REAL>
delfem2::CCylinder<REAL>::CCylinder(
    double r, const double cnt[3], const double dir[3], bool is_out) {
  cent_[0] = cnt[0];
  cent_[1] = cnt[1];
  cent_[2] = cnt[2];
  dir_[0] = dir[0];
  dir_[1] = dir[1];
  dir_[2] = dir[2];
  radius_ = r;
  this->is_out_ = is_out;
}

// return penetration depth (inside is positive)
template<typename REAL>
double delfem2::CCylinder<REAL>::Projection(
    double n[3],
    double px, double py, double pz) const // normal outward
{
  double dd = dir_[0] * dir_[0] + dir_[1] * dir_[1] + dir_[2] * dir_[2];
  double pod = (px - cent_[0]) * dir_[0] + (py - cent_[1]) * dir_[1] + (pz - cent_[2]) * dir_[2];
  double a = pod / dd;
  double hp[3] = {cent_[0] + a * dir_[0], cent_[1] + a * dir_[1], cent_[2] + a * dir_[2]};
  double d[3] = {px - hp[0], py - hp[1], pz - hp[2]};
  const double len = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  const double invlen = 1.0 / len;
  if (!is_out_) {
    n[0] = -d[0] * invlen;
    n[1] = -d[1] * invlen;
    n[2] = -d[2] * invlen;
    return +len - radius_;
  }
  n[0] = d[0] * invlen;
  n[1] = d[1] * invlen;
  n[2] = d[2] * invlen;
  return radius_ - len;
}

template<typename REAL>
unsigned int delfem2::CCylinder<REAL>::FindInOut(
    double px, double py, double pz) const {
  double n[3];
  double pd = this->Projection(n,
                               px, py, pz);
  if (!is_out_) pd *= -1.0;
  if (pd > 0) { return 0; }
  return 1;
}

template<typename REAL>
bool delfem2::CCylinder<REAL>::IntersectionPoint(
    double p[3],
    const double o[3], const double d[3]) const {
  const double q[3] = {o[0] - cent_[0], o[1] - cent_[1], o[2] - cent_[2]};
  const double a = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  const double b = q[0] * d[0] + q[1] * d[1] + q[2] * d[2];
  const double c = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] - radius_ * radius_;
  const double det = b * b - a * c;
  if (det < 0) return false;
  const double t = (-b + sqrt(det)) / a;
  p[0] = o[0] + t * d[0];
  p[1] = o[1] + t * d[1];
  p[2] = o[2] + t * d[2];
  return true;
}

// --------------------------------------------------------

// return penetration depth (inside is positive)
template<typename REAL>
double delfem2::CTorus<REAL>::Projection(
    double n[3],
    double px, double py, double pz) const // normal outward
{
  double dir[3] = {px - cent_[0], py - cent_[1], pz - cent_[2]};
  const double t = dir[2];
//	dir[0] -= t*norm_[0];
//	dir[1] -= t*norm_[1];
  dir[2] -= t;
  const double len = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  if (len < 1.0e-20) {    // vertical to the center of torus
    return radius_tube_ - radius_;
  }
  const double invlen = 1.0 / len;
  dir[0] *= invlen;
  dir[1] *= invlen;
  dir[2] *= invlen;
  double p[3];
  p[0] = cent_[0] + radius_ * dir[0];
  p[1] = cent_[1] + radius_ * dir[1];
  p[2] = cent_[2] + radius_ * dir[2];
  double dir2[3] = {px - p[0], py - p[1], pz - p[2]};
  const double len2 = sqrt(dir2[0] * dir2[0] + dir2[1] * dir2[1] + dir2[2] * dir2[2]);
  const double invlen2 = 1.0 / len2;
  n[0] = dir2[0] * invlen2;
  n[1] = dir2[1] * invlen2;
  n[2] = dir2[2] * invlen2;
  //		std::cout << len << " " << len2 << std::endl;
  return radius_tube_ - len2;
}
template double delfem2::CTorus<double>::Projection(
    double n[3],
    double px, double py, double pz) const; // normal outward

// ------------------

template<typename REAL>
unsigned int delfem2::CTorus<REAL>::FindInOut(
    double px, double py, double pz) const {
  double n[3];
  const double pd = this->Projection(n,
                                     px, py, pz);
  if (pd > 0) { return 0; }
  return 1;
}

