/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @details the file depends on OpenGL 2.x.
 * There are a lot of legacy commands such as glBegin(),glEnd()
 */

#ifndef DFM2_FUNCS_GLOLD_H
#define DFM2_FUNCS_GLOLD_H

#include "delfem2/dfm2_inline.h"
#include <string>
#include <vector>

namespace delfem2{
namespace opengl{

// x = ax*[x] + bx
// y = ay*[y] + by
DFM2_INLINE void DrawCharacter
 (int* pChr,
  double ax, double bx,
  double ay, double by);

// x = ax*[x] + bx
// y = ay*[y] + by
DFM2_INLINE void DrawCharacter
 (char ic,
  double ax, double bx,
  double ay, double by);


DFM2_INLINE void setSomeLighting();
DFM2_INLINE void setSomeLighting2();
DFM2_INLINE void setSomeLighting3();
DFM2_INLINE void drawFloorShadow
 (void (*DrawObject)(), float yfloor, float wfloor);


DFM2_INLINE void DrawRectangle_FullCanvas();
DFM2_INLINE void showdepth();

// --------------------------------------------------
DFM2_INLINE void getPosOnScreen_Camera2D
 (double& x, double& y,
  int i, int j);

DFM2_INLINE void setGL_Camera2D();


// -------------------
// draw functions

DFM2_INLINE void DrawAxis(double s);

DFM2_INLINE void DrawSphere(int nla, int nlo);
DFM2_INLINE void DrawSphereAt
 (int nla, int nlo, double rad, double x, double y, double z);
DFM2_INLINE void DrawSphere_Edge(double radius_);

DFM2_INLINE void DrawTorus_Edge
 (double radius_, double radius_tube_);

DFM2_INLINE void DrawCylinder_Face
 (const double* dir_, double radius_, const double* cent_);

DFM2_INLINE void DrawCylinder_Edge
 (const double* dir_, double radius_, const double* cent_);

DFM2_INLINE void DrawPlane_Edge
 (const double* origin_, const double* normal_);

DFM2_INLINE void DrawBox_MinMaxXYZ
 (double x_min, double x_max,
  double y_min, double y_max,
  double z_min, double z_max);

DFM2_INLINE void DrawBox_MinMaxXYZ(double aabbMinMaxXYZ[6]);
DFM2_INLINE void DrawAABB3D_Edge
 (double cx, double cy, double cz,
  double wx, double wy, double wz);

DFM2_INLINE void DrawAABB3D_Edge
 (const double cw[6]);

DFM2_INLINE void Draw_AABB3D_MinMaxXYZ_Edge
 (double x_min, double x_max,
  double y_min, double y_max,
  double z_min, double z_max);

// -------------
// Draw Point

DFM2_INLINE void DrawPoints2D_Vectors(
    const double* aXY,
    unsigned int nXY,
    const double* aVal,
    int nstride,
    int noffset,
    double mag);

DFM2_INLINE void DrawPoints2d_Points(const std::vector<double>& aXY);

DFM2_INLINE void DrawPoints3d_Points(const std::vector<double>& aXYZ);

DFM2_INLINE void DrawPoints3d_NormVtx
 (const std::vector<double>& aXYZ,
  const std::vector<double>& aNrm,
  double scale);


// --------------
// Draw Line
DFM2_INLINE void DrawMeshLine3D_Edge
 (const double* aXYZ,
  unsigned int nXYZ,
  const unsigned int* aLine,
  unsigned int nLine);

// ------------
// Draw Tri

DFM2_INLINE void DrawMeshTri2D_Face
 (const std::vector<unsigned int>& aTri,
  const std::vector<double>& aXY);

DFM2_INLINE void DrawMeshTri2D_FaceDisp2D(
    const double* aXY,
    unsigned int nXY,
    const unsigned int* aTri,
    unsigned int nTri,
    const double* aDisp,
    int nstride);

DFM2_INLINE void DrawMeshTri2D_Edge
 (const double* aXY, unsigned int nXY,
  const unsigned int* aTri, unsigned int nTri);

DFM2_INLINE void DrawMeshTri2D_Edge
 (const std::vector<unsigned int>& aTri,
  const std::vector<double>& aXY);

//
DFM2_INLINE void DrawMeshTri3D_FaceEdge
 (const std::vector<double>& aXYZ,
  const std::vector<unsigned int>& aTri);

DFM2_INLINE void DrawMeshTri3D_FaceNorm
 (const double* paXYZ,
  const unsigned int* paTri, unsigned int nTri);

DFM2_INLINE void DrawMeshTri3D_FaceNorm
 (const std::vector<double>& aXYZ,
  const std::vector<unsigned int>& aTri);

DFM2_INLINE void DrawMeshTri3D_FaceNorm
 (const std::vector<double>& aXYZ,
  const std::vector<unsigned int>& aTri,
  const std::vector<double>& aNorm);

DFM2_INLINE void DrawMeshTri3D_FaceNorm
 (const std::vector<double>& aXYZ,
  const std::vector<unsigned int>& aTriVtx,
  const std::vector<double>& aNorm,
  const std::vector<unsigned int>& aTriNrm);

//
DFM2_INLINE void DrawMeshTri3DPart_FaceNorm
 (const std::vector<double>& aXYZ,
  const std::vector<int>& aTri,
  const std::vector<int>& aIndTri);

DFM2_INLINE void DrawMeshTri3D_FaceNorm_Flg
 (const std::vector<double>& aXYZ,
  const std::vector<int>& aTri,
  int iflg,
  const std::vector<int>& aFlgTri);

DFM2_INLINE void DrawMeshTri3D_FaceNorm_XYsym
 (const std::vector<double>& aXYZ,
  const std::vector<unsigned int>& aTri);

DFM2_INLINE void DrawMeshTri3D_FaceNormEdge(const std::vector<double>& aXYZ,
                            const std::vector<unsigned int>& aTri);
DFM2_INLINE void DrawMeshTri3D_FaceNorm_TexFace(const std::vector<double>& aXYZ,
                                    const std::vector<unsigned int>& aTri,
                                    const std::vector<double>& aTex);
DFM2_INLINE void DrawMeshTri3D_FaceNorm_TexVtx(const std::vector<double>& aXYZ,
                                   const std::vector<unsigned int>& aTri,
                                   const std::vector<double>& aTex);
// Edge 3D
DFM2_INLINE void DrawMeshTri3D_Edge(const double* aXYZ, unsigned int nXYZ,
                        const unsigned int* aTri, unsigned int nTri);
DFM2_INLINE void DrawMeshTri3D_Edge(const std::vector<double>& aXYZ,
                        const std::vector<unsigned int>& aTri);
DFM2_INLINE void DrawMeshTriMap3D_Edge(const std::vector<double>& aXYZ,
                           const std::vector<unsigned int>& aTri,
                           const std::vector<int>& map);

// -----------------
// Draw Quad

DFM2_INLINE void DrawMeshQuad3D_Edge
 (const double* aXYZ, unsigned int nXYZ,
  const unsigned int* aQuad, unsigned int nQuad);

DFM2_INLINE void DrawMeshQuad3D_Edge
 (const std::vector<double>& aXYZ,
  const std::vector<unsigned int>& aQuad);

DFM2_INLINE void DrawMeshQuad2D_Edge
 (const double* aXY, unsigned int nXY,
  const unsigned int* aQuad, unsigned int nQuad);

DFM2_INLINE void DrawMeshQuad2D_Edge
 (const std::vector<double>& aXY,
  const std::vector<unsigned int>& aQuad);

DFM2_INLINE void DrawMeshQuad3D_FaceNorm
 (const double* aXYZ,
  const unsigned int* aQuad, unsigned int nQuad);

DFM2_INLINE void DrawMeshQuad3D_FaceNorm
 (const std::vector<double>& aXYZ,
  const std::vector<unsigned int>& aQuad);

// -------------------
// Draw Tet

DFM2_INLINE void DrawMeshTet3D_Edge(
    const double* aXYZ,
    unsigned int nXYZ,
    const unsigned int* aTet,
    unsigned int nTet);

DFM2_INLINE void DrawMeshTet3D_EdgeDisp(const double* aXYZ,
                            const unsigned int* aTet, unsigned int nTet,
                            const double* aDisp,
                            double s0);
DFM2_INLINE void DrawMeshTet3D_FaceNorm(const double* aXYZ,
                            const unsigned int* aTet, unsigned int nTet);
DFM2_INLINE void DrawMeshTet3DSurface_FaceNorm(const std::vector<double>& aXYZ,
                               const std::vector<unsigned int>& aTet,
                               const std::vector<unsigned int>& aTetFace);
DFM2_INLINE void DrawMeshTet3DSurface_Edge(const std::vector<double>& aXYZ,
                           const std::vector<unsigned int>& aTet,
                           const std::vector<unsigned int>& aTetFace);
DFM2_INLINE void DrawMeshTet3D_FaceNormDisp(const double* aXYZ, int nXYZ,
                                const unsigned int* aTet, int nTet,
                                const double* aDisp);

// -------------
// Draw Hex

DFM2_INLINE void DrawMeshHex3D_Edge(const double* aXYZ, int nXYZ,
                        const unsigned int* aHex, int nHex);
DFM2_INLINE void DrawMeshHex3D_FaceNorm(const double* aXYZ,
                            const unsigned int* aHex, unsigned int nHex);
DFM2_INLINE void DrawHex3D_FaceNormDisp(const std::vector<double>& aXYZ,
                           const std::vector<int>& aHex,
                           const std::vector<double>& aDisp);
DFM2_INLINE void Draw_HexMeshFaceDisp(const std::vector<double>& aXYZ,
                          const std::vector<unsigned int>& aHex,
                          const std::vector<double>& aDisp);


// -----------
// Draw Mix

DFM2_INLINE void DrawMeshElem3D_FaceNorm(const std::vector<double>& aXYZ,
                             const std::vector<unsigned int>& aElemInd,
                             const std::vector<unsigned int>& aElem);

DFM2_INLINE void DrawMeshElem3D_FaceNorm(const std::vector<double>& aXYZ,
                             const std::vector<unsigned int>& aElemInd,
                             const std::vector<unsigned int>& aElem,
                             const std::vector<double>& aUV);

DFM2_INLINE void DrawMeshElemPart3D_FaceNorm_TexPoEl(const std::vector<double>& aXYZ,
                                         const std::vector<unsigned int>& aElemInd,
                                         const std::vector<unsigned int>& aElem,
                                         const std::vector<int>& aIndElem,
                                         const std::vector<double>& aUV);


// ------------------------

class CAxisXYZ {
public:
  CAxisXYZ(): len(1.0){
    line_width = 1.0;
  }
  CAxisXYZ(double len): len(len){
    line_width=1.0;
  }
  void Draw() const;
  std::vector<double> MinMaxXYZ() const{
    std::vector<double> mm(6,0);
    mm[1] = len;  mm[3] = len;  mm[5] = len;
    return mm;
  }
public:
  double len;
  double line_width;
};

} // namespace opengl
} // namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/opengl/funcs_glold.cpp"
#endif

#endif /* utility_gl_h */
