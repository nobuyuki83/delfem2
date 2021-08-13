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

#ifndef DFM2_OPENGL_OLD_FUNCS_H
#define DFM2_OPENGL_OLD_FUNCS_H

#include <string>
#include <vector>

#include "delfem2/dfm2_inline.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846264338327950288
#endif

namespace delfem2 {
namespace opengl {

// x = ax*[x] + bx
// y = ay*[y] + by
DFM2_INLINE void DrawCharacter(
    int *pChr,
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


//DFM2_INLINE void DrawRectangle_FullCanvas();
DFM2_INLINE void showdepth();

// --------------------------------------------------
DFM2_INLINE void getPosOnScreen_Camera2D(
    double &x, double &y,
    int i, int j);

DFM2_INLINE void setGL_Camera2D();


// ===============================================
// draw primitives

DFM2_INLINE void DrawAxis(double s);

DFM2_INLINE void DrawSphere(int nla, int nlo);
DFM2_INLINE void DrawSphereAt(
    int nla,
    int nlo,
    double rad,
    double x,
    double y,
    double z);

DFM2_INLINE void DrawSphere_Edge(
    double radius_);

DFM2_INLINE void DrawTorus_Edge(
    double radius_,
    double radius_tube_);

/**
 * @param rad_longtitude longer radius of doughnuts
 * @param rad_meridian shoter radius of  doughnuts
 */
DFM2_INLINE void DrawTorus_Solid(
    double rad_longtitude,
    double rad_meridian,
    double scale_tex);

DFM2_INLINE void DrawCylinder_Face(
    const double *dir_, double radius_, const double *cent_);

DFM2_INLINE void DrawCylinder_Edge(
    const double *dir_, double radius_, const double *cent_);

DFM2_INLINE void DrawPlane_Edge(
    const double *origin_, const double *normal_);

// ========================================
// Draw Axis-Aligned Box

DFM2_INLINE void DrawBox_MinMaxXYZ(
    double x_min, double x_max,
    double y_min, double y_max,
    double z_min, double z_max);

DFM2_INLINE void DrawBox_MinMaxXYZ(
    double aabbMinMaxXYZ[6]);

DFM2_INLINE void DrawAABB3D_Edge(
    double cx, double cy, double cz,
    double wx, double wy, double wz);

DFM2_INLINE void DrawAABB3D_Edge(
    const double cw[6]);

/**
 * @brief draw bounding box with edge
 * @details this function do nothing when pmin[0]  > pmin[1] 
 */
DFM2_INLINE void DrawBox3_Edge(
    const double *pmin,
    const double *pmax);

DFM2_INLINE void DrawBox2_Edge(
    const double *pmin,
    const double *pmax);

DFM2_INLINE void DrawBox3_Face(
    const double *pmin,
    const double *pmax);


// ------------------------

class CAxisXYZ {
 public:
  CAxisXYZ() : len(1.0) {
    line_width = 1.0;
  }
  CAxisXYZ(double len) : len(len) {
    line_width = 1.0;
  }
  void Draw() const;
  std::vector<double> MinMaxXYZ() const {
    std::vector<double> mm(6, 0);
    mm[1] = len;
    mm[3] = len;
    mm[5] = len;
    return mm;
  }
 public:
  double len;
  double line_width;
};

} // namespace opengl
} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/old/funcs.cpp"
#endif

#endif
