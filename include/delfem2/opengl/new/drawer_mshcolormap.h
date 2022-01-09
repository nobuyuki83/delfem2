/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_OPENGL_NEW_DRAWER_MSHCOLORMAP_H
#define DFM2_OPENGL_NEW_DRAWER_MSHCOLORMAP_H

#include <cstdio>
#include <vector>
#include <tuple>

#include "delfem2/dfm2_inline.h"
#include "delfem2/opengl/new/funcs.h" // CGL4_VAO_Mesh
#include "delfem2/color.h"

// -------------------------------------

namespace delfem2::opengl {

class Drawer_MeshColormap {
 public:
  Drawer_MeshColormap() = default;

  void AddConnectivity(
      std::vector<unsigned int> &aTri,
      int gl_primitive_type);

  template<typename REAL>
  void SetCoordinates(
      std::vector<REAL> &aPosD,
      unsigned int ndim);

  template<typename REAL>
  void SetValues(
      std::vector<REAL> &aValD);

  void InitGL();

  void Draw(
      const float mP[16],
      const float mMV[16],
      float val_min,
      float val_max) const;

 public:
  VertexArrayObject vao; // gl4
  int shaderProgram = -1;
  int Loc_MatrixProjection = -1;
  int Loc_MatrixModelView = -1;
  int Loc_UniformColor = -1;
  int Loc_UseUniformColor = -1;
  int Loc_ValMin = -1, Loc_ValMax = -1;
  std::vector<std::tuple<double, double, double>> colors = {
      {0,0,0},
      {1,1,1},
  };

};

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/new/drawer_mshcolormap.cpp"
#endif

#endif /* DFM2_OPENGL_NEW_DRAWER_MSHCOLORMAP_H */
