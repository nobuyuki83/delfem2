//
//  gizmo_glold.h
//  examples_glfwold_hdronly
//
//  Created by Nobuyuki Umetani on 2020-04-19.
//

#ifndef DFM2_OPENGL_OLD_GIZMO_H
#define DFM2_OPENGL_OLD_GIZMO_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/gizmo_geo3.h"

namespace delfem2{
namespace opengl{

template <typename REAL>
DFM2_INLINE void DrawAxisHandler(
    REAL s,
    const CVec3<REAL>& p,
    int ielem_picked);

template <typename REAL>
DFM2_INLINE void DrawHandlerRotation_PosQuat(
    const CVec3<REAL>& pos,
    const REAL quat[4],
    REAL size,
    int ielem_picked);

template <typename REAL>
void Draw(
    const CGizmo_Rotation<REAL>& gizmo_rot);

template <typename REAL>
void Draw(
    const CGizmo_Transl<REAL>& gizmo_trnsl);


template <typename REAL>
void Draw(
    const CGizmo_Affine<REAL>& ga);

}
}

#ifndef DFM2_STATIC_LIBRARY
  #include "delfem2/opengl/old/gizmo.cpp"
#endif

#endif /* gizmo_glold_h */
