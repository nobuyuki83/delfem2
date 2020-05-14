/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_RENDER2TEX_CEREAL_H
#define DFM2_RENDER2TEX_CEREAL_H

#include "cereal/cereal.hpp"
#include "cereal/archives/json.hpp"
#include "cereal/types/vector.hpp"
//
#include "delfem2/opengl/r2tglo_glold.h" // header itself doesn't depend on OpenGLa

namespace cereal{

template <class Archive>
void serialize
(Archive & ar,
 delfem2::opengl::CRender2Tex& smplr)
{
  ar(cereal::make_nvp("nResX", smplr.nResX),
     cereal::make_nvp("nResY", smplr.nResY),
     cereal::make_nvp("is_rgba_8ui", smplr.is_rgba_8ui),
     cereal::make_nvp("lengrid", smplr.lengrid),
     cereal::make_nvp("z_range", smplr.z_range),
     cereal::make_nvp("z_axis", smplr.z_axis ),
     cereal::make_nvp("x_axis", smplr.x_axis ),
     cereal::make_nvp("origin", smplr.origin ));
}


template <class Archive>
void serialize
 (Archive & ar,
  delfem2::opengl::CRender2Tex_DrawOldGL& smplr)
{
  ar(cereal::make_nvp("draw_len_axis", smplr.draw_len_axis),
     cereal::make_nvp("aZ", smplr.aZ),
     cereal::make_nvp("aRGBA_8ui", smplr.aRGBA_8ui),
     cereal::make_nvp("aRGBA_32f", smplr.aRGBA_32f),
     cereal::make_nvp("pointSize", smplr.pointSize),
     cereal::make_nvp("isDrawTex", smplr.isDrawTex),
     cereal::make_nvp("isDrawOnlyHitPoints",smplr.isDrawOnlyHitPoints),
     cereal::make_nvp("colorPoint", smplr.colorPoint) );
  ar( cereal::make_nvp("CRender2Tex", (delfem2::opengl::CRender2Tex&)smplr) );
}

template <class Archive>
void serialize
(Archive & ar,
 delfem2::opengl::CRender2Tex_DrawOldGL_BOX& smplr_box)
{
  ar( cereal::make_nvp("aSampler", smplr_box.aSampler) );
}

} // cereal


#endif /* render2tex_h */
