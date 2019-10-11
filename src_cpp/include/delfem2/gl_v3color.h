/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef funcs_glcolorvec3_hpp
#define funcs_glcolorvec3_hpp

#include "delfem2/gl2_color.h"
#include "delfem2/vec3.h"

#include <stdio.h>

void DrawQuad_ScalarQ1
(const CVector3& p0, const CVector3& p1, const CVector3& p2, const CVector3& p3,
 double v0, double v1, double v2, double v3,
 const std::vector<std::pair<double, CColor> >& colorMap);

#endif /* funcs_glcolorvec3_hpp */
