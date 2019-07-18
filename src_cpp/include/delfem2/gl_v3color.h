#ifndef funcs_glcolorvec3_hpp
#define funcs_glcolorvec3_hpp

#include "delfem2/gl_color.h"
#include "delfem2/vec3.h"

#include <stdio.h>

void DrawQuad_ScalarQ1
(const CVector3& p0, const CVector3& p1, const CVector3& p2, const CVector3& p3,
 double v0, double v1, double v2, double v3,
 const std::vector<std::pair<double, CColor> >& colorMap);

#endif /* funcs_glcolorvec3_hpp */
