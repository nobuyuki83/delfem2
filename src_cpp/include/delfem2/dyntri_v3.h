#ifndef DYNTRI_V3_H
#define DYNTRI_V3_H

#include <map>
#include <algorithm>
#include <stack>

#include "delfem2/vec3.h"
#include "delfem2/dyntri.h"

CVector3 normalTri(int itri0,
                   const std::vector<CEPo2>& aPo3D,
                   const std::vector<ETri>& aSTri);

bool CheckTri(const std::vector<CEPo2>& aPo3D,
              const std::vector<ETri>& aSTri,
              bool is_assert=true);

bool FindRayTriangleMeshIntersections(const CVector3 &line0,
                                      const CVector3 &line1,
                                      const std::vector<ETri>& aTri,
                                      const std::vector<CEPo2> &aPoint3D,
                                      std::vector<CVector3> &intersectionPoints);

//////////////////////////////////////////////////////////////////////////////////////////////////


#ifdef USE_GL

void DrawMeshDynTri_FaceNorm(const std::vector< CEPo2>& aPo3D,
                             const std::vector<ETri>& aSTri );

void DrawMeshDynTri_Edge(const std::vector< CEPo2>& aPo3D,
                         const std::vector<ETri>& aSTri );

#endif


#endif // #endif SURFACE_MESH_H
