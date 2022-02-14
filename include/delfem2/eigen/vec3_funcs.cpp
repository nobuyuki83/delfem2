//
// Created by Nobuyuki Umetani on 2022/02/05.
//

// this file is for explicit instanciation
#undef DFM2_STATIC_LIBRARY
#include "delfem2/vec3_funcs.h"

#include <Eigen/Dense>

namespace delfem2 {
using f0 = Eigen::Vector3f;
using d0 = Eigen::Vector3d;
//
template float Dot3(const f0 &, const f0 &);
template double Dot3(const d0 &, const d0 &);
//
template float Distance3(const f0 &, const f0 &);
template double Distance3(const d0 &, const d0 &);
//
template float Length3(const f0 &);
template double Length3(const d0 &);
//
template void Normalize3(f0 &);
template void Normalize3(d0 &);
//
template f0 Cross(const f0 &, const f0 &);
template d0 Cross(const d0 &, const d0 &);
//
template void Cross(f0 &, const f0 &, const f0 &);
template void Cross(d0 &, const d0 &, const d0 &);
//
template float ScalarTripleProduct(const f0 &, const f0 &, const f0 &);
template double ScalarTripleProduct(const d0 &, const d0 &, const d0 &);
//
template f0 RotateVec3WithAxisAngleVector(const f0 &, const f0 &);
template d0 RotateVec3WithAxisAngleVector(const d0 &, const d0 &);
//
template void FrameFromVectorZ(f0 &, f0 &, const f0 &);
template void FrameFromVectorZ(d0 &, d0 &, const d0 &);
//
template std::ostream &operator<<(std::ostream &output, const std::vector<f0> &);
}


