//
// Created by Nobuyuki Umetani on 2022/01/17.
//

#ifndef GJK2_GEO_GJK2_H_
#define GJK2_GEO_GJK2_H_

namespace delfem2 {

// finds furthest point from poly in direction dir
template <typename VEC>
auto FindFurthest(const std::vector<VEC>& poly, const VEC& dir) -> VEC {
  double max_dist = std::numeric_limits<double>::lowest();
  VEC ret = *poly.begin();
  for (auto vert : poly) {
    double dist = dir.dot(vert);
    if (dist > max_dist) {
      ret = vert;
      max_dist = dist;
    }
  }
  return ret;
}

/**
 * The functions returns true if the convex hull of a 2D point set A
 * and the convex hull of a 2D point set B intersects.
 * @tparam VEC
 * @param vtxsA points of the point set A
 * @param vtxsB points of the point set B
 * @return
 */
template <typename VEC>
bool IsIntersect_Points2_Points2_Gjk(
    const std::vector<VEC>& vtxsA,
    const std::vector<VEC>& vtxsB) {

  // compute support on Minkowski difference
  auto support = [&vtxsA, &vtxsB](const VEC& dir) -> VEC {
    VEC ndir = -dir;
    return FindFurthest(vtxsA, dir) - FindFurthest(vtxsB, ndir);
  };

  // add first and second points of simplex
  std::vector<VEC> simplex;
  VEC init_dir = vtxsB[0] - vtxsA[0];
  simplex.push_back(support(init_dir));
  VEC O = VEC{0, 0};
  VEC d = O - simplex[0];

  VEC AB, AC, AO, ABperp, ACperp; // all temporary variables
  while (true) {
    VEC P = support(d);
    if (P.dot(d) < 0) {
      return false;
    } else {
      simplex.push_back(P);
    }

    if (simplex.size() == 2) { // line case
      AB = simplex[0] - simplex[1];
      AO = O - simplex[1];
      d = AO - AB.dot(AO) / AB.squaredNorm() * AB; // ABperp
    } else if (simplex.size() == 3) { // triangle case
      AB = simplex[1] - simplex[2];
      AC = simplex[0] - simplex[2];
      AO = O - simplex[2];
      ABperp = -(AC - AB.dot(AC) / AB.squaredNorm() * AB);
      ACperp = -(AB - AC.dot(AB) / AC.squaredNorm() * AC);
      if (ABperp.dot(AO) > 0) {
        simplex.erase(simplex.begin()); // remove C
        d = ABperp;
      } else if (ACperp.dot(AO) > 0) {
        simplex.erase(simplex.begin() + 1); // remove B
        d = ACperp;
      } else {
        return true;
      }
    }
  }
}


/**
 *
 * @tparam VEC Eigen::Vector2, Eigen::Vector3, dfm2::CVec2
 * @param[in] vtxs points
 * @param[in] a axis of projection
 * @return range of the projection
 */
template<typename VEC>
std::pair<double, double> Range_ProjectionPointsOnAxis(
    const std::vector<VEC> &vtxs,
    const VEC &a) {
  double min0 = a.dot(vtxs[0]);
  double max0 = min0;
  for (unsigned int ivtx = 1; ivtx < vtxs.size(); ++ivtx) {
    double d = a.dot(vtxs[ivtx]);
    min0 = (d < min0) ? d : min0;
    max0 = (d > max0) ? d : max0;
  }
  assert(min0 <= max0);
  return {min0, max0};
}

template<typename VEC>
bool IsIntersect_Points2_Points2_Sat(
    const std::vector<VEC> &vtxsA,
    const std::vector<VEC> &vtxsB) {
  for (unsigned int iB = 0; iB < vtxsB.size(); ++iB) {
    for (unsigned int jB = iB + 1; jB < vtxsB.size(); ++jB) {
      VEC a = vtxsB[iB] - vtxsB[jB];
      a = {a[1], -a[0]}; // rotate 90 degree
      auto rangeA = Range_ProjectionPointsOnAxis(vtxsA, a);
      auto rangeB = Range_ProjectionPointsOnAxis(vtxsB, a);
      if (rangeA.second < rangeB.first) { return false; } // not intersect
      if (rangeB.second < rangeA.first) { return false; } // not intersect
    }
  }
  for (unsigned int iA = 0; iA < vtxsA.size(); ++iA) {
    for (unsigned int jA = iA + 1; jA < vtxsA.size(); ++jA) {
      VEC a = vtxsA[iA] - vtxsA[jA];
      a = {a[1], -a[0]}; // rotate 90 degree
      auto rangeA = Range_ProjectionPointsOnAxis(vtxsA, a);
      auto rangeB = Range_ProjectionPointsOnAxis(vtxsB, a);
      if (rangeA.second < rangeB.first) { return false; } // not intersect
      if (rangeB.second < rangeA.first) { return false; } // not intersect
    }
  }
  return true;
}



}

#endif //GJK2_GEO_GJK2_H_
