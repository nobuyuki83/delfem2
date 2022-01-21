//
// Created by Nobuyuki Umetani on 2022/01/17.
//

#ifndef GJK2_GEO_GJK2_H_
#define GJK2_GEO_GJK2_H_

namespace delfem2 {

// finds furthest point from poly in direction dir
template<typename VEC, typename SCALAR = typename VEC::Scalar>
VEC FindFurthest(
    const std::vector<VEC> &poly,
    const VEC &dir) {
  SCALAR max_dist = std::numeric_limits<SCALAR>::lowest();
  assert(!poly.empty());
  VEC ret = poly[0];
  for (auto vert: poly) {
    SCALAR dist = dir.dot(vert);
    if (dist < max_dist) { continue; }
    ret = vert;
    max_dist = dist;
  }
  return ret;
}

/**
 * The functions returns the triangle simplex that GJK finds if there is intersection
 * @tparam VEC
 * @param vtxsA points of the point set A
 * @param vtxsB points of the point set B
 */
template<typename VEC>
auto FindIntersectingSimplexForGjk2(
    const std::vector<VEC> &vtxsA,
    const std::vector<VEC> &vtxsB) -> std::vector<VEC> {
  assert(!vtxsA.empty() && !vtxsB.empty());

  // compute support on Minkowski difference
  auto support = [&vtxsA, &vtxsB](const VEC &dir) -> VEC {
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
      return std::vector<VEC>();
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
        return simplex;
      }
    }
  }
}

/**
 * The functions returns true if the convex hull of a 2D point set A
 * and the convex hull of a 2D point set B intersects.
 * @tparam VEC
 * @param vtxsA points of the point set A
 * @param vtxsB points of the point set B
 * @return
 */
template<typename VEC>
bool IsIntersect_Points2_Points2_Gjk(
    const std::vector<VEC> &vtxsA,
    const std::vector<VEC> &vtxsB) {
  assert(!vtxsA.empty() && !vtxsB.empty());
  return !FindIntersectingSimplexForGjk2(vtxsA, vtxsB).empty();
}

/**
 * @brief return normal vector towards origin of the closest edge and the starting index
 * @param simplex
 */
template<typename VEC, typename SCALAR = typename VEC::Scalar>
auto FindClosestEdge(
    const std::vector<VEC> &simplex) -> std::pair<unsigned int, VEC> {
  assert(!simplex.empty());
  SCALAR min_dist = std::numeric_limits<SCALAR>::max();
  unsigned int min_idx = UINT_MAX;
  VEC ret_normal;

  for (unsigned int i = 0; i < simplex.size(); i++) {
    unsigned int j = (i + 1 + simplex.size()) % simplex.size();
    VEC edge = simplex[j] - simplex[i];
    VEC n{edge[1], -edge[0]}; // we know the simplex is counterclockwise, so origin is always at left
    n.normalize();
    SCALAR dist = n.dot(simplex[i]);
    if (dist >= min_dist) { continue; }
    min_dist = dist;
    min_idx = i;
    ret_normal = n;
  }

  return {min_idx, ret_normal};
}

/**
 * computing maximum penetration depth and its normal for the intersection of convex hulls
 * @tparam VEC Eigen::Vector2x
 * @param[out] normalA if we move all the vertices of vtxB with normalA, there is no collision
 * @param[in] vtxsA coordinates of point set A
 * @param[in] vtxsB coordinates of point set B
 * @return Direction (Vector), depth
 */
template<typename VEC, typename SCALAR = typename VEC::Scalar>
bool Penetration_Points2_Points2_Epa(
    VEC &normalA,
    const std::vector<VEC> &vtxsA,
    const std::vector<VEC> &vtxsB,
    SCALAR tolerance) {
  std::vector<VEC> simplex = FindIntersectingSimplexForGjk2(vtxsA, vtxsB);
  if (simplex.empty()) { return false; } // no intersection
  assert(simplex.size() == 3);

  { // make the simplex counterclockwise
    const VEC v01 = simplex[1] - simplex[0];
    const VEC v02 = simplex[2] - simplex[0];
    if (v01[0] * v02[1] - v01[1] * v02[0] < 0) {  // check area of triangle
      const VEC temp = simplex[2];
      simplex[2] = simplex[1];
      simplex[1] = temp;
    }
  }

  auto support = [&vtxsA, &vtxsB](const VEC &dir) -> VEC {
    const VEC ndir = -dir;
    return FindFurthest(vtxsA, dir) - FindFurthest(vtxsB, ndir);
  };

  while (true) {
    const std::pair<int, VEC> ret = FindClosestEdge(simplex);
    int v0 = ret.first;
    const VEC n = ret.second;
    const SCALAR dist = n.dot(simplex[v0]);
    const VEC p = support(n);
    const SCALAR d = p.dot(n);
    if (d - dist < tolerance) {
      normalA = n * dist;
      return true;
    } else {
      simplex.insert(simplex.begin() + v0 + 1, p);
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
template<typename VEC, typename SCALAR=typename VEC::Scalar>
std::pair<SCALAR, SCALAR> Range_ProjectionPointsOnAxis(
    const std::vector<VEC> &vtxs,
    const VEC &a) {
  SCALAR min0 = a.dot(vtxs[0]);
  SCALAR max0 = min0;
  for (unsigned int ivtx = 1; ivtx < vtxs.size(); ++ivtx) {
    const SCALAR d = a.dot(vtxs[ivtx]);
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
      const auto rangeA = Range_ProjectionPointsOnAxis(vtxsA, a);
      const auto rangeB = Range_ProjectionPointsOnAxis(vtxsB, a);
      if (rangeA.second < rangeB.first) { return false; } // not intersect
      if (rangeB.second < rangeA.first) { return false; } // not intersect
    }
  }
  for (unsigned int iA = 0; iA < vtxsA.size(); ++iA) {
    for (unsigned int jA = iA + 1; jA < vtxsA.size(); ++jA) {
      VEC a = vtxsA[iA] - vtxsA[jA];
      a = {a[1], -a[0]}; // rotate 90 degree
      const auto rangeA = Range_ProjectionPointsOnAxis(vtxsA, a);
      const auto rangeB = Range_ProjectionPointsOnAxis(vtxsB, a);
      if (rangeA.second < rangeB.first) { return false; } // not intersect
      if (rangeB.second < rangeA.first) { return false; } // not intersect
    }
  }
  return true;
}

}

#endif //GJK2_GEO_GJK2_H_
