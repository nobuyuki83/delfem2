#ifndef GJK3_GEO_GJK3_H_
#define GJK3_GEO_GJK3_H_

#include <cassert>
#include <vector>
#include <limits>
#include <climits>
#include <utility>
#include <cmath>

namespace delfem2 {

template <typename VEC, typename SCALAR = typename VEC::Scalar>
VEC FindFurthest(const std::vector<VEC>& poly, const VEC& dir)
{
	SCALAR max_dist = std::numeric_limits<SCALAR>::lowest();
	assert(!poly.empty());
	VEC ret = poly[0];
	for (const auto& vert : poly) {
		SCALAR dist = dir.dot(vert);
		if (dist < max_dist)
			continue;
		ret	 = vert;
		max_dist = dist;
	}
	return ret;
}

/**
 * The functions returns the tetrahedra(tetrahedron) simplex that GJK finds if there is intersection
 * @tparam VEC
 * @param vtxsA points of the point set A
 * @param vtxsB points of the point set B
 */
template <typename VEC>
auto FindIntersectingSimplexForGjk3(const std::vector<VEC>& vtxsA, const std::vector<VEC>& vtxsB) -> std::vector<VEC>
{
	assert(!vtxsA.empty() && !vtxsB.empty());

	// compute support on Minkowski difference
	auto support = [&vtxsA, &vtxsB](const VEC& dir) -> VEC {
		VEC ndir = -dir;
		return FindFurthest(vtxsA, dir) - FindFurthest(vtxsB, ndir);
	};

	std::vector<VEC> simplex;

	{ // add first point of simplex
		VEC init_dir = vtxsB[0] - vtxsA[0];
		simplex.push_back(support(init_dir));
	}

	const VEC O = VEC { 0.0, 0.0, 0.0 };
	VEC d	    = O - simplex[0];

	{ // add second point
		VEC P = support(d);
		simplex.push_back(P);
	}

	{ // add thierd point
		VEC AB = simplex[0] - simplex[1];
		VEC AO = O - simplex[1];
		d      = AO - AB.dot(AO) / AB.squaredNorm() * AB; // ABperp
		VEC P  = support(d);
		simplex.push_back(P);
	}

	{ //set direction as the triangle normal (looking Origin)
		VEC AB = simplex[0] - simplex[1];
		VEC AC = simplex[0] - simplex[2];
		d      = AB.cross(AC);
		if ((O - simplex[0]).dot(d) < 0.0)
			d = -d;
	}

	VEC AB, AC, AD, AO, ABCperp, ACDperp, ADBperp; // all temporary variables
	while (true) {
		VEC P = support(d);
		if (P.dot(d) < 0) { // Minkowski difference does not contain Origin
			return std::vector<VEC>();
		} else {
			simplex.push_back(P);
		}

		AB = simplex[2] - simplex[3];
		AC = simplex[1] - simplex[3];
		AD = simplex[0] - simplex[3];
		AO = O - simplex[3];

		if (std::abs(AB.dot(AC.cross(AD))) < 0.000001) {
			return std::vector<VEC>();
		}

		// compute normal of each faces (looking outside of the tetrahedron)
		ABCperp = AB.cross(AC);
		if (ABCperp.dot(AD) > 0.0)
			ABCperp = -ABCperp;
		ACDperp = AC.cross(AD);
		if (ACDperp.dot(AB) > 0.0)
			ACDperp = -ACDperp;
		ADBperp = AD.cross(AB);
		if (ADBperp.dot(AC) > 0.0)
			ADBperp = -ADBperp;

		if (ABCperp.dot(AO) > 0) {
			simplex.erase(simplex.begin()); // remove D
			d = ABCperp;
		} else if (ACDperp.dot(AO) > 0) {
			simplex.erase(simplex.begin() + 2); // remove B
			d = ACDperp;
		} else if (ADBperp.dot(AO) > 0) {
			simplex.erase(simplex.begin() + 1); // remove C
			d = ADBperp;
		} else { // the tetrahedron contains Origin
			return simplex;
		}
	}
}

/**
 * The functions returns true if the convex hull of a 3D point set A
 * and the convex hull of a 3D point set B intersects.
 * @tparam VEC
 * @param vtxsA points of the point set A
 * @param vtxsB points of the point set B
 * @return
 */
template <typename VEC>
bool IsIntersect_Points3_Points3_Gjk(
    const std::vector<VEC>& vtxsA,
    const std::vector<VEC>& vtxsB)
{
	assert(!vtxsA.empty() && !vtxsB.empty());
	return !FindIntersectingSimplexForGjk3(vtxsA, vtxsB).empty();
}

/**
 * @brief return normal vector towards origin of the closest triangle and the starting index (if it returns i, 3*i+0,3*i+1,3*i+2 make face)
 * @param simplex
 */
template <typename VEC, typename SCALAR = typename VEC::Scalar>
auto FindClosestFace(
    const std::vector<VEC>& simplex, const std::vector<unsigned int>& facelist) -> std::pair<unsigned int, VEC>
{
	assert(!simplex.empty());
	SCALAR min_dist	     = std::numeric_limits<SCALAR>::max();
	unsigned int min_idx = UINT_MAX;
	VEC ret_normal;

	for (unsigned int i = 0; i < facelist.size() / 3; i++) {
		VEC v0 = simplex[facelist[3 * i + 0]];
		VEC v1 = simplex[facelist[3 * i + 1]];
		VEC v2 = simplex[facelist[3 * i + 2]];

		VEC n = (v1 - v0).cross(v2 - v0);
		n.normalize();

		SCALAR dist = n.dot(v0);
		if (dist >= min_dist)
			continue;

		min_dist   = dist;
		min_idx	   = i;
		ret_normal = n;
	}

	return { min_idx, ret_normal };
}

/**
 * computing maximum penetration depth and its normal for the intersection of convex hulls
 * @tparam VEC Eigen::Vector3x
 * @param[out] normalA if we move all the vertices of vtxB with normalA, there is no collision
 * @param[in] vtxsA coordinates of point set A
 * @param[in] vtxsB coordinates of point set B
 * @return Direction (Vector), depth
 */
template <typename VEC, typename SCALAR = typename VEC::Scalar>
bool Penetration_Points3_Points3_Epa(
    VEC& normalA,
    const std::vector<VEC>& vtxsA,
    const std::vector<VEC>& vtxsB,
    SCALAR tolerance)
{
	//the order of VEC does not make sense.
	std::vector<VEC> simplex = FindIntersectingSimplexForGjk3(vtxsA, vtxsB);
	if (simplex.empty()) {
		return false;
	} // no intersection
	assert(simplex.size() == 4);

	{ //make sure that v01.cross(v02) looking outside

		const VEC v01 = simplex[1] - simplex[0];
		const VEC v02 = simplex[2] - simplex[0];
		const VEC v03 = simplex[3] - simplex[0];

		const VEC n012 = v01.cross(v02);

		if (n012.dot(v03) > 0.0) {
			const VEC temp = simplex[3];
			simplex[3]     = simplex[0];
			simplex[0]     = temp;
		}
	}

	std::vector<unsigned int> facelist;
	facelist.reserve(3000);
	// facelist.size() = the number of faces x 3
	// vertex indies of the ith triangle are faceilst[3*i+0], facelist[3*i+1] and facelist[3*i+2] , where i is in [0,facesize/3]
	// each face look outside of polytope

	facelist.push_back(0);
	facelist.push_back(1);
	facelist.push_back(2);

	facelist.push_back(1);
	facelist.push_back(3);
	facelist.push_back(2);

	facelist.push_back(2);
	facelist.push_back(3);
	facelist.push_back(0);

	facelist.push_back(3);
	facelist.push_back(1);
	facelist.push_back(0);

	//temporary buffers
	std::vector<unsigned int> facebuffer;
	facebuffer.reserve(3000);
	std::vector<unsigned int> edgebuffer0;
	edgebuffer0.reserve(2000);
	std::vector<unsigned int> edgebuffer1;
	edgebuffer1.reserve(2000);

	auto support = [&vtxsA, &vtxsB](const VEC& dir) -> VEC {
		const VEC ndir = -dir;
		return FindFurthest(vtxsA, dir) - FindFurthest(vtxsB, ndir);
	};

	for (unsigned int iteration_counter = 0; iteration_counter < 1000; iteration_counter++) {

		const std::pair<unsigned int, VEC> ret = FindClosestFace(simplex, facelist);
		const unsigned int faceindex	       = ret.first;
		const VEC n			       = ret.second;
		const SCALAR dist		       = n.dot(simplex[facelist[3 * faceindex + 0]]);
		const VEC p			       = support(n);
		const SCALAR d			       = p.dot(n);
		if (d - dist < tolerance) {
			normalA = n * dist;
			return true;
		} else {

			// clear buffers
			facebuffer.clear();
			edgebuffer0.clear();
			edgebuffer1.clear();

			// check whether each face looks p or not
			// if the face looks p, it will be removed from polytope
			for (unsigned int i = 0; i < facelist.size() / 3; i++) {
				const VEC vert0 = simplex[facelist[3 * i + 0]];
				const VEC vert1 = simplex[facelist[3 * i + 1]];
				const VEC vert2 = simplex[facelist[3 * i + 2]];

				const VEC fn = (vert1 - vert0).cross(vert2 - vert0);

				if (fn.dot(p - vert0) < 0.0) {
					// the face will be contained in next polytope

					facebuffer.push_back(facelist[3 * i + 0]);
					facebuffer.push_back(facelist[3 * i + 1]);
					facebuffer.push_back(facelist[3 * i + 2]);
				} else {
					// the face will be removed
					// collect edges of removed faces

					edgebuffer0.push_back(facelist[3 * i + 0]);
					edgebuffer0.push_back(facelist[3 * i + 1]);

					edgebuffer0.push_back(facelist[3 * i + 1]);
					edgebuffer0.push_back(facelist[3 * i + 2]);

					edgebuffer0.push_back(facelist[3 * i + 2]);
					edgebuffer0.push_back(facelist[3 * i + 0]);
				}
			}

			// remove inner edge
			for (unsigned int i = 0; i < edgebuffer0.size() / 2; i++) {

				bool is_inner = false;

				// check whether the edge is shared by two removed faces (inner edge)
				for (unsigned int j = 0; j < edgebuffer0.size() / 2; j++) {
					if (i != j) {
						if ((edgebuffer0[2 * i + 0] == edgebuffer0[2 * j + 0] && edgebuffer0[2 * i + 1] == edgebuffer0[2 * j + 1]) || (edgebuffer0[2 * i + 1] == edgebuffer0[2 * j + 0] && edgebuffer0[2 * i + 0] == edgebuffer0[2 * j + 1])) {
							is_inner = true;
						}
					}
				}

				if (!is_inner) {
					edgebuffer1.push_back(edgebuffer0[2 * i + 0]);
					edgebuffer1.push_back(edgebuffer0[2 * i + 1]);
				}
			}

			//add new vertex
			simplex.push_back(p);

			const unsigned int vinew = simplex.size() - 1;

			for (unsigned int i = 0; i < edgebuffer1.size() / 2; i++) {
				const unsigned int vi0 = edgebuffer1[2 * i + 0];
				const unsigned int vi1 = edgebuffer1[2 * i + 1];

				facebuffer.push_back(vinew);
				facebuffer.push_back(vi0);
				facebuffer.push_back(vi1);
			}

			facelist = std::move(facebuffer);
		}
	}

	return false;
}

/**
 *
 * @tparam VEC Eigen::Vector3
 * @param[in] vtxs points
 * @param[in] a axis of projection
 * @return range of the projection
 */
template <typename VEC, typename SCALAR = typename VEC::Scalar>
std::pair<SCALAR, SCALAR> Range_ProjectionPointsOnAxis(
    const std::vector<VEC>& vtxs,
    const VEC& a)
{
	SCALAR min0 = a.dot(vtxs[0]);
	SCALAR max0 = min0;
	for (unsigned int ivtx = 1; ivtx < vtxs.size(); ++ivtx) {
		const SCALAR d = a.dot(vtxs[ivtx]);
		min0	       = (d < min0) ? d : min0;
		max0	       = (d > max0) ? d : max0;
	}
	assert(min0 <= max0);
	return { min0, max0 };
}

template <typename VEC>
bool IsIntersect_Points3_Points3_Sat(
    const std::vector<VEC>& vtxsA,
    const std::vector<unsigned int>& faceA,
    const std::vector<VEC>& vtxsB,
    const std::vector<unsigned int>& faceB)
{

	for (unsigned int fiB = 0; fiB < faceB.size() / 3; ++fiB) {
		const VEC vert0 = vtxsB[faceB[3 * fiB + 0]];
		const VEC vert1 = vtxsB[faceB[3 * fiB + 1]];
		const VEC vert2 = vtxsB[faceB[3 * fiB + 2]];

		const VEC normal = (vert1 - vert0).cross(vert2 - vert0);

		const auto rangeA = Range_ProjectionPointsOnAxis(vtxsA, normal);
		const auto rangeB = Range_ProjectionPointsOnAxis(vtxsB, normal);

		if (rangeA.second < rangeB.first) {
			return false;
		} // not intersect
		if (rangeB.second < rangeA.first) {
			return false;
		} // not intersect
	}

	for (unsigned int fiA = 0; fiA < faceA.size() / 3; ++fiA) {
		const VEC vert0 = vtxsA[faceA[3 * fiA + 0]];
		const VEC vert1 = vtxsA[faceA[3 * fiA + 1]];
		const VEC vert2 = vtxsA[faceA[3 * fiA + 2]];

		const VEC normal = (vert1 - vert0).cross(vert2 - vert0);

		const auto rangeA = Range_ProjectionPointsOnAxis(vtxsA, normal);
		const auto rangeB = Range_ProjectionPointsOnAxis(vtxsB, normal);

		if (rangeA.second < rangeB.first) {
			return false;
		} // not intersect
		if (rangeB.second < rangeA.first) {
			return false;
		} // not intersect
	}

	for (unsigned int fiB = 0; fiB < faceB.size() / 3; fiB++) {
		VEC edgeB[3];
		edgeB[0] = vtxsB[faceB[3 * fiB + 1]] - vtxsB[faceB[3 * fiB + 0]];
		edgeB[1] = vtxsB[faceB[3 * fiB + 2]] - vtxsB[faceB[3 * fiB + 0]];
		edgeB[2] = vtxsB[faceB[3 * fiB + 2]] - vtxsB[faceB[3 * fiB + 1]];
		for (unsigned int fiA = 0; fiA < faceA.size() / 3; fiA++) {
			VEC edgeA[3];
			edgeA[0] = vtxsA[faceA[3 * fiA + 1]] - vtxsA[faceA[3 * fiA + 0]];
			edgeA[1] = vtxsA[faceA[3 * fiA + 2]] - vtxsA[faceA[3 * fiA + 0]];
			edgeA[2] = vtxsA[faceA[3 * fiA + 2]] - vtxsA[faceA[3 * fiA + 1]];

			for (unsigned int i = 0; i < 6; i++) {

				constexpr unsigned int FaceToEdge[12] = {
					0, 0,
					0, 1,
					0, 2,
					1, 1,
					1, 2,
					2, 2
				};

				const VEC normal = (edgeB[FaceToEdge[2 * i + 1]]).cross(edgeA[FaceToEdge[2 * i + 0]]);

				const auto rangeA = Range_ProjectionPointsOnAxis(vtxsA, normal);
				const auto rangeB = Range_ProjectionPointsOnAxis(vtxsB, normal);

				if (rangeA.second < rangeB.first) {
					return false;
				} // not intersect
				if (rangeB.second < rangeA.first) {
					return false;
				} // not intersect
			}
		}
	}

	return true;
}

}

#endif //GJK3_GEO_GJK3_H_
