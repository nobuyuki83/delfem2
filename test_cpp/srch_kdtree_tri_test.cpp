#include <random>

#include "gtest/gtest.h"

#include "delfem2/vec3.h"
#include "delfem2/geo_tri.h"
#include "delfem2/msh_primitive.h"
#include "delfem2/msh_affine_transformation.h"

#include "delfem2/srch_kdtree_tri.h"

TEST(kdtree_tri, test0)
{
	using VEC    = delfem2::CVec3d;
	using SCALAR = double;

	std::mt19937 rngeng(std::random_device {}());
	std::uniform_real_distribution<SCALAR> dist_m1p1(-1, +1);

	std::vector<double> XYZ;
	std::vector<unsigned int> Tri;

	delfem2::MeshTri3D_Sphere(XYZ, Tri, 1.0, 64, 64);
	delfem2::Rotate_Points3(XYZ, 0.1, 0.2, 0.3);

	std::vector<VEC> Points(XYZ.size() / 3);
	for (unsigned int i = 0; i < XYZ.size() / 3; i++) {
		Points[i] = VEC { XYZ[3 * i + 0], XYZ[3 * i + 1], XYZ[3 * i + 2] };
	}

	auto KDTree = delfem2::KDTTriangle<VEC>(Points, Tri);

	for (int iitr = 0; iitr < 1000; iitr++) {
		const VEC Query = { dist_m1p1(rngeng), dist_m1p1(rngeng), dist_m1p1(rngeng) };

		VEC ClosestBruteForce = Points[Tri[0]];
		for (unsigned int i = 0; i < Tri.size() / 3; i++) {
			const VEC v0 = Points[Tri[3 * i + 0]] - Query;
			const VEC v1 = Points[Tri[3 * i + 1]] - Query;
			const VEC v2 = Points[Tri[3 * i + 2]] - Query;

			SCALAR hoge0, hoge1;

			VEC temp = Nearest_Origin3_Tri3(hoge0, hoge1, v0, v1, v2);
			if (temp.squaredNorm() < (ClosestBruteForce - Query).squaredNorm())
				ClosestBruteForce = temp + Query;
		}

		auto Candidate	  = KDTree.SearchNearest(Query);
		VEC ClosestKDTree = Points[Tri[0]];
		for (const auto& TriInd : Candidate) {
			const VEC v0 = Points[Tri[3 * TriInd + 0]] - Query;
			const VEC v1 = Points[Tri[3 * TriInd + 1]] - Query;
			const VEC v2 = Points[Tri[3 * TriInd + 2]] - Query;

			SCALAR hoge0, hoge1;

			VEC temp = Nearest_Origin3_Tri3(hoge0, hoge1, v0, v1, v2);
			if (temp.squaredNorm() < (ClosestKDTree - Query).squaredNorm())
				ClosestKDTree = temp + Query;
		}

		EXPECT_EQ((ClosestBruteForce - ClosestKDTree).squaredNorm() < 0.0000001, true);
	}
}
