#ifndef GEO_GJKSAT3_TEST_H_
#define GEO_GJKSAT3_TEST_H_

#include <random>

#include "delfem2/geo_gjksat3.h"
#include "delfem2/geo_convhull3.h"

template <typename VEC, typename SCALAR = typename VEC::Scalar>
void TestGjkSat3Test0(unsigned int nitr)
{
	std::mt19937 rngeng(std::random_device {}());
	std::uniform_real_distribution<SCALAR> dist_m1p1(-1, +1);
	for (int iitr = 0; iitr < nitr; iitr++) {
		std::vector<VEC> vtxA_xyz(50);
		for (auto& p : vtxA_xyz) {
			p = { dist_m1p1(rngeng), dist_m1p1(rngeng), dist_m1p1(rngeng) };
		}
		std::vector<unsigned int> faceA_vtx;
		delfem2::ConvexHull3<VEC>(faceA_vtx, vtxA_xyz);
		//
		std::vector<VEC> vtxB_xyz(50);
		for (auto& p : vtxB_xyz) {
			p = { dist_m1p1(rngeng), dist_m1p1(rngeng), dist_m1p1(rngeng) };
		}
		std::vector<unsigned int> faceB_vtx;
		delfem2::ConvexHull3<VEC>(faceB_vtx, vtxB_xyz);

		const bool is_intersect_gjk = delfem2::IsIntersect_Points3_Points3_Gjk(vtxA_xyz, vtxB_xyz);
		EXPECT_EQ(is_intersect_gjk,
		    delfem2::IsIntersect_Points3_Points3_Sat(vtxA_xyz, faceA_vtx, vtxB_xyz, faceB_vtx));

		if (!is_intersect_gjk) {
			continue;
		}

		VEC normalA;
		delfem2::Penetration_Points3_Points3_Epa(normalA, vtxA_xyz, vtxB_xyz, 1.0e-5);

		std::vector<VEC> vtxB_xyz_translated;

		for (unsigned int ivtxB = 0; ivtxB < vtxB_xyz.size(); ++ivtxB) {
			vtxB_xyz_translated.push_back(vtxB_xyz[ivtxB] + normalA * 1.002);
		}
		EXPECT_EQ(delfem2::IsIntersect_Points3_Points3_Gjk(vtxA_xyz, vtxB_xyz_translated), false);

		vtxB_xyz_translated.clear();
		for (unsigned int ivtxB = 0; ivtxB < vtxB_xyz.size(); ++ivtxB) {
			vtxB_xyz_translated.push_back(vtxB_xyz[ivtxB] + normalA * 0.998);
		}
		EXPECT_EQ(delfem2::IsIntersect_Points3_Points3_Gjk(vtxA_xyz, vtxB_xyz_translated), true);
	}
}

#endif //GEO_GJKSAT2_TEST_H_
