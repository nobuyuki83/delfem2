#ifndef SRCH_KDTREE_TEST_H_
#define SRCH_KDTREE_TEST_H_

#include <random>

#include "delfem2/srch_kdtree.h"

template <typename VEC, typename SCALAR = typename VEC::Scalar>
void TestKdtreeTest0(unsigned int nitr)
{
	std::mt19937 rngeng(std::random_device {}());
	std::uniform_real_distribution<SCALAR> dist_m1p1(-1, +1);
	for (int iitr = 0; iitr < nitr; iitr++) {
		std::vector<VEC> Points(1000);
		for (auto& p : Points) {
			p = { dist_m1p1(rngeng), dist_m1p1(rngeng), dist_m1p1(rngeng) };
		}

		VEC Query = { dist_m1p1(rngeng), dist_m1p1(rngeng), dist_m1p1(rngeng) };

		VEC ClosestBruteForse = Points[0];
		for (const auto v : Points) {
			if ((ClosestBruteForse - Query).norm() > (v - Query).norm())
				ClosestBruteForse = v;
		}

		delfem2::KDTNodePoint<VEC>* KDT = delfem2::ConstructKDT<VEC>(Points);
		VEC ClosestKDTree		= delfem2::SearchNearstPoint<VEC>(KDT, Query);

		EXPECT_EQ((ClosestBruteForse - ClosestKDTree).squaredNorm() < 0.00000001, true);
	}
}

#endif
