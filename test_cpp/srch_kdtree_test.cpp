#include <random>

#include "gtest/gtest.h"

#include "delfem2/vec3.h"
#include "srch_kdtree_test.hpp"

TEST(kdtree, test0)
{
	TestKdtreeTest0<delfem2::CVec3d>(20);
}
