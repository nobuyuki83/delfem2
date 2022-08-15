#include <random>

#include "gtest/gtest.h"

#include "delfem2/vec3.h"
#include "geo_gjksat3_test.h"

TEST(gjksat3, test0)
{
	TestGjkSat3Test0<delfem2::CVec3d>(1000);
}
