/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "gtest/gtest.h"

#include "delfem2/mshmisc.h"
#include "delfem2/msh_io_obj.h"

#ifndef M_PI
#  define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// ------------------------------------------

TEST(mshio, load_obj) {
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  delfem2::Read_Obj(
      aXYZ, aTri,
      std::filesystem::path(PATH_INPUT_DIR) / "bunny_1k.obj");
  EXPECT_EQ(aTri.size(), 1000 * 3);
}
