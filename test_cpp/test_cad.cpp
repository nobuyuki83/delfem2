/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <string>
#include <random>

#include "gtest/gtest.h"

#include "delfem2/cad2.h"
#include "delfem2/cad2_io_svg.h"
#include "delfem2/geo_bezier_cubic.h"

namespace dfm2 = delfem2;

TEST(cad,read_svg) {
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/shape0.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/shape1.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/shape2.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/tshirt.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/ltshirt.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/raglan.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/raglan2.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
}

