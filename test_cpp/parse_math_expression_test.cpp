/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstring>
#include <random>

#include "gtest/gtest.h"
#include "delfem2/parse_math_expression.h"

#ifndef M_PI
#  define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// ------------------------------------------


TEST(mathexpeval,test1){
  delfem2::CMathExpressionEvaluator e;
  e.SetKey("x", 3.0);
  e.SetExp("x+3.0");
  EXPECT_DOUBLE_EQ(e.Eval(),6);
  e.SetKey("x", 5.0);
  EXPECT_DOUBLE_EQ(e.Eval(),8.0);
  //
  e.SetKey("x", 1.0);
  e.SetKey("y", 2.0);
  e.SetExp("x+y");
  EXPECT_DOUBLE_EQ(e.Eval(),3.0);
  //
  e.SetExp("sin(PI*0.5*x)");
  EXPECT_DOUBLE_EQ(e.Eval(),1.0);
}