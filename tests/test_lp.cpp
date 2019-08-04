/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>

#include "gtest/gtest.h"
#include "delfem2/lp.h"

TEST(linpro,test1)
{
  // example in https://people.richland.edu/james/ictcm/2006/simplex.html
  // http://www.me.titech.ac.jp/~mizu_lab/text/PDF-LP/LP1-problem.pdf
  CLinPro lp;
  lp.AddEqn({ 1.0,  2.0}, 16.0, CLinPro::LE);
  lp.AddEqn({ 1.0,  1.0},  9.0, CLinPro::LE);
  lp.AddEqn({ 3.0,  2.0}, 24.0, CLinPro::LE);
  int nitr = 10;
  lp.Precomp(nitr);
  EXPECT_LT(nitr,10);
  //////////////////////////////
  std::vector<double> solution;
  double opt_val;
  nitr = 10;
  lp.Solve(solution,opt_val,nitr,
           {40.0, 30.0});
  EXPECT_LT(nitr,10);
  EXPECT_EQ(solution.size(),2);
  EXPECT_FLOAT_EQ(solution[0],   6);
  EXPECT_FLOAT_EQ(solution[1],   3);
  EXPECT_FLOAT_EQ(opt_val,     330);
}


TEST(linpro,test2)
{
  // example in http://www.bunkyo.ac.jp/~nemoto/lecture/or/99/simplex.pdf
  CLinPro lp;
  lp.AddEqn({ 1.0,  2.0},  800, CLinPro::LE);
  lp.AddEqn({ 3.0,  4.0}, 1800, CLinPro::LE);
  lp.AddEqn({ 3.0,  1.0}, 1500, CLinPro::LE);
  int nitr = 10;
  lp.Precomp(nitr);
  //////////////////////////////
  std::vector<double> solution;
  double opt_val;
  nitr = 10;
  lp.Solve(solution,opt_val,nitr,
           {20.0, 30.0});
  EXPECT_LT(nitr,10);
  EXPECT_EQ(solution.size(),2);
  EXPECT_FLOAT_EQ(solution[0],   200);
  EXPECT_FLOAT_EQ(solution[1],   300);
  EXPECT_FLOAT_EQ(opt_val,     13000);
}


TEST(linpro,test3)
{
  // example in https://people.richland.edu/james/ictcm/2006/simplex.html
  // http://www.me.titech.ac.jp/~mizu_lab/text/PDF-LP/LP1-problem.pdf
  CLinPro lp;
  lp.AddEqn({+0.0, +1.0}, +1.0, CLinPro::LE);
  lp.AddEqn({+1.0, +0.0}, +1.0, CLinPro::LE);
  lp.AddEqn({-2.0, -1.0}, -1.0, CLinPro::LE);
  int nitr = 10;
  lp.Precomp(nitr);
  EXPECT_LT(nitr, 10);
  std::vector<double> sol_ini;
  double opt_val;
  nitr = 10;
  lp.Solve(sol_ini,opt_val,nitr,
           {1.0, 1.0});
  EXPECT_LT(nitr, 10);
  EXPECT_EQ(sol_ini.size(),2);
  EXPECT_FLOAT_EQ(sol_ini[0], 1.0);
  EXPECT_FLOAT_EQ(sol_ini[1], 1.0);
  EXPECT_FLOAT_EQ(opt_val,    2.0);
}


TEST(linpro,test4)
{
  // http://zeus.mech.kyushu-u.ac.jp/~tsuji/java_edu/TwoPhase.html
  CLinPro lp;
  lp.AddEqn({+2.0, +1.0}, +8.0, CLinPro::GE);
  lp.AddEqn({+1.0, +1.0}, +6.0, CLinPro::GE);
  lp.AddEqn({+1.0, +2.0}, +8.0, CLinPro::GE);
  int nitr = 10;
  lp.Precomp(nitr);
  EXPECT_LT(nitr, 10);
  std::vector<double> sol_ini;
  double opt_val;
  nitr = 10;
  lp.Solve(sol_ini,opt_val,nitr,
           {-4.0, -3.0});
  EXPECT_LT(nitr, 10);
  EXPECT_EQ(sol_ini.size(),2);
  EXPECT_FLOAT_EQ(sol_ini[0], 2.0);
  EXPECT_FLOAT_EQ(sol_ini[1], 4.0);
  EXPECT_FLOAT_EQ(opt_val,    -20.0);
}


TEST(linpro,test5)
{
  // http://zeus.mech.kyushu-u.ac.jp/~tsuji/java_edu/TwoPhase.html
  CLinPro lp;
  lp.AddEqn({+1.0, +2.0, +0.0}, +12.0, CLinPro::EQ);
  lp.AddEqn({+1.0, +4.0, +3.0}, +20.0, CLinPro::EQ);
  int nitr = 10;
  lp.Precomp(nitr);
  EXPECT_LT(nitr, 10);
  std::vector<double> sol_ini;
  double opt_val;
  nitr = 10;
  lp.Solve(sol_ini,opt_val,nitr,
           {2.0, 1.0, 1.0});
  EXPECT_LT(nitr, 10);
  EXPECT_EQ(sol_ini.size(),3);
  EXPECT_FLOAT_EQ(sol_ini[0], 12.0);
  EXPECT_FLOAT_EQ(sol_ini[1], 0.0);
  EXPECT_FLOAT_EQ(sol_ini[2], 8.0/3.0);
  EXPECT_FLOAT_EQ(opt_val,    +80.0/3.0);
}

TEST(linpro,test6)
{
  // http://www.fujilab.dnj.ynu.ac.jp/lecture/system4.pdf
  CLinPro lp;
  lp.AddEqn({+1.0, +3.0}, +4.0, CLinPro::GE);
  lp.AddEqn({+2.0, +1.0}, +3.0, CLinPro::GE);
  int nitr = 10;
  bool res = lp.Precomp(nitr);
  EXPECT_TRUE(res);
  EXPECT_LT(nitr, 10);
  std::vector<double> sol_ini;
  double opt_val;
  nitr = 10;
  lp.Solve(sol_ini,opt_val,nitr,
           {-4.0, -1.0});
  EXPECT_LT(nitr, 10);
  EXPECT_EQ(sol_ini.size(),2);
  EXPECT_FLOAT_EQ(sol_ini[0], +0.0);
  EXPECT_FLOAT_EQ(sol_ini[1], +3.0);
  EXPECT_FLOAT_EQ(opt_val,    -3.0);
}

TEST(linpro,test7)
{
  // http://www.bunkyo.ac.jp/~nemoto/lecture/mathpro/2002/2stage-simplex.pdf
  CLinPro lp;
  lp.AddEqn({+2.0, +3.0},  +6.0, CLinPro::LE);
  lp.AddEqn({-5.0, +9.0}, +15.0, CLinPro::EQ);
  lp.AddEqn({-6.0, +3.0},  +3.0, CLinPro::GE);
  int nitr = 10;
  bool res = lp.Precomp(nitr);
  EXPECT_TRUE(res);
  EXPECT_LT(nitr, 10);
  std::vector<double> sol_ini;
  double opt_val;
  nitr = 10;
  lp.Solve(sol_ini,opt_val,nitr,
           {-6.0, +6.0});
  EXPECT_LT(nitr, 10);
  EXPECT_EQ(sol_ini.size(),2);
  EXPECT_FLOAT_EQ(sol_ini[0], +0.0);
  EXPECT_FLOAT_EQ(sol_ini[1], +5.0/3.0);
  EXPECT_FLOAT_EQ(opt_val,    +10);
}

TEST(linpro,test8)
{
  // http://www.bunkyo.ac.jp/~nemoto/lecture/mathpro/2002/2stage-simplex.pdf
  CLinPro lp;
  lp.AddEqn({-1.0, +1.0, +1.0}, +2.0, CLinPro::LE);
  lp.AddEqn({+2.0, +1.0, -1.0}, +8.0, CLinPro::EQ);
  lp.AddEqn({+1.0, +2.0, -1.0}, +1.0, CLinPro::GE);
  int nitr = 10;
  bool res = lp.Precomp(nitr);
  EXPECT_TRUE(res);
  EXPECT_LT(nitr, 10);
  std::vector<double> sol_ini;
  double opt_val;
  nitr = 10;
  lp.Solve(sol_ini,opt_val,nitr,
           {+1.0, +3.0, +5});
  EXPECT_LT(nitr, 10);
  EXPECT_EQ(sol_ini.size(),3);
  EXPECT_FLOAT_EQ(sol_ini[0], +8.0);
  EXPECT_FLOAT_EQ(sol_ini[1], +1.0);
  EXPECT_FLOAT_EQ(sol_ini[2], +9.0);
  EXPECT_FLOAT_EQ(opt_val,    +56.0);
}

TEST(linpro,test9)
{
  // http://www.bunkyo.ac.jp/~nemoto/lecture/mathpro/2002/2stage-simplex.pdf
  CLinPro lp;
  lp.AddEqn({+1.0, +1.0, +2.0}, +10.0, CLinPro::GE);
  lp.AddEqn({+3.0, +1.0, +1.0}, +20.0, CLinPro::GE);
  int nitr = 10;
  bool res = lp.Precomp(nitr);
  EXPECT_TRUE(res);
  EXPECT_LT(nitr, 10);
  std::vector<double> sol_ini;
  double opt_val;
  nitr = 10;
  lp.Solve(sol_ini,opt_val,nitr,
           {-12.0, -6.0, -10.0});
  EXPECT_LT(nitr, 10);
  EXPECT_EQ(sol_ini.size(),3);
  EXPECT_FLOAT_EQ(sol_ini[0], +5.0);
  EXPECT_FLOAT_EQ(sol_ini[1], +5.0);
  EXPECT_FLOAT_EQ(sol_ini[2], +0.0);
  EXPECT_FLOAT_EQ(opt_val,    -90.0);
}
