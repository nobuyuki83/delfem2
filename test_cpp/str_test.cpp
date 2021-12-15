/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "gtest/gtest.h"
#include "delfem2/str.h"
#include "delfem2/mshmisc.h"

#ifndef M_PI
#  define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// ------------------------------------------

TEST(str,split_parentheses){
  {
    std::string str = "(a,b),c,(d,e)";
    std::vector<std::string> aS = dfm2::Split_Parentheses(str, ',', "()");
    EXPECT_EQ(aS.size(), 3);
    EXPECT_EQ(aS[0],"(a,b)");
    EXPECT_EQ(aS[1],"c");
    EXPECT_EQ(aS[2],"(d,e)");
  }
  {
    std::string str = "(a,b),c";
    std::vector<std::string> aS = dfm2::Split_Parentheses(str, ',', "()");
    EXPECT_EQ(aS.size(), 2);
    EXPECT_EQ(aS[0],"(a,b)");
    EXPECT_EQ(aS[1],"c");
  }
  {
    std::string str = "a,(b,c)";
    std::vector<std::string> aS = dfm2::Split_Parentheses(str, ',', "()");
    EXPECT_EQ(aS.size(), 2);
    EXPECT_EQ(aS[0],"a");
    EXPECT_EQ(aS[1],"(b,c)");
  }
}

TEST(str,split_quote){
  {
    std::string str = R"("a,b",c,"d,e")";
    std::vector<std::string> aS = dfm2::Split_Quote(str, ',', '\"' );
    EXPECT_EQ(aS.size(), 3);
    EXPECT_EQ(aS[0],"\"a,b\"");
    EXPECT_EQ(aS[1],"c");
    EXPECT_EQ(aS[2],"\"d,e\"");
  }
  {
    std::string str = R"("a,b",,c,"d,e")";
    std::vector<std::string> aS = dfm2::Split_Quote(str, ',', '\"' );
    EXPECT_EQ(aS.size(), 3);
    EXPECT_EQ(aS[0],"\"a,b\"");
    EXPECT_EQ(aS[1],"c");
    EXPECT_EQ(aS[2],"\"d,e\"");
  }
}

TEST(str,split){
  std::vector<std::string> aToken;
  aToken = dfm2::Split("chr count=80"," =");
  EXPECT_EQ(aToken.size(), 3);
  EXPECT_EQ(aToken[0],"chr");
  EXPECT_EQ(aToken[1],"count");
  EXPECT_EQ(aToken[2],"80");
  //
  aToken = dfm2::Split("chr = 80"," =");
  EXPECT_EQ(aToken.size(), 2);
  EXPECT_EQ(aToken[0],"chr");
  EXPECT_EQ(aToken[1],"80");
  //
  aToken = dfm2::Split("=chr = 80="," =");
  EXPECT_EQ(aToken.size(), 2);
  EXPECT_EQ(aToken[0],"chr");
  EXPECT_EQ(aToken[1],"80");
}