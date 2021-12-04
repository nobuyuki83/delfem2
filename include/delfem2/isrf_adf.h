/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_ISRF_ADF_H
#define DFM2_ISRF_ADF_H

#include <vector>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

/**
 * @brief virtual input class
 */
class Input_AdaptiveDistanceField3 {
 public:
  virtual double sdf(double px, double py, double pz) const = 0;
};

/**
 * @brief Adaptive distance field
 */
class AdaptiveDistanceField3 {
 public:
    AdaptiveDistanceField3() = default;
    
    ~AdaptiveDistanceField3() = default;
    
  void SetUp(const Input_AdaptiveDistanceField3 &ct, double bb[6]);
    
  virtual double Projection(
      double px, double py, double pz,
      double n[3]) const;
    
  void BuildIsoSurface_MarchingCube(std::vector<double> &aTri);
    
 public:
  class CNode {
   public:
    CNode();
    CNode(const CNode &) = default;
    void SetCornerDist(const Input_AdaptiveDistanceField3 &ct);
    void MakeChildTree(const Input_AdaptiveDistanceField3 &ct, std::vector<CNode> &aNo, double min_hw, double max_hw);
    double FindDistNormal
        (double px, double py, double pz,
         double n[3],
         const std::vector<CNode> &aNo) const;
    void GenerateIsoSurface
        (std::vector<double> &aTri,
         const std::vector<CNode> &aNo) const;
   public:
    double cent_[3];
    double hw_;
    int ichilds_[8];
    double dists_[8];
  };
 public:
  std::vector<CNode> aNode;
  double dist_min, dist_max;
};



} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/isrf_adf.cpp"
#endif

#endif
