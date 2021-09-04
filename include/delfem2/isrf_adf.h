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
class CInput_ADF3 {
 public:
  virtual double sdf(double px, double py, double pz) const = 0;
};

/**
 * @brief Adaptive distance field
 */
class CADF3 {
 public:
  CADF3();
  ~CADF3();
  void SetUp(const CInput_ADF3 &ct, double bb[6]);
  void SetFaceColor(double r, double g, double b) {
    color_[0] = r;
    color_[1] = g;
    color_[2] = b;
  }
  virtual double Projection(
      double px, double py, double pz,
      double n[3]) const;
  virtual bool IntersectionPoint(
      [[maybe_unused]] double p[3],
      [[maybe_unused]] const double org[3], 
	  [[maybe_unused]] const double dir[3]) const { return false; }
  virtual void GetMesh(
	  [[maybe_unused]] std::vector<unsigned int> &aTri, 
	  [[maybe_unused]] std::vector<double> &aXYZ, 
	  [[maybe_unused]] double elen) const {}
  //
  void BuildIsoSurface_MarchingCube();
  void BuildMarchingCubeEdge();
  void SetShowCage(bool is_show) { this->is_show_cage = is_show; }
 public:
  class CNode {
   public:
    CNode();
    CNode(const CNode &no);
    void SetCornerDist(const CInput_ADF3 &ct);
    void MakeChildTree(const CInput_ADF3 &ct, std::vector<CNode> &aNo, double min_hw, double max_hw);
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
  unsigned int nIsoTri_;
  double *aIsoTri_;
  double *aIsoEdge_;
  bool is_show_cage;
  double color_[3];
};

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/isrf_adf.cpp"
#endif

#endif
