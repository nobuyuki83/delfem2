/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <vector>

class CInputIsosurfaceStuffing
{
public:
  virtual double SignedDistance(double px, double py, double pz) const = 0;
  
  virtual void Level(int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
                     double px, double py, double pz) const
  {
    sdf = this->SignedDistance(px, py, pz);
    ilevel_vol = -1;
    ilevel_srf = -1;
    nlayer = 1;
  }
};


bool IsoSurfaceStuffing
(std::vector<double>& aXYZ,
 std::vector<unsigned int>& aTet,
 std::vector<int>& aIsOnSurfXYZ,
 ////
 const CInputIsosurfaceStuffing& input,
 double elen_in,
 double width,
 const double center[3]);

class CPointLattice
{
public:
  CPointLattice(){
    pos[0] = 0;
    pos[1] = 0;
    pos[2] = 0;
    sdf = 0;
  }
  CPointLattice(double px, double py, double pz, double sdf){
    pos[0] = px;
    pos[1] = py;
    pos[2] = pz;
    this->sdf = sdf;
  }
public:
  double pos[3];
  double sdf;
  int iflg;
};

void makeBackgroundLattice(std::vector<CPointLattice>& aPoint, std::vector<unsigned int>& aTet,
                           const CInputIsosurfaceStuffing& input,
                           double elen, int  ndiv, const double org[3]);

