/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_SRCHBV2AABB_H
#define DFM2_SRCHBV2AABB_H

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>

#include "delfem2/dfm2_inline.h"
#include "delfem2/vec2.h"

// -----------------------------------------------------

namespace delfem2 {

/**
 * @brief 2D bounding box
 * @details inactive if x_min > x_max
 */
template <typename REAL>
class CBoundingBox2
{
public:
  CBoundingBox2(){
    x_min=1;	x_max=-1;
    y_min=0;	y_max=0;
  }
  CBoundingBox2(double x_min0,double x_max0,  double y_min0,double y_max0)
  : x_min(x_min0),x_max(x_max0),  y_min(y_min0),y_max(y_max0)
  {}
  CBoundingBox2( const CBoundingBox2& bb ) = default;
  
  // -------------------------
  // const functions from here
  
  [[nodiscard]] bool isActive() const { return x_min <= x_max; }
  [[nodiscard]] bool IsIntersectSphere(const CVec2<double>& vec, const double radius ) const
  {
    if( !isActive() ) return false;
    if( vec.p[0] < x_min-radius || vec.p[0] > x_max+radius ||
       vec.p[1] < y_min-radius || vec.p[1] > y_max+radius ) return false;
    return true;
  }
  [[nodiscard]] bool IsIntersect(const CBoundingBox2& bb_j, double clearance) const
  {
    if( bb_j.x_min > x_max+clearance || bb_j.x_max < x_min-clearance ) return false;
    if( bb_j.y_min > y_max+clearance || bb_j.y_max < y_min-clearance ) return false;
    return true;
  }
  [[nodiscard]] std::vector<double> MinMaxXYZ() const {
    const double tmp[6] = {x_min,x_max, y_min,y_max, 0.0, 0.0};
    std::vector<double> bb(tmp,tmp+6);
    return bb;
  }
  [[nodiscard]] double LengthDiagonal() const {
    return sqrt( (x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) );
  }
  
  // -------------------------------
  // non const functions from here
  
  void Add(double x0, double y0){
    if( !isActive() ){
      x_min = x_max = x0;
      y_min = y_max = y0;
      return;
    }
    x_max = ( x_max > x0 ) ? x_max : x0;
    x_min = ( x_min < x0 ) ? x_min : x0;
    y_max = ( y_max > y0 ) ? y_max : y0;
    y_min = ( y_min < y0 ) ? y_min : y0;
  }
  CBoundingBox2& operator+=(const CBoundingBox2& bb)
  {
    if( !bb.isActive() ){ return *this; }
    if( !isActive() ){
      x_max = bb.x_max;	x_min = bb.x_min;
      y_max = bb.y_max;	y_min = bb.y_min;
      return *this;
    }
    x_max = ( x_max > bb.x_max ) ? x_max : bb.x_max;
    x_min = ( x_min < bb.x_min ) ? x_min : bb.x_min;
    y_max = ( y_max > bb.y_max ) ? y_max : bb.y_max;
    y_min = ( y_min < bb.y_min ) ? y_min : bb.y_min;
    return *this;
  }
  /*
  bool IsInside(const CVec2<double>& vec)
  {
    if( !isActive() ) return false;
    if(   vec.p[0] >= x_min && vec.p[0] <= x_max
       && vec.p[1] >= y_min && vec.p[1] <= y_max ) return true;
    return false;
  }
   */
public:
  REAL x_min,x_max,  y_min,y_max;
};
    
}

#endif /* DFM2_SRCHBV2AABB_H */


