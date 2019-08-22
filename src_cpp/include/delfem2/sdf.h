/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef SDF_H
#define SDF_H

class CSDF3
{
public:
  virtual ~CSDF3(){};
  // return signed distance function (inside is positive, outside is negative)
  // the normal direction (facing outward) is set to n[3]
  // "sdf * n" gives the nearest point on the surface
  virtual double Projection(double n[3],
                            double px, double py, double pz) const = 0; 
};

#endif
