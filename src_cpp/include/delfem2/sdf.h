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
  virtual double Projection(double n[3],
                            double px, double py, double pz) const = 0; // normal
};

#endif
