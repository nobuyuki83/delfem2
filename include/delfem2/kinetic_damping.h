//
// Created by Nobuyuki Umetani on 2022/02/16.
//

#ifndef DFM2_KINETIC_DAMPING_H_
#define DFM2_KINETIC_DAMPING_H_

#include <vector>

namespace delfem2 {

double EnergyKinetic(
    const double *aUVW,
    size_t np) {
  double E = 0.0;
  for (unsigned int ip = 0; ip < np; ++ip) {
    double u0 = aUVW[ip * 3 + 0];
    double v0 = aUVW[ip * 3 + 1];
    double w0 = aUVW[ip * 3 + 2];
    E += u0 * u0 + v0 * v0 + w0 * w0;
  }
  return E;
}

class CKineticDamper {
 public:
  void Damp(std::vector<double> &aUVW) {
    aEnergy.push_back(EnergyKinetic(aUVW.data(), aUVW.size() / 3));
    if (aEnergy.size() > 3) {
      aEnergy.erase(aEnergy.begin());
      const double g0 = aEnergy[1] - aEnergy[0];
      const double g1 = aEnergy[2] - aEnergy[1];
      if (g0 > 0 && g1 < 0) { aUVW.assign(aUVW.size(), 0.0); }
    }
  }
 public:
  std::vector<double> aEnergy;
};

}

#endif //DFM2_KINETIC_DAMPING_H_
