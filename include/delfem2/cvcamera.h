#ifndef CV_CAMERA_H
#define CV_CAMERA_H

#include <vector>
#include <fstream>

namespace delfem2 {

class CCvCamera {
public:
  char name[32];
  float t[3];
  float R[9];
  float K[9];
};

bool ReadCamera(
    std::vector<CCvCamera> &aCamera,
    const std::string &path_file) {
  std::ifstream fin(path_file);
  if (fin.fail()) { return false; }
  int ncam;
  fin >> ncam;
  aCamera.resize(ncam);
  for (int icam = 0; icam < ncam; ++icam) {
    int icam0;
    fin >> icam0;
    assert(icam0 == icam);
    fin >> aCamera[icam].name;
    for (float &k : aCamera[icam].K) { fin >> k; }
    for (float &t : aCamera[icam].t) { fin >> t; }
    for (float &r : aCamera[icam].R) { fin >> r; }
  }
  return true;
}

void SetCameraInteriorMat(
    float *K,
    unsigned int width,
    unsigned int height,
    float sensor_width,
    float focal_length) {
  const float pix_size = sensor_width / float(width);
  const float fcllenpix = focal_length / pix_size;
  K[0] = +fcllenpix;
  K[1] = 0.f;
  K[2] = -float(width) * 0.5f;
  //
  K[3] = 0.f;
  K[4] = -fcllenpix;
  K[5] = -float(height) * 0.5f;
  //
  K[6] = 0.f;
  K[7] = 0.f;
  K[8] = -1.f;
}

}

#endif