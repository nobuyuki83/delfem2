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

/**
 * Based on the Metashape's document:
 *   https://www.agisoft.com/pdf/metashape-pro_1_5_en.pdf
 * @tparam T
 * @param coord_image
 * @param pos_xyz
 * @param cx
 * @param cy
 * @param f
 * @param k1
 * @param k2
 * @param k3
 * @param p1
 * @param p2
 * @param img_width
 * @param img_height
 */
template <typename T>
void ProjectionCameraWithCalibration(
    T coord_image[2],
    const T pos_xyz[3],
    T cx, T cy, T f,
    T k1, T k2, T k3, T p1, T p2,
    unsigned int img_width,
    unsigned int img_height)
{
  T x0 = pos_xyz[0] / pos_xyz[2];
  T y0 = pos_xyz[1] / pos_xyz[2];
  T r1 = sqrt(x0 * x0 + y0 * y0);
  T r2 = r1 * r1;
  T r4 = r2 * r2;
  T r6 = r2 * r4;
  T t0 = 1. + k1 * r2 + k2 * r4 + k3 * r6;
  T x1 = x0 * t0 + p1 * (r2 + 2 * x0 * x0) + 2 * p2 * x0 * y0;
  T y1 = y0 * t0 + p2 * (r2 + 2 * y0 * y0) + 2 * p1 * x0 * y0;
  coord_image[0] = img_width * 0.5 + cx + f * x1;
  coord_image[1] = img_height * 0.5 + cy + f * y1;
}

}

#endif