//
// Created by Nobuyuki Umetani on 2021-08-20.
//

#ifndef METASHAPE_CAMERA_H_
#define METASHAPE_CAMERA_H_

#include <string>

#include "delfem2/cvcamera.h"
#include "delfem2/mat4.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"

namespace delfem2 {

class MetashapeCamConfig {
 public:
  class Camera {
   public:
    Camera() : has_transform(false) {}
   public:
    int id_sensor = -1;
    bool has_transform;
    delfem2::CMat4f transform;
    std::string label;
  };
  class Sensor {
   public:
    float f;
    float cx, cy;
    float k1, k2, k3;
    float p1, p2;
    unsigned int width, height;
  };
 public:
  std::vector<Sensor> aSensor;
  std::vector<Camera> aCamera;
  delfem2::CVec3f region_center;
  delfem2::CVec3f region_size;
  delfem2::CMat3f region_R;
};

std::pair<delfem2::CMat4f, delfem2::CMat4f>
CameraTransformationMVP(
    const MetashapeCamConfig &cam_config,
    unsigned int idx_camera) {
  assert(idx_camera < cam_config.aCamera.size());
  const MetashapeCamConfig::Camera &cam = cam_config.aCamera[idx_camera];
  assert(cam.id_sensor >= 0 && cam.id_sensor < static_cast<int>(cam_config.aSensor.size()));
  const MetashapeCamConfig::Sensor &snsr = cam_config.aSensor[cam.id_sensor];
  const auto width0 = static_cast<float>(snsr.width);
  const auto height0 = static_cast<float>(snsr.height);
  const delfem2::CMat4f glb2img(
      delfem2::Mat4_CameraInternal_MetashapePinhole(
          snsr.f, snsr.cx, snsr.cy, width0, height0).data());
  const delfem2::CMat4f img2scr(
      delfem2::Mat4_Image2Screen(
          width0, height0, -0.1).data());
  delfem2::CMat4f zinv(
      1, 0, 0, 0,
      0, -1, 0, 0,
      0, 0, -1, 0,
      0, 0, 0, 1);
  return {zinv * cam.transform.Inverse(), img2scr * glb2img * zinv};
//  return {cam.transform.Inverse(), img2scr * glb2img};
}

unsigned int SnapCamera(
    const delfem2::CMat4f &model_view_matrix,
    const MetashapeCamConfig &cam_config) {
  delfem2::CMat4f mv0 = model_view_matrix; // GetModelViewMatrix();
  delfem2::CVec3f dir0 = mv0.Inverse().MultVec3(std::array<float, 3>{0, 0, 1}.data());
  const float phi0 = std::asin(dir0.y); // [-pi/2, +pi/2]
  const float theta0 = std::atan2(dir0.x, dir0.z);  // [-pi,+pi]
  unsigned int icam_best = UINT_MAX;
  double dist_best = -1;
  for (unsigned int icam = 0; icam < cam_config.aCamera.size(); ++icam) {
    {  // skip inactive camera
      const int idx_sensor = cam_config.aCamera[icam].id_sensor;
      if (idx_sensor < 0 || idx_sensor >= static_cast<int>(cam_config.aSensor.size())) {
        continue;
      }
    }
    auto[mv1, mp] = CameraTransformationMVP(cam_config, icam);
    delfem2::CVec3f dir1 = mv1.Inverse().MultVec3(std::array<float, 3>{0, 0, 1}.data());
    const float phi1 = std::asin(dir1.y); // [-pi/2, +pi/2]
    const float theta1 = std::atan2(dir1.x, dir1.z);  // [-pi,+pi]
    double dt0 = theta0 - theta1;
    double dt1 = theta0 - theta1 + 2 * M_PI;
    double dt2 = theta0 - theta1 - 2 * M_PI;
    double dist0 = (phi0 - phi1) * (phi0 - phi1) + dt0 * dt0;
    double dist1 = (phi0 - phi1) * (phi0 - phi1) + dt1 * dt1;
    double dist2 = (phi0 - phi1) * (phi0 - phi1) + dt2 * dt2;
    if (dist0 < dist_best || dist_best < 0) {
      icam_best = icam;
      dist_best = dist0;
    }
    if (dist1 < dist_best || dist_best < 0) {
      icam_best = icam;
      dist_best = dist1;
    }
    if (dist2 < dist_best || dist_best < 0) {
      icam_best = icam;
      dist_best = dist2;
    }
  }
  return icam_best;
}

}

#endif /* METASHAPE_CAMERA_H_ */
