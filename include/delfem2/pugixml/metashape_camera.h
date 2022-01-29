//
// Created by Nobuyuki Umetani on 2021-08-20.
//

#ifndef CAMERA_COFING_H_
#define CAMERA_COFING_H_

#include <sstream>
#include <optional>

#include "pugixml.hpp"

#include "delfem2/mat4.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"

#include "delfem2/metashape_camera.h"

namespace delfem2 {

namespace pugixml::metashape_camera {

std::vector<float> SplitTextIntoVectorFloat(const char *pc) {
  std::stringstream ss(pc);
  std::vector<float> vec;
  std::string buf_split;
  while (std::getline(ss, buf_split, ' ')) {
    vec.push_back(std::stof(buf_split));
  }
  return vec;
}

}

bool ReadMetashapeCameraXml(
    MetashapeCamConfig &cc,
    const char *path) {
  namespace lcl = pugixml::metashape_camera;
  pugi::xml_document file;
  const auto res = file.load_file(path);
  if (!res) {
    return false;
  }
  {  // region
    const auto &region = file.child("document").child("chunk").child("region");
    cc.region_center = delfem2::CVec3f(
        lcl::SplitTextIntoVectorFloat(region.child_value("center")).data());
    cc.region_size = delfem2::CVec3f(
        lcl::SplitTextIntoVectorFloat(region.child_value("size")).data());
    cc.region_R = delfem2::CMat3f(
        lcl::SplitTextIntoVectorFloat(region.child_value("R")).data());
  }
  for (const pugi::xml_node sensor: file.child("document").child("chunk").child("sensors")) {
    unsigned int id = sensor.attribute("id").as_int();
    if (cc.aSensor.size() < id + 1) {
      cc.aSensor.resize(id + 1);
    }
    const pugi::xml_node calibration = sensor.child("calibration");
    if (calibration.empty()) { continue; }  // this node does not have <calibration> tag
    cc.aSensor[id].width = calibration.child("resolution").attribute("width").as_uint();
    cc.aSensor[id].height = calibration.child("resolution").attribute("height").as_uint();
    cc.aSensor[id].f = std::stof(calibration.child_value("f"));
    {  // cx
      const pugi::char_t *sval = calibration.child_value("cx");
      cc.aSensor[id].cx = (*sval == '\0') ? 0.f : std::stof(sval);
    }
    {  // cy
      const pugi::char_t *sval = calibration.child_value("cy");
      cc.aSensor[id].cy = (*sval == '\0') ? 0.f : std::stof(sval);
    }
    cc.aSensor[id].k1 = std::stof(calibration.child_value("k1"));
    cc.aSensor[id].k2 = std::stof(calibration.child_value("k2"));
    {  // k3
      const pugi::char_t *sval = calibration.child_value("k3");
      cc.aSensor[id].k3 = (*sval == '\0') ? 0.f : std::stof(sval);
    }
    cc.aSensor[id].p1 = std::stof(calibration.child_value("p1"));
    cc.aSensor[id].p2 = std::stof(calibration.child_value("p2"));
  }
  for (pugi::xml_node camera: file.child("document").child("chunk").child("cameras")) {
    unsigned int id = camera.attribute("id").as_int();
    if (cc.aCamera.size() < id + 1) {
      cc.aCamera.resize(id + 1);
    }
    cc.aCamera[id].label = camera.attribute("label").as_string();
    cc.aCamera[id].id_sensor = camera.attribute("sensor_id").as_int();
    if (camera.child("transform").empty()) {
      continue;
    }
    std::vector<float> vec = lcl::SplitTextIntoVectorFloat(camera.child_value("transform"));
    assert(vec.size() == 16);
    cc.aCamera[id].transform = delfem2::CMat4f(vec.data());
    cc.aCamera[id].has_transform = true;
  }
  return true;
}

}

#endif
