//
// Created by Nobuyuki Umetani on 2021-10-30.
//

#ifndef DFM2_RIG_BVH_H_
#define DFM2_RIG_BVH_H_

#include "delfem2/rig_geo3.h"

namespace delfem2 {

class CChannel_BioVisionHierarchy {
 public:
  CChannel_BioVisionHierarchy(unsigned int ib, int ia, bool br)
  : ibone(ib), iaxis(ia), isrot(br) {}

 public:
  unsigned int ibone;
  int iaxis;
  bool isrot;
};

DFM2_INLINE void Read_BioVisionHierarchy(
    std::vector<CRigBone> &bones,
    std::vector<CChannel_BioVisionHierarchy> &channels,
    int &nframe,
    std::vector<double> &frame_channel,
    const std::string &path_bvh);

/**
 * @brief set value to CRigBone.rot (bone rotation from parent bone)
 */
DFM2_INLINE void SetPose_BioVisionHierarchy(
    std::vector<CRigBone> &bones,
    const std::vector<CChannel_BioVisionHierarchy> &channels,
    const double *values);

class BioVisionHierarchy {
 public:
  BioVisionHierarchy() = default;

  void Open(const std::string &file_path) {
    delfem2::Read_BioVisionHierarchy(
        bones, channels, nframe, frame_channel,
        file_path);
    assert(frame_channel.size() == nframe * channels.size());
  }

  void SetFrame(int iframe) {
    if (iframe < 0 || iframe >= nframe) { return; }
    const size_t nch = channels.size();
    delfem2::SetPose_BioVisionHierarchy(
        bones, channels,
        frame_channel.data() + nch * iframe);
  }

  void ClearPose() {
    const size_t nch = channels.size();
    std::vector<double> val(nch, 0.0);
    delfem2::SetPose_BioVisionHierarchy(
        bones, channels,
        val.data());
  }

  [[nodiscard]] std::vector<double> MinMaxXYZ() const {
    std::vector<double> bb{1., 1., 1., -1., -1., -1.};
    for (const auto &bone: bones) {
      const auto p = bone.RootPosition();
      UpdateMinMax(bb[0], bb[3], p[0]);
      UpdateMinMax(bb[1], bb[4], p[1]);
      UpdateMinMax(bb[2], bb[5], p[2]);
    }
    return bb;
  }

  [[nodiscard]] const delfem2::CRigBone &GetBone(unsigned int bone_idx) const {
    return bones[bone_idx];
  }

 private:
  template<typename T>
  void UpdateMinMax(
      T &vmin, T &vmax,
      T v) const {
    if (vmin > vmax) {
      vmin = vmax = v;
      return;
    }
    vmin = (vmin < v) ? vmin : v;
    vmax = (vmax > v) ? vmax : v;
  }

 public:
  std::vector<delfem2::CRigBone> bones;
  std::vector<delfem2::CChannel_BioVisionHierarchy> channels;
  int nframe = 0;
  std::vector<double> frame_channel;
};

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/rig_bvh.cpp"
#endif

#endif //DFM2_RIG_BVH_H_
