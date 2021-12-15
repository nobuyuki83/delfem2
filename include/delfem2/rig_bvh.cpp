//
// Created by Nobuyuki Umetani on 2021-10-30.
//

#include "delfem2/rig_bvh.h"

#include <fstream>
#include <sstream>

namespace delfem2::rig_bvh {

DFM2_INLINE std::vector<std::string> MySplit(
  const std::string &str,
  char delimiter) {
  std::vector<std::string> aToken;
  aToken.clear();
  std::stringstream data(str);
  std::string line;
  while (std::getline(data, line, delimiter)) {
    if (line.empty()) { continue; }
    aToken.push_back(line);
  }
  return aToken;
}

DFM2_INLINE double myStod(const std::string &str) {
  char *e;
  double d = std::strtod(str.c_str(), &e);
  return d;
}

// probably std::stroi is safer to use but it is only for C++11
DFM2_INLINE int myStoi(const std::string &str) {
  char *e;
  long d = std::strtol(str.c_str(), &e, 0);
  return (int) d;
}

DFM2_INLINE std::string MyReplace(
  const std::string &str,
  const char cf,
  const char ct) {
  const size_t n = str.size();
  //
  std::string ss(str);
  for (unsigned int i = 0; i < n; ++i) {
    if (ss[i] != cf) { continue; }
    ss[i] = ct;
  }
  return ss;
}

DFM2_INLINE void CalcInvMat(double *a, const int n, int &info) {
  double tmp1;

  info = 0;
  int i, j, k;
  for (i = 0; i < n; i++) {
    if (fabs(a[i * n + i]) < 1.0e-30) {
      info = 1;
      return;
    }
    if (a[i * n + i] < 0.0) {
      info--;
    }
    tmp1 = 1.0 / a[i * n + i];
    a[i * n + i] = 1.0;
    for (k = 0; k < n; k++) {
      a[i * n + k] *= tmp1;
    }
    for (j = 0; j < n; j++) {
      if (j != i) {
        tmp1 = a[j * n + i];
        a[j * n + i] = 0.0;
        for (k = 0; k < n; k++) {
          a[j * n + k] -= tmp1 * a[i * n + k];
        }
      }
    }
  }
}

}


// ------------------------------------
// from here BioVisionHierarchy

DFM2_INLINE void delfem2::Read_BioVisionHierarchy(
  std::vector<CRigBone> &bones,
  std::vector<CChannel_BioVisionHierarchy> &channels,
  size_t &nframe,
  double &frame_time,
  std::vector<double> &frame_channel,
  std::string &bvh_header,
  const std::string &path_bvh) {
  std::ifstream fin;
  fin.open(path_bvh.c_str());
  if (!fin.is_open()) {
    std::cout << "cannot open file" << std::endl;
    return;
  }
  bvh_header.clear();
  bones.clear();
  channels.clear();
  //
  std::string line;
  std::vector<int> stackIndBone;
  while (std::getline(fin, line)) {
    if (line[line.size() - 1] == '\n') line.erase(line.size() - 1); // remove the newline code
    if (line[line.size() - 1] == '\r') line.erase(line.size() - 1); // remove the newline code
    bvh_header += (line + '\n');
    line = rig_bvh::MyReplace(line, '\t', ' ');
    std::vector<std::string> aToken = rig_bvh::MySplit(line, ' ');
//    std::cout << aToken[0] << std::endl;
    if (aToken[0] == "HIERARCHY") {
      assert(bones.empty());
    } else if (aToken[0] == "ROOT") {
      assert(bones.empty());
      CRigBone br;
      assert(aToken.size() == 2);
      br.name = aToken[1];
      bones.push_back(br);
    } else if (aToken[0] == "{") {
      stackIndBone.push_back(static_cast<int>(bones.size() - 1));
      if (stackIndBone.size() > 1) {
        int ibp = stackIndBone[stackIndBone.size() - 2];
        auto ib = static_cast<unsigned int>(bones.size() - 1);
        bones[ib].ibone_parent = ibp;
      }
    } else if (aToken[0] == "}") {
      stackIndBone.resize(stackIndBone.size() - 1);
    } else if (aToken[0] == "OFFSET") {
      assert(aToken.size() == 4);
      size_t ib = bones.size() - 1;
      double org_x = rig_bvh::myStod(aToken[1]);
      double org_y = rig_bvh::myStod(aToken[2]);
      double org_z = rig_bvh::myStod(aToken[3]);
      bones[ib].invBindMat[3] = -org_x;
      bones[ib].invBindMat[7] = -org_y;
      bones[ib].invBindMat[11] = -org_z;
      if (stackIndBone.size() > 1) {
        const int ibp = stackIndBone[stackIndBone.size() - 2];
        assert(ibp < (int) bones.size());
        bones[ib].invBindMat[3] += bones[ibp].invBindMat[3];
        bones[ib].invBindMat[7] += bones[ibp].invBindMat[7];
        bones[ib].invBindMat[11] += bones[ibp].invBindMat[11];
      }
    } else if (aToken[0] == "CHANNELS") {
      assert(aToken.size() >= 2);
      int nch = rig_bvh::myStoi(aToken[1]);
      assert((int) aToken.size() == nch + 2);
      assert(!bones.empty());
      const auto ib = static_cast<unsigned int>(bones.size() - 1);
      for (int ich = 0; ich < nch; ++ich) {
        const std::string &type_ch = aToken[ich + 2];
        if (type_ch == "Xposition") { channels.emplace_back(ib, 0, false); }
        else if (type_ch == "Yposition") { channels.emplace_back(ib, 1, false); }
        else if (type_ch == "Zposition") { channels.emplace_back(ib, 2, false); }
        else if (type_ch == "Xrotation") { channels.emplace_back(ib, 0, true); }
        else if (type_ch == "Yrotation") { channels.emplace_back(ib, 1, true); }
        else if (type_ch == "Zrotation") { channels.emplace_back(ib, 2, true); }
        else {
          std::cout << "ERROR-->undefiend type" << std::endl;
        }
      }
    } else if (aToken[0] == "JOINT") {
      CRigBone br;
      assert(aToken.size() == 2);
      br.name = aToken[1];
      bones.push_back(br);
    } else if (aToken[0] == "End") {
      assert(aToken[1] == "Site");
      CRigBone br;
      assert(aToken.size() == 2);
      br.name = aToken[1];
      bones.push_back(br);
    } else if (aToken[0] == "MOTION") {
      break;
    }
  }
  nframe = 0;
  {
    std::string stmp0;
    std::getline(fin, line);  // Frames: ***
    std::stringstream ss(line);
    ss >> stmp0 >> nframe;
    // std::cout << "frame: " << nframe << std::endl;
  }
  {
    std::string stmp0, stmp1;
    std::getline(fin, line);  // Frame Time: ***
    std::stringstream ss(line);
    ss >> stmp0 >> stmp1 >> frame_time;
    // std::cout << "frametime: " << frame_time << std::endl;
  }

  const size_t nchannel = channels.size();
  frame_channel.resize(nframe * nchannel);
  for (unsigned int iframe = 0; iframe < nframe; ++iframe) {
    std::getline(fin, line);
    line = rig_bvh::MyReplace(line, '\t', ' ');
    if (line[line.size() - 1] == '\n') line.erase(line.size() - 1); // remove the newline code
    if (line[line.size() - 1] == '\r') line.erase(line.size() - 1); // remove the newline code
    std::vector<std::string> aToken = rig_bvh::MySplit(line, ' ');
//    std::cout << aToken.size() << " " << aChannelRotTransBone.size() << std::endl;
    assert(aToken.size() == channels.size());
    for (unsigned int ich = 0; ich < nchannel; ++ich) {
      frame_channel[iframe * nchannel + ich] = rig_bvh::myStod(aToken[ich]);
    }
  }
  // ---------------
  for (unsigned int ibone = 0; ibone < bones.size(); ++ibone) {
    CRigBone &bone = bones[ibone];
    bone.scale = 1.0;
    bone.quatRelativeRot[0] = 0.0;
    bone.quatRelativeRot[1] = 0.0;
    bone.quatRelativeRot[2] = 0.0;
    bone.quatRelativeRot[3] = 1.0;
    bone.transRelative[0] = 0.0;
    bone.transRelative[1] = 0.0;
    bone.transRelative[2] = 0.0;
    if (bone.ibone_parent != -1) {
      const CRigBone &bone_p = bones[bone.ibone_parent];
      bone.transRelative[0] = (-bone.invBindMat[3]) - (-bone_p.invBindMat[3]);
      bone.transRelative[1] = (-bone.invBindMat[7]) - (-bone_p.invBindMat[7]);
      bone.transRelative[2] = (-bone.invBindMat[11]) - (-bone_p.invBindMat[11]);
    }
  }
  for (auto &bone: bones) {
    for (int i = 0; i < 16; ++i) { bone.affmat3Global[i] = bone.invBindMat[i]; }
    int info;
    rig_bvh::CalcInvMat(bone.affmat3Global, 4, info);
  }
}

DFM2_INLINE void delfem2::SetPose_BioVisionHierarchy(
  std::vector<CRigBone> &bones,
  const std::vector<CChannel_BioVisionHierarchy> &channels,
  const double *values) {
  for (auto &bone: bones) {
    bone.quatRelativeRot[0] = 0.0;
    bone.quatRelativeRot[1] = 0.0;
    bone.quatRelativeRot[2] = 0.0;
    bone.quatRelativeRot[3] = 1.0;
  }
  const size_t nch = channels.size();
  for (unsigned int ich = 0; ich < nch; ++ich) {
    const int ibone = channels[ich].ibone;
    const int iaxis = channels[ich].iaxis;
    const bool isrot = channels[ich].isrot;
    const double val = values[ich];
    assert(ibone < (int) bones.size());
    assert(iaxis >= 0 && iaxis < 3);
    if (!isrot) {
      bones[ibone].transRelative[iaxis] = val;
    } else {
      const double ar = val * M_PI / 180.0;
      double v0[3] = {0, 0, 0};
      v0[iaxis] = 1.0;
      double dq[4] = {v0[0] * sin(ar * 0.5), v0[1] * sin(ar * 0.5), v0[2] * sin(ar * 0.5), cos(ar * 0.5)};
      double qtmp[4];
      QuatQuat(qtmp,
               bones[ibone].quatRelativeRot, dq);
      Copy_Quat(bones[ibone].quatRelativeRot, qtmp);
    }
  }
  UpdateBoneRotTrans(bones);
}