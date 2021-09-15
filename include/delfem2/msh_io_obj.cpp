/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/msh_io_obj.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <string>
#include <cstring>
#include <climits>
#include <map>

#include "delfem2/points.h"

namespace delfem2::msh_ioobj {

DFM2_INLINE void ParseVtxObj_(
    int &ip,
    int &it,
    int &in,
    char *str) {
  ip = -1;
  it = -1;
  in = -1;
  std::size_t n = strlen(str);
  int icnt = 0;
  unsigned int aloc[3];
  for (unsigned int i = 0; i < n; ++i) {
    if (str[i] != '/') { continue; }
    str[i] = '\0';
    aloc[icnt] = i;
    icnt++;
  }
  ip = std::stoi(str);
  ip--;
  if (icnt == 0) {
    return;
  }
  if (icnt == 1) {
    it = std::stoi(str + aloc[0] + 1);
    it--;
    return;
  }
  if (icnt == 2) {
    if (aloc[1] - aloc[0] != 1) {
      it = std::stoi(str + aloc[0] + 1);
      it--;
    }
    in = std::stoi(str + aloc[1] + 1);
    in--;
    return;
  }
}

}

DFM2_INLINE void delfem2::Write_Obj(
    const std::string &str,
    const double *aXYZ, int nXYZ,
    const unsigned int *aTri, int nTri) {
  int np = nXYZ;
  int nt = nTri;
  std::ofstream fout(str.c_str(), std::ofstream::out);
  for (int ip = 0; ip < np; ip++) {
    fout << "v " << aXYZ[ip * 3 + 0] << " " << aXYZ[ip * 3 + 1] << " " << aXYZ[ip * 3 + 2] << std::endl;
  }
  for (int itri = 0; itri < nt; itri++) {
    fout << "f " << aTri[itri * 3 + 0] + 1 << " " << aTri[itri * 3 + 1] + 1 << " " << aTri[itri * 3 + 2] + 1
         << std::endl;
  }
}

DFM2_INLINE void delfem2::Write_Obj_Quad(
    const std::string &str,
    const std::vector<double> &vec_xyz,
    const std::vector<int> &vec_quad) {
  int np = (int) vec_xyz.size() / 3;
  int nq = (int) vec_quad.size() / 4;
  std::ofstream fout(str.c_str(), std::ofstream::out);
  for (int ip = 0; ip < np; ip++) {
    fout << "v " << vec_xyz[ip * 3 + 0] << " " << vec_xyz[ip * 3 + 1] << " " << vec_xyz[ip * 3 + 2] << std::endl;
  }
  for (int iq = 0; iq < nq; iq++) {
    fout << "f " << vec_quad[iq * 4 + 0] + 1 << " " << vec_quad[iq * 4 + 1] + 1 << " " << vec_quad[iq * 4 + 2] + 1 << " "
         << vec_quad[iq * 4 + 3] + 1 << std::endl;
  }
}

DFM2_INLINE void delfem2::Write_Obj_ElemJArray(
    const std::string &str,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aElemInd,
    const std::vector<int> &aElem) {
  assert(!aElemInd.empty());
  const size_t np = aXYZ.size() / 3;
  std::ofstream fout(str.c_str(), std::ofstream::out);
  for (unsigned int ip = 0; ip < np; ip++) {
    fout << "v " << aXYZ[ip * 3 + 0] << " " << aXYZ[ip * 3 + 1] << " " << aXYZ[ip * 3 + 2] << std::endl;
  }
  const size_t ne = aElemInd.size() - 1;
  for (unsigned int iie = 0; iie < ne; iie++) {
    const int ie0 = aElemInd[iie];
    const int nnoel = aElemInd[iie + 1] - ie0;
    assert(nnoel == 3 || nnoel == 4);
    if (nnoel == 3) {
      fout << "f " << aElem[ie0 + 0] + 1 << " " << aElem[ie0 + 1] + 1 << " " << aElem[ie0 + 2] + 1 << std::endl;
    } else if (nnoel == 4) {
      fout << "f " << aElem[ie0 + 0] + 1 << " " << aElem[ie0 + 1] + 1 << " " << aElem[ie0 + 2] + 1 << " "
           << aElem[ie0 + 3] + 1 << std::endl;
    }
  }
}

DFM2_INLINE void delfem2::Write_Obj_TriFlag(
    const std::string &pathf,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    std::vector<unsigned int> &aFlgTri) {
  const size_t nt = aTri.size() / 3;
//  std::cout << nt << " " << aFlgTri.size() << std::endl;
  assert(aFlgTri.size() == nt);
  unsigned int flgmax = 0;
  for (unsigned int it = 0; it < nt; ++it) {
    if (aFlgTri[it] > flgmax) { flgmax = aFlgTri[it]; }
  }
//  std::cout << flgmax << std::endl;
  std::ofstream fout(pathf.c_str(), std::ofstream::out);
  const size_t np = aXYZ.size() / 3;
  for (unsigned int ip = 0; ip < np; ip++) {
    fout << "v " << aXYZ[ip * 3 + 0] << " " << aXYZ[ip * 3 + 1] << " " << aXYZ[ip * 3 + 2] << std::endl;
  }
  for (unsigned int iflg = 0; iflg < flgmax + 1; ++iflg) {
    fout << "g flag" << std::to_string(iflg) << std::endl;
    for (unsigned int it = 0; it < aTri.size() / 3; ++it) {
      if (aFlgTri[it] != iflg) { continue; }
      fout << "f " << aTri[it * 3 + 0] + 1 << " " << aTri[it * 3 + 1] + 1 << " " << aTri[it * 3 + 2] + 1 << std::endl;
    }
  }
}

DFM2_INLINE void delfem2::Write_Obj(
    const std::string &str,
    const std::vector<std::pair<std::vector<double>, std::vector<unsigned int> > > &aMesh) {
  std::ofstream fout(str.c_str(), std::ofstream::out);
  int ipsum = 0;
  for (int im = 0; im < (int) aMesh.size(); im++) {
    const std::vector<double> &aXYZ = aMesh[im].first;
    const std::vector<unsigned int> &aTri = aMesh[im].second;
    int np = (int) aXYZ.size() / 3;
    int nt = (int) aTri.size() / 3;
    { // group id
      fout << "g " << im << std::endl;
    }
    for (int ip = 0; ip < np; ip++) {
      fout << "v " << aXYZ[ip * 3 + 0] << " " << aXYZ[ip * 3 + 1] << " " << aXYZ[ip * 3 + 2] << std::endl;
    }
    for (int itri = 0; itri < nt; itri++) {
      fout << "f " << aTri[itri * 3 + 0] + 1 + ipsum << " " << aTri[itri * 3 + 1] + 1 + ipsum << " "
           << aTri[itri * 3 + 2] + 1 + ipsum << std::endl;
    }
    ipsum += np;
  }
}

DFM2_INLINE void delfem2::Read_Obj(
    const std::string &fname,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri) {
  std::ifstream fin;
  fin.open(fname.c_str());
  if (fin.fail()) {
    std::cout << "File Read Fail" << std::endl;
    return;
  }
  aXYZ.clear();
  aTri.clear();
  aXYZ.reserve(256 * 16);
  aTri.reserve(256 * 16);
  const int BUFF_SIZE = 256;
  char buff[BUFF_SIZE];
  while (fin.getline(buff, BUFF_SIZE)) {
    if (buff[0] == '#') { continue; }
    if (buff[0] == 'v' && buff[1] == ' ') {
      char str[256];
      double x, y, z;
      {
        std::istringstream is(buff);
        is >> str >> x >> y >> z;
//        sscanf(buff, "%s %lf %lf %lf", str, &x, &y, &z);
      }
      aXYZ.push_back(x);
      aXYZ.push_back(y);
      aXYZ.push_back(z);
    }
    if (buff[0] == 'f') {
      char str[256];
      int i0, i1, i2;
      {
        std::istringstream is(buff);
        is >> str >> i0 >> i1 >> i2;
//       sscanf(buff, "%s %d %d %d", str, &i0, &i1, &i2);
      }
      aTri.push_back(i0 - 1);
      aTri.push_back(i1 - 1);
      aTri.push_back(i2 - 1);
    }
  }
}

DFM2_INLINE void delfem2::Read_Obj_MeshQuad3(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aQuad,
    const std::filesystem::path &fname) {
  std::ifstream fin;
  fin.open(fname.c_str());
  if (fin.fail()) {
    std::cout << "File Read Fail" << std::endl;
    return;
  }
  aXYZ.clear();
  aQuad.clear();
  aXYZ.reserve(256 * 16);
  aQuad.reserve(256 * 16);
  const int BUFF_SIZE = 256;
  char buff[BUFF_SIZE];
  while (fin.getline(buff, BUFF_SIZE)) {
    if (buff[0] == '#') { continue; }
    if (buff[0] == 'v' && buff[1] == ' ') {
      char str[256];
      double x, y, z;
      std::istringstream is(buff);
      is >> str >> x >> y >> z;
//      sscanf(buff, "%s %lf %lf %lf", str, &x, &y, &z);
      aXYZ.push_back(x);
      aXYZ.push_back(y);
      aXYZ.push_back(z);
    }
    if (buff[0] == 'f') {
      char str[256];
      int i0, i1, i2, i3;
      std::istringstream is(buff);
      is >> str >> i0 >> i1 >> i2 >> i3;
//      sscanf(buff, "%s %d %d %d %d", str, &i0, &i1, &i2, &i3);
      aQuad.push_back(i0 - 1);
      aQuad.push_back(i1 - 1);
      aQuad.push_back(i2 - 1);
      aQuad.push_back(i3 - 1);
    }
  }
}

DFM2_INLINE void delfem2::Read_Obj2(
    const std::string &fname,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri) {
  std::ifstream fin;
  fin.open(fname.c_str());
  if (fin.fail()) {
    std::cout << "File Read Fail" << std::endl;
    return;
  }
  aXYZ.clear();
  aTri.clear();
  aXYZ.reserve(256 * 16);
  aTri.reserve(256 * 16);
  const int BUFF_SIZE = 256;
  char buff[BUFF_SIZE];
  while (fin.getline(buff, BUFF_SIZE)) {
    if (buff[0] == '#') { continue; }
    if (buff[0] == 'v' && buff[1] == ' ') {
      char str[256];
      double x, y, z;
      std::istringstream is(buff);
      is >> str >> x >> y >> z;
//      sscanf(buff, "%s %lf %lf %lf", str, &x, &y, &z);
      aXYZ.push_back(x);
      aXYZ.push_back(y);
      aXYZ.push_back(z);
    }
    if (buff[0] == 'f') {
      char str[256], str0[256], str1[256], str2[256];
      {
        std::istringstream is(buff);
        is >> str >> str0 >> str1 >> str2;
//        sscanf(buff, "%s %s %s %s", str, str0, str1, str2);
      }
      for (unsigned int i = 0; i < strlen(str0); ++i) { if (str0[i] == '/') { str0[i] = '\0'; }}
      for (unsigned int i = 0; i < strlen(str1); ++i) { if (str1[i] == '/') { str1[i] = '\0'; }}
      for (unsigned int i = 0; i < strlen(str2); ++i) { if (str2[i] == '/') { str2[i] = '\0'; }}
      const int i0 = std::stoi(str0);
      const int i1 = std::stoi(str1);
      const int i2 = std::stoi(str2);
//      sscanf(str0,"%d",&i0);
//      sscanf(str1,"%d",&i1);
//      sscanf(str2,"%d",&i2);
      aTri.push_back(i0 - 1);
      aTri.push_back(i1 - 1);
      aTri.push_back(i2 - 1);
    }
  }
}

DFM2_INLINE void delfem2::Read_Obj3(
    const std::string &fname,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri) {
  std::ifstream fin;
  fin.open(fname.c_str());
  if (fin.fail()) {
    std::cout << "File Read Fail" << std::endl;
    return;
  }
  aXYZ.clear();
  aTri.clear();
  aXYZ.reserve(256 * 16);
  aTri.reserve(256 * 16);
  const int BUFF_SIZE = 256;
  char buff[BUFF_SIZE];
  while (fin.getline(buff, BUFF_SIZE)) {
    if (buff[0] == '#') { continue; }
    if (buff[0] == 'v' && buff[1] == ' ') {
      char str[256];
      double x, y, z;
      std::istringstream is(buff);
      is >> str >> x >> y >> z;
//      sscanf(buff, "%s %lf %lf %lf", str, &x, &y, &z);
      aXYZ.push_back(x);
      aXYZ.push_back(y);
      aXYZ.push_back(z);
    }
    if (buff[0] == 'f') {
      int nv = 0;
      char *aPtr[4] = {nullptr, nullptr, nullptr, nullptr};
      for (int i = 0; i < BUFF_SIZE; ++i) {
        if (buff[i] == '\0') { break; }
        if (buff[i] == ' ') {
          aPtr[nv] = buff + i + 1;
          nv++;
        }
      }
      int aI[4] = {-1, -1, -1, -1};
      for (int iv = 0; iv < nv; ++iv) {
        for (int i = 0; i < BUFF_SIZE; ++i) {
          if (aPtr[iv][i] == '/') {
            aPtr[iv][i] = '\0';
            break;
          }
          if (aPtr[iv][i] == ' ') {
            break;
          }
        }
        const int i0 = std::stoi(aPtr[iv]);
//        std::istringstream is(aPtr[iv]);
//        is >> i0;
//        sscanf(aPtr[iv],"%d",&i0);
        aI[iv] = i0 - 1;
      }
      if (nv == 3) {
        aTri.push_back(aI[0]);
        aTri.push_back(aI[1]);
        aTri.push_back(aI[2]);
      }
      if (nv == 4) {
        aTri.push_back(aI[0]);
        aTri.push_back(aI[1]);
        aTri.push_back(aI[2]);
        ///
        aTri.push_back(aI[0]);
        aTri.push_back(aI[2]);
        aTri.push_back(aI[3]);
      }
    }
  }
}

DFM2_INLINE void delfem2::Read_WavefrontObjWithSurfaceAttributes(
    const std::string &fname,
    std::string &fname_mtl,
    std::vector<double> &vec_xyz,
    std::vector<double> &vec_tex,
    std::vector<double> &vec_nrm,
    std::vector<TriGroupWavefrontObj> &vec_tri_group) {
  std::ifstream fin;
  fin.open(fname.c_str());
  if (fin.fail()) {
    std::cout << "File Read Fail: " << fname << std::endl;
    return;
  }
  vec_xyz.clear();
  vec_tri_group.clear();
  vec_nrm.clear();
  vec_tex.clear();
  vec_xyz.reserve(256 * 16);
  vec_tri_group.reserve(100);
  const int BUFF_SIZE = 256;
  char buff[BUFF_SIZE];
  fname_mtl.clear();
  std::string current_material_name = "";
  while (fin.getline(buff, BUFF_SIZE)) {
    if (buff[0] == '#') { continue; }
    if (buff[0] == 'm') {
      std::stringstream ss(buff);
      std::string str0, str1;
      ss >> str0 >> str1;
      fname_mtl = str1;
      continue;
    }
    if (buff[0] == 'v') {
      char str[256];
      double x, y, z;
      std::istringstream is(buff);
      if (buff[1] == ' ') { // vertex
        is >> str >> x >> y >> z;
        vec_xyz.push_back(x);
        vec_xyz.push_back(y);
        vec_xyz.push_back(z);
      } else if (buff[1] == 'n') { // noraml
        is >> str >> x >> y >> z;
        double len = sqrt(x * x + y * y + z * z);
        x /= len;
        y /= len;
        z /= len;
        vec_nrm.push_back(x);
        vec_nrm.push_back(y);
        vec_nrm.push_back(z);
      } else if(buff[1] == 't'){ // tex
        is >> str >> x >> y;
        vec_tex.push_back(x);
        vec_tex.push_back(y);
      }
    }
    if (buff[0] == 'g') { // group
      vec_tri_group.resize(vec_tri_group.size() + 1);
      std::stringstream ss(buff);
      std::string str0, str1;
      ss >> str0 >> str1;
      const std::size_t iogt0 = vec_tri_group.size() - 1;
      vec_tri_group[iogt0].name_group = str1;
      vec_tri_group[iogt0].name_mtl = current_material_name;
      continue;
    }
    if (buff[0] == 'u') { // usemtl
      std::stringstream ss(buff);
      std::string str0, str1;
      ss >> str0 >> str1;
      current_material_name = str1;
      continue;
    }
    if (buff[0] == 'f') {
      if( vec_tri_group.empty() ){
        vec_tri_group.resize(vec_tri_group.size() + 1);
        const std::size_t iogt0 = vec_tri_group.size() - 1;
        vec_tri_group[iogt0].name_group = "";
        vec_tri_group[iogt0].name_mtl = current_material_name;
      }
      const std::size_t iogt0 = vec_tri_group.size() - 1;
      char str[256], str0[256], str1[256], str2[256];
      {
        std::istringstream is(buff);
        is >> str >> str0 >> str1 >> str2;
      }
      int ip0, it0, in0;
      msh_ioobj::ParseVtxObj_(ip0, it0, in0, str0);
      int ip1, it1, in1;
      msh_ioobj::ParseVtxObj_(ip1, it1, in1, str1);
      int ip2, it2, in2;
      msh_ioobj::ParseVtxObj_(ip2, it2, in2, str2);
      vec_tri_group[iogt0].vec_idx_vtx.push_back(ip0);
      vec_tri_group[iogt0].vec_idx_vtx.push_back(ip1);
      vec_tri_group[iogt0].vec_idx_vtx.push_back(ip2);
      {
        vec_tri_group[iogt0].vec_idx_tex.push_back(it0);
        vec_tri_group[iogt0].vec_idx_tex.push_back(it1);
        vec_tri_group[iogt0].vec_idx_tex.push_back(it2);
      }
      {
        vec_tri_group[iogt0].vec_idx_nrm.push_back(in0);
        vec_tri_group[iogt0].vec_idx_nrm.push_back(in1);
        vec_tri_group[iogt0].vec_idx_nrm.push_back(in2);
      }
    }
  }
}

void delfem2::Read_WavefrontMaterial(
    const std::string &fname,
    std::vector<MaterialWavefrontObj> &aMtl) {
  std::ifstream fin;
  fin.open(fname.c_str());
  if (fin.fail()) {
    std::cout << "File Read Fail" << std::endl;
    return;
  }
  aMtl.clear();
  const int BUFF_SIZE = 256;
  char buff[BUFF_SIZE];
  while (fin.getline(buff, BUFF_SIZE)) {
    if (buff[0] == '#') { continue; }
    if (buff[0] == '\n') { continue; }
    std::stringstream ss(buff);
    std::string str0, str1, str2, str3, str4;
    ss >> str0;
    if (str0 == "newmtl") {
      aMtl.resize(aMtl.size() + 1);
      const int imtl0 = static_cast<int>(aMtl.size()) - 1;
      ss >> str1;
      aMtl[imtl0].name_mtl = str1;
    }
	const int imtl0 = static_cast<int>(aMtl.size()) - 1;
    if (str0 == "Kd") {
      ss >> str1 >> str2 >> str3;
      aMtl[imtl0].Kd[0] = std::stof(str1);
      aMtl[imtl0].Kd[1] = std::stof(str2);
      aMtl[imtl0].Kd[2] = std::stof(str3);
      aMtl[imtl0].Kd[3] = 1.0;
    }
    if (str0 == "Ka") {
      ss >> str1 >> str2 >> str3;
      aMtl[imtl0].Ka[0] = std::stof(str1);
      aMtl[imtl0].Ka[1] = std::stof(str2);
      aMtl[imtl0].Ka[2] = std::stof(str3);
      aMtl[imtl0].Ka[3] = 1.0;
    }
    if (str0 == "Ks") {
      ss >> str1 >> str2 >> str3;
      aMtl[imtl0].Ks[0] = std::stof(str1);
      aMtl[imtl0].Ks[1] = std::stof(str2);
      aMtl[imtl0].Ks[2] = std::stof(str3);
      aMtl[imtl0].Ks[3] = 1.0;
    }
    if (str0 == "Ke") {
      ss >> str1 >> str2 >> str3;
      aMtl[imtl0].Ke[0] = std::stof(str1);
      aMtl[imtl0].Ke[1] = std::stof(str2);
      aMtl[imtl0].Ke[2] = std::stof(str3);
      aMtl[imtl0].Ke[3] = std::stof(str3);
    }
    if (str0 == "Ns") {
      ss >> str1;
      aMtl[imtl0].Ns = std::stof(str1);
    }
    if (str0 == "illum") {
      ss >> str1;
      aMtl[imtl0].illum = std::stoi(str1);
    }
    if (str0 == "map_Kd") {      
      ss >> str1;
      aMtl[imtl0].map_Kd = str1;
    }
  }
}

// ----------------------

void delfem2::Shape3_WavefrontObj::ReadObj(
                                      const std::string &path_obj) {
  std::string fname_mtl;
  Read_WavefrontObjWithSurfaceAttributes(
           path_obj,
           fname_mtl, aXYZ, aTex, aNorm, aObjGroupTri);
  std::string path_dir = std::string(path_obj.begin(), path_obj.begin() + path_obj.rfind("/"));
  Read_WavefrontMaterial(path_dir + "/" + fname_mtl,
           aMaterial);
  //  std::cout << aObjGroupTri.size() << " " << aMaterial.size() << std::endl;
  { //
    std::map<std::string, int> mapMtlName2Ind;
    for (int imtl = 0; imtl < (int) aMaterial.size(); ++imtl) {
      mapMtlName2Ind.insert(std::make_pair(aMaterial[imtl].name_mtl, imtl));
    }
    for (auto &iogt : aObjGroupTri) {
      std::string name_mtl = iogt.name_mtl;
      auto itr = mapMtlName2Ind.find(name_mtl);
      if (name_mtl.empty() || itr == mapMtlName2Ind.end()) {
        iogt.idx_material = -1;
        continue;
      }
      iogt.idx_material = itr->second;
    }
  }
}

std::vector<double> delfem2::Shape3_WavefrontObj::AABB3_MinMax() const {
  double c[3], w[3];
  delfem2::CenterWidth_Points3(c, w,
                               aXYZ);
  std::vector<double> aabb(6);
  aabb[0] = c[0] - 0.5 * w[0];
  aabb[1] = c[0] + 0.5 * w[0];
  aabb[2] = c[1] - 0.5 * w[1];
  aabb[3] = c[1] + 0.5 * w[1];
  aabb[4] = c[2] - 0.5 * w[2];
  aabb[5] = c[2] + 0.5 * w[2];
  return aabb;
}

void delfem2::Shape3_WavefrontObj::ScaleXYZ(
                                       double s) {
  delfem2::Scale_PointsX(aXYZ,
                         s);
}

void delfem2::Shape3_WavefrontObj::TranslateXYZ(
                                           double x, double y, double z) {
  delfem2::Translate_Points3(aXYZ,
                             x, y, z);
}
