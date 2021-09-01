/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <string>
#include <cstring>
#include <climits>

#include "delfem2/msh_ioobj.h"

namespace delfem2::msh_ioobj {

DFM2_INLINE void UnitNormalAreaTri3D_(
    double *n, double &a,
    const double *v1, const double *v2, const double *v3) {
  n[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
  n[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
  n[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
  a = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) * 0.5;
  const double invlen = 0.5 / a;
  n[0] *= invlen;
  n[1] *= invlen;
  n[2] *= invlen;
}

// probably std::stroi is safer to use but it is only for C++11
DFM2_INLINE int MyStoi_(const std::string &str) {
  char *e;
  long d = std::strtol(str.c_str(), &e, 0);
  return (int) d;
}

DFM2_INLINE double MyStod_(const std::string &str) {
  char *e;
  double d = std::strtod(str.c_str(), &e);
  return d;
}

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
  int aloc[3];
  for (unsigned int i = 0; i < n; ++i) {
    if (str[i] != '/') { continue; }
    str[i] = '\0';
    aloc[icnt] = i;
    icnt++;
  }
  ip = std::stoi(str);
//  sscanf(str,"%d",&ip);
  ip--;
  if (icnt == 0) {
    return;
  }
  if (icnt == 1) {
    it = std::stoi(str + aloc[0] + 1);
//    sscanf(str+aloc[0]+1,"%d",&it);
    it--;
    return;
  }
  if (icnt == 2) {
    if (aloc[1] - aloc[0] != 1) {
      it = std::stoi(str + aloc[0] + 1);
//      sscanf(str+aloc[0]+1,"%d",&it);
      it--;
    }
    in = std::stoi(str + aloc[1] + 1);
//    sscanf(str+aloc[1]+1,"%d",&in);
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
    const std::vector<double> &aXYZ,
    const std::vector<int> &aQuad) {
  int np = (int) aXYZ.size() / 3;
  int nq = (int) aQuad.size() / 4;
  std::ofstream fout(str.c_str(), std::ofstream::out);
  for (int ip = 0; ip < np; ip++) {
    fout << "v " << aXYZ[ip * 3 + 0] << " " << aXYZ[ip * 3 + 1] << " " << aXYZ[ip * 3 + 2] << std::endl;
  }
  for (int iq = 0; iq < nq; iq++) {
    fout << "f " << aQuad[iq * 4 + 0] + 1 << " " << aQuad[iq * 4 + 1] + 1 << " " << aQuad[iq * 4 + 2] + 1 << " "
         << aQuad[iq * 4 + 3] + 1 << std::endl;
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
    const std::string &fname) {
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
      char *aPtr[4] = {0, 0, 0, 0};
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

DFM2_INLINE void delfem2::Load_Obj(
    const std::string &fname,
    std::string &fname_mtl,
    std::vector<double> &aXYZ,
    std::vector<double> &aNorm,
    std::vector<CTriGroup> &aTriGroup) {
  std::ifstream fin;
  fin.open(fname.c_str());
  if (fin.fail()) {
    std::cout << "File Read Fail: " << fname << std::endl;
    return;
  }
  aXYZ.clear();
  aTriGroup.clear();
  aNorm.clear();
  aXYZ.reserve(256 * 16);
  aTriGroup.reserve(100);
  const int BUFF_SIZE = 256;
  char buff[BUFF_SIZE];
  fname_mtl.clear();
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
      is >> str >> x >> y >> z;
//      sscanf(buff, "%s %lf %lf %lf", str, &x, &y, &z);
      if (buff[1] == ' ') { // vertex
        aXYZ.push_back(x);
        aXYZ.push_back(y);
        aXYZ.push_back(z);
      } else if (buff[1] == 'n') { // noraml
        double len = sqrt(x * x + y * y + z * z);
        x /= len;
        y /= len;
        z /= len;
        aNorm.push_back(x);
        aNorm.push_back(y);
        aNorm.push_back(z);
      }
    }
    if (buff[0] == 'g') { // group
      aTriGroup.resize(aTriGroup.size() + 1);
      std::stringstream ss(buff);
      std::string str0, str1;
      ss >> str0 >> str1;
      const std::size_t iogt0 = aTriGroup.size() - 1;
      aTriGroup[iogt0].name_group = str1;
      continue;
    }
    if (buff[0] == 'u') { // usemtl
      std::stringstream ss(buff);
      std::string str0, str1;
      ss >> str0 >> str1;
      const std::size_t iogt0 = aTriGroup.size() - 1;
      aTriGroup[iogt0].name_mtl = str1;
      continue;
    }
    if (buff[0] == 'f') {
      if (aTriGroup.empty()) {
        std::cout << fname << std::endl;
      }
      const std::size_t iogt0 = aTriGroup.size() - 1;
      char str[256], str0[256], str1[256], str2[256];
      {
//      sscanf(buff, "%s %s %s %s", str, str0, str1, str2);
        std::istringstream is(buff);
        is >> str >> str0 >> str1 >> str2;
        //      std::cout << str0 << " " << str1 << " " << str2 << std::endl;
      }
      int ip0, it0, in0;
      msh_ioobj::ParseVtxObj_(ip0, it0, in0, str0);
      int ip1, it1, in1;
      msh_ioobj::ParseVtxObj_(ip1, it1, in1, str1);
      int ip2, it2, in2;
      msh_ioobj::ParseVtxObj_(ip2, it2, in2, str2);
      aTriGroup[iogt0].aTriVtx.push_back(ip0);
      aTriGroup[iogt0].aTriVtx.push_back(ip1);
      aTriGroup[iogt0].aTriVtx.push_back(ip2);
//     if( in0 != -1 && in1 != -1 && in2 != -1 ){
      {
        aTriGroup[iogt0].aTriNrm.push_back(in0);
        aTriGroup[iogt0].aTriNrm.push_back(in1);
        aTriGroup[iogt0].aTriNrm.push_back(in2);
      }
    }
  }
}
