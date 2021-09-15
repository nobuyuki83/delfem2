/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/msh_io_ply.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>

// ----------------------------------------------------

DFM2_INLINE void delfem2::Read_Ply(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    const std::filesystem::path &file_path) {
  std::ifstream fin;
  fin.open(file_path.c_str());
  if (fin.fail()) {
    std::cout << "Fail Read Fail" << std::endl;
    return;
  }
  const unsigned int nbuff = 256;
  char buff[nbuff], buff1[nbuff], buff2[nbuff];
  std::string str1, str2;
  fin.getline(buff, nbuff);  // ply
  fin.getline(buff, nbuff);  // format asi 1.0
  for (;;) {
    fin.getline(buff, nbuff);
    if (strncmp(buff, "comment ", 8) != 0) { break; }
  }
  // ----
  int nno;
  {
    std::istringstream is(buff);
    is >> buff1 >> buff2 >> nno;
//    sscanf(buff, "%s %s %d", buff1, buff2, &nno);
  }
  // ----
  for (;;) {
    fin.getline(buff, nbuff);
    if (strncmp(buff, "property ", 9) != 0) { break; }
  }
  int ntri;
  {
    std::istringstream is(buff);
    is >> buff1 >> buff2 >> ntri;
//    sscanf(buff, "%s %s %d", buff1, buff2, &ntri);
  }
  aTri.resize(ntri * 3);
  // ----
  fin.getline(buff, nbuff);  // property list int int vertex_indices
  fin.getline(buff, nbuff);  // end header
  // ----
  aXYZ.resize(nno * 3);
  for (int ino = 0; ino < nno; ++ino) {
    double x, y, z;
    fin >> x >> y >> z;
    //    std::cout << ino << " " << x << " " << y << " " << z << std::endl;
    aXYZ[ino * 3 + 0] = x;
    aXYZ[ino * 3 + 1] = y;
    aXYZ[ino * 3 + 2] = z;
  }
  for (int itri = 0; itri < ntri; ++itri) {
    int itmp, i1, i2, i3;
    fin >> itmp >> i1 >> i2 >> i3;
    aTri[itri * 3 + 0] = i1;
    aTri[itri * 3 + 1] = i2;
    aTri[itri * 3 + 2] = i3;
    //    std::cout << itri << " " << itmp << " " << i1 << " " << i2 << " " << i3 << std::endl;
  }
  //  if( is_norm_ ){ this->MakeNormal(); }
}

void delfem2::Write_Ply(
    const std::string &fname,
    unsigned int nXYZ,
    double *paXYZ,
    unsigned int nTri,
    unsigned int *paTri) {
  std::cout << "File load " << fname << std::endl;
  std::ofstream fout;
  fout.open(fname.c_str());
  if (fout.fail()) {
    std::cout << "Fail Read Fail" << std::endl;
    return;
  }
  fout << "ply\n";
  fout << "format ascii 1.0\n";
  fout << "element vertex " << nXYZ << "\n";
  fout << "property float x\n";
  fout << "property float y\n";
  fout << "property float z\n";
  fout << "element face " << nTri << "\n";
  fout << "property list uchar int vertex_indices" << "\n";
  fout << "end_header\n";
  //
  for (unsigned int ixyz = 0; ixyz < nXYZ; ixyz++) {
    fout << paXYZ[ixyz * 3 + 0] << " " << paXYZ[ixyz * 3 + 1] << " " << paXYZ[ixyz * 3 + 2] << "\n";
  }
  for (unsigned int itri = 0; itri < nTri; itri++) {
    fout << "3 " << paTri[itri * 3 + 0] << " " << paTri[itri * 3 + 1] << " " << paTri[itri * 3 + 2] << "\n";
  }
}

DFM2_INLINE void delfem2::Write_Ply(
    const std::string &fname,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTri) {
  std::cout << "Open file for writing: " << fname << std::endl;
  std::ofstream fout;
  fout.open(fname.c_str(), std::ios::out);
  if (fout.fail()) {
    std::cout << "File Open Fail" << std::endl;
    return;
  }
  int nXYZ = (int) aXYZ.size() / 3;
  int nTri = (int) aTri.size() / 3;
  fout << "ply\n";
  fout << "format ascii 1.0\n";
  fout << "element vertex " << nXYZ << "\n";
  fout << "property float x\n";
  fout << "property float y\n";
  fout << "property float z\n";
  fout << "element face " << nTri << "\n";
  fout << "property list uchar int vertex_indices" << "\n";
  fout << "end_header\n";
  // ----
  for (int ixyz = 0; ixyz < nXYZ; ixyz++) {
    fout << aXYZ[ixyz * 3 + 0] << " " << aXYZ[ixyz * 3 + 1] << " " << aXYZ[ixyz * 3 + 2] << "\n";
  }
  for (int itri = 0; itri < nTri; itri++) {
    fout << "3 " << aTri[itri * 3 + 0] << " " << aTri[itri * 3 + 1] << " " << aTri[itri * 3 + 2] << "\n";
  }
}

DFM2_INLINE void delfem2::Read_Ply(
    const std::string &fname,
    int &nnode_,
    double *&pXYZs_,
    int &ntri_,
    unsigned int *&aTri_) {
  std::cout << "File load " << fname << std::endl;
  std::ifstream fin;
  fin.open(fname.c_str());
  if (fin.fail()) {
    std::cout << "Fail Read Fail" << std::endl;
    return;
  }
  const unsigned int nbuff = 256;
  char buff[nbuff];
  //  char buff1[nbuff], buff2[nbuff];
  std::string sbuff1, sbuff2;
  std::string str1, str2;
  fin.getline(buff, nbuff);  // ply
  fin.getline(buff, nbuff);  // format asi 1.0
  for (;;) {
    fin.getline(buff, nbuff);
    if (strncmp(buff, "comment ", 8) != 0) { break; }
  }
  /////
  { // read number of points
    //    sscanf(buff, "%s %s %d", buff1, buff2, &nnode_);
    std::stringstream ss(buff);
    ss >> sbuff1 >> sbuff2 >> nnode_;
    std::cout << "Nnode " << nnode_ << std::endl;
  }
  ////
  for (;;) {
    fin.getline(buff, nbuff);
    if (strncmp(buff, "property ", 9) != 0) { break; }
  }
  { // read number of triangles
    //    sscanf(buff, "%s %s %d", buff1, buff2, &ntri_);
    std::stringstream ss(buff);
    ss >> sbuff1 >> sbuff2 >> ntri_;
    std::cout << "NTri " << ntri_ << std::endl;
  }
  /////
  fin.getline(buff, nbuff);  // property list int int vertex_indices
  fin.getline(buff, nbuff);  // end header
  ////
  pXYZs_ = new double[nnode_ * 3];
  for (int ino = 0; ino < (int) nnode_; ++ino) {
    double x, y, z;
    fin >> x >> y >> z;
    //    std::cout << ino << " " << x << " " << y << " " << z << std::endl;
    pXYZs_[ino * 3 + 0] = x;
    pXYZs_[ino * 3 + 1] = y;
    pXYZs_[ino * 3 + 2] = z;
  }
  aTri_ = new unsigned int[ntri_ * 3];
  for (int itri = 0; itri < (int) ntri_; ++itri) {
    int itmp, i1, i2, i3;
    fin >> itmp >> i1 >> i2 >> i3;
    aTri_[itri * 3 + 0] = i1;
    aTri_[itri * 3 + 1] = i2;
    aTri_[itri * 3 + 2] = i3;
    //    std::cout << itri << " " << itmp << " " << i1 << " " << i2 << " " << i3 << std::endl;
  }
  //  if( is_norm_ ){ this->MakeNormal(); }
}

