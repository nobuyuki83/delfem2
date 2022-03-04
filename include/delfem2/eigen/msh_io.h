//
// Created by Nobuyuki Umetani on 2021-08-14.
//

#ifndef DFM2_EIGEN_MSH_IO_H
#define DFM2_EIGEN_MSH_IO_H

#include <Eigen/Dense>

#include "delfem2/msh_io_obj.h"

std::tuple<
    Eigen::Matrix<double, -1, 3, Eigen::RowMajor>,
    Eigen::Matrix<unsigned int, -1, 3, Eigen::RowMajor> >
ReadTriangleMeshObj(
    const std::string &fpath) {
  std::vector<double> vec_xyz;
  std::vector<unsigned int> vec_tri;
  delfem2::Read_Obj3(
      vec_xyz, vec_tri,
      fpath);
  const auto V =
      Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(
          vec_xyz.data(),
          static_cast<unsigned int>(vec_xyz.size() / 3), 3);
  const auto F =
      Eigen::Map<Eigen::Matrix<unsigned int, -1, -1, Eigen::RowMajor>>(
          vec_tri.data(),
          static_cast<unsigned int>(vec_tri.size() / 3), 3);
  return std::make_tuple(V, F);
}

// -----------------------------------------------

void WriteUniformMesh(
    const std::string &file_path,
    const Eigen::Matrix<double, -1, -1, Eigen::RowMajor> &v,
    const Eigen::Matrix<unsigned int, -1, -1, Eigen::RowMajor> &f) {
  const std::filesystem::path path(file_path);
  if (path.extension() == ".obj" || path.extension() == ".OBJ") {
    delfem2::Write_Obj_UniformMesh(
        file_path,
        v.data(), v.rows(),
        f.data(), f.rows(), f.cols());
  } else {
    std::cerr << "Error->there's no support for extension: " << path.extension() << std::endl;
  }
}

#endif
