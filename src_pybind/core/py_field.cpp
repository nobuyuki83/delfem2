/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdio>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/mshtopoio.h"
#include "delfem2/mshio.h"

namespace py = pybind11;
namespace dfm2 = delfem2;

/*
{
  std::string fpath = path_dir + std::string("/tetvelo_") + std::to_string(istep) + ".vtk";
  std::ofstream fout(fpath);
  WriteVTK_Mesh(fout, "Cd_"+std::to_string(Cd),
                aXYZ, FEMELEM_TET, aTet);
  WriteVTK_PointVec(fout, aXYZ.size()/3, aVal, 4);
}
{
  std::string fpath = path_dir + std::string("/tripress_") + std::to_string(istep) + ".vtk";
  WriteVTK_MapTriScalar(fpath, "Cd_"+std::to_string(CdP),
                        aXYZ,mapPsrf2Ptet, aTriSrf, aVal,4,3);
}
 */


void PyWrite_VTK_MeshPoint
(const std::string& file_path,
 const std::string& message,
 const py::array_t<double>& aXYZ)
{
  std::ofstream fout(file_path);
  delfem2::WriteVTK_Points(fout, message,
                           aXYZ.data(), aXYZ.shape()[0], aXYZ.shape()[1]);
}


void PyWrite_VTK_MeshElem
(const std::string& file_path,
 const py::array_t<int>& aElem,
 dfm2::MESHELEM_TYPE meshelem_type)
{
  std::ofstream fout(file_path, std::ios_base::app);
  int vtk_elem_type = 0;
  if( meshelem_type == dfm2::MESHELEM_TRI ){  vtk_elem_type = 5;  }
  if( meshelem_type == dfm2::MESHELEM_TET ){  vtk_elem_type = 10;  }
  delfem2::WriteVTK_Cells(fout, vtk_elem_type, aElem.data(), aElem.shape()[0]);
}

void PyWrite_VTK_PointScalar
(const std::string& file_path,
 const py::array_t<double>& aVal)
{
  std::ofstream fout(file_path, std::ios_base::app);
  delfem2::WriteVTK_Data_PointScalar(fout,
                            aVal.data(), aVal.shape()[0]);
}

void PyWrite_VTK_PointVector
(const std::string& file_path,
 const py::array_t<double>& aVal)
{
  std::ofstream fout(file_path, std::ios_base::app);
  delfem2::WriteVTK_Data_PointVec(fout,
                         aVal.data(),
                         aVal.shape()[0], aVal.shape()[1], aVal.shape()[1]);
}

void init_field(py::module &m){
  m.def("write_vtk_meshelem",   &PyWrite_VTK_MeshElem);
  m.def("write_vtk_meshpoint",  &PyWrite_VTK_MeshPoint);
  m.def("write_vtk_pointscalar",&PyWrite_VTK_PointScalar);
  m.def("write_vtk_pointvector",&PyWrite_VTK_PointVector);
}
