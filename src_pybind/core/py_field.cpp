#include <stdio.h>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/funcs_gl.h"
#include "delfem2/mshtopoio_gl.h"
#include "delfem2/mshio.h"

namespace py = pybind11;


void DrawField_ColorMap
(const py::array_t<double>& pos,
 const py::array_t<unsigned int>& elm,
 const py::array_t<double>& val,
 const CColorMap& color_map)
{
//  DrawMeshTri2D_ScalarP1(me.aPos,me.aElem,a.data(),1,0,colorMap);
  const int np = pos.shape()[0];
  const int ndim = pos.shape()[1];
  const int nelm = elm.shape()[0];
  assert( val.shape()[0] == np);
  const int nstride = val.strides()[0] / sizeof(double);
  if( elm.shape()[1] == 3 ){
    if( ndim == 3 ){
      DrawMeshTri3D_ScalarP1(pos.data(), np,
                             elm.data(), nelm,
                             val.data(),
                             color_map.aColor);
    }
    else if( ndim == 2 ){
      DrawMeshTri2D_ScalarP1(pos.data(), np,
                             elm.data(), nelm,
                             val.data(), nstride,
                             color_map.aColor);
    }
  }
  if( elm.shape()[1] == 4 ){
    if( ndim == 3 ){
      DrawMeshTet3D_ScalarP1(pos.data(), np,
                             elm.data(), nelm,
                             val.data(),
                             color_map.aColor);
    }
  }
}

void DrawField_Disp
(const py::array_t<double>& pos,
 const py::array_t<unsigned int>& elm,
 MESHELEM_TYPE meshelem_type,
 const py::array_t<double>& disp)
{
  //  DrawMeshTri2D_ScalarP1(me.aPos,me.aElem,a.data(),1,0,colorMap);
  const int np = pos.shape()[0];
  const int ndim = pos.shape()[1];
  const int nelm = elm.shape()[0];
  assert( disp.shape()[0] == np );
  assert( disp.shape()[1] == ndim );
  const int nstride = disp.strides()[0] / sizeof(double);
  if( ndim == 3 ){
    if( meshelem_type == MESHELEM_TET ){
      DrawMeshTet3D_FaceNormDisp(pos.data(), np,
                                 elm.data(), nelm,
                                 disp.data());
    }
  }
  else if( ndim == 2 ){
    if( meshelem_type == MESHELEM_TRI ){
      DrawMeshTri2D_FaceDisp2D(pos.data(), np,
                               elm.data(), nelm,
                               disp.data(), nstride);
    }
  }
}

void DrawField_Hedgehog
(const py::array_t<double>& pos,
 const py::array_t<double>& disp,
 double mag)
{
  //  DrawMeshTri2D_ScalarP1(me.aPos,me.aElem,a.data(),1,0,colorMap);
  const int np = pos.shape()[0];
  const int ndim = pos.shape()[1];
  assert( disp.shape()[0] == np );
  assert( disp.shape()[1] == ndim );
  const int nstride = disp.strides()[0] / sizeof(double);
  if( ndim == 3 ){
  }
  else if( ndim == 2 ){
    DrawPoints2D_Vectors(pos.data(), np,
                         disp.data(), nstride, 0, mag);
  }
}


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
  WriteVTK_Points(fout, message,
                  aXYZ.data(), aXYZ.shape()[0], aXYZ.shape()[1]);
}


void PyWrite_VTK_MeshElem
(const std::string& file_path,
 const py::array_t<int>& aElem,
 MESHELEM_TYPE meshelem_type)
{
  std::ofstream fout(file_path, std::ios_base::app);
  int vtk_elem_type = 0;
  if( meshelem_type == MESHELEM_TRI ){  vtk_elem_type = 5;  }
  if( meshelem_type == MESHELEM_TET ){  vtk_elem_type = 10;  }
  WriteVTK_Cells(fout, vtk_elem_type, aElem.data(), aElem.shape()[0]);
}

void PyWrite_VTK_PointScalar
(const std::string& file_path,
 const py::array_t<double>& aVal)
{
  std::ofstream fout(file_path, std::ios_base::app);
  WriteVTK_Data_PointScalar(fout,
                            aVal.data(), aVal.shape()[0]);
}

void PyWrite_VTK_PointVector
(const std::string& file_path,
 const py::array_t<double>& aVal)
{
  std::ofstream fout(file_path, std::ios_base::app);
  WriteVTK_Data_PointVec(fout,
                         aVal.data(),
                         aVal.shape()[0], aVal.shape()[1], aVal.shape()[1]);
}

void init_field(py::module &m){
  m.def("drawField_colorMap",   &DrawField_ColorMap);
  m.def("drawField_disp",       &DrawField_Disp);
  m.def("drawField_hedgehog",   &DrawField_Hedgehog);
  m.def("write_vtk_meshelem",   &PyWrite_VTK_MeshElem);
  m.def("write_vtk_meshpoint",  &PyWrite_VTK_MeshPoint);
  m.def("write_vtk_pointscalar",&PyWrite_VTK_PointScalar);
  m.def("write_vtk_pointvector",&PyWrite_VTK_PointVector);
}
