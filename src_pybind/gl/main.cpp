#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vector>
#include <map>
#include <deque>

#include "delfem2/dyntri_v2.h"
#include "delfem2/dyntri_v3.h"
#include "delfem2/mshtopo.h"

#include "delfem2/glew_funcs.h" // have to be included in the beginning
#include "delfem2/gl_funcs.h"
#include "delfem2/gl_cad_dyntri_v23.h"
#include "delfem2/gl_voxbv.h"
#include "delfem2/gl_color.h"

namespace py = pybind11;

//////////////////////////////////////////////////////////////////////////////////////////

void init_sampler(py::module &m);
void init_texture(py::module &m);


void PyDrawEdge_CMeshDynTri3D(const CMeshDynTri3D& dmsh){
  DrawMeshDynTri_Edge(dmsh.aETri,dmsh.aVec3);
}

void PyDrawEdge_CMeshDynTri2D(const CMeshDynTri2D& dmsh){
  DrawMeshDynTri_Edge(dmsh.aETri,dmsh.aVec2);
}

void PyDraw_CCad2D(const CCad2D& cad){
  Draw_CCad2D(cad);
}

void PyDrawMesh_FaceNorm
(const py::array_t<double>& pos,
 const py::array_t<unsigned int>& elm,
 const MESHELEM_TYPE type)
{
  assert(pos.ndim()==2);
  assert(elm.ndim()==2);
  const auto shape_pos = pos.shape();
  const auto shape_elm = elm.shape();
  if( shape_pos[1] == 3 ){ // 3D Mesh
    if( type == MESHELEM_TRI  ){  DrawMeshTri3D_FaceNorm( pos.data(), elm.data(), shape_elm[0]); }
    if( type == MESHELEM_QUAD ){  DrawMeshQuad3D_FaceNorm(pos.data(), elm.data(), shape_elm[0]); }
    if( type == MESHELEM_HEX  ){  DrawMeshHex3D_FaceNorm( pos.data(), elm.data(), shape_elm[0]); }
    if( type == MESHELEM_TET  ){  DrawMeshTet3D_FaceNorm( pos.data(), elm.data(), shape_elm[0]); }
  }
}

void PyDrawMesh_Edge
(const py::array_t<double>& pos,
 const py::array_t<unsigned int>& elm,
 const MESHELEM_TYPE type)
{
  assert(pos.ndim()==2);
  assert(elm.ndim()==2);
  const auto shape_pos = pos.shape();
  const auto shape_elm = elm.shape();
  if( shape_pos[1] == 3 ){ // 3D Mesh
    if( type == MESHELEM_TRI  ){ DrawMeshTri3D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == MESHELEM_QUAD ){ DrawMeshQuad3D_Edge(pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == MESHELEM_HEX  ){ DrawMeshHex3D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == MESHELEM_TET  ){ DrawMeshTet3D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == MESHELEM_LINE ){ DrawMeshLine3D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
  }
  if( shape_pos[1] == 2 ){ // 2D Mesh
    if( type == MESHELEM_TRI  ){  DrawMeshTri2D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == MESHELEM_QUAD ){  DrawMeshQuad2D_Edge(pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
  }
}



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


PYBIND11_MODULE(c_gl, m) {
  m.doc() = "pybind11 delfem2 binding";
  ///////////////////////////////////
  init_sampler(m);
  init_texture(m);
   
 ////////////////////////////////////

  py::class_<CColorMap>(m,"ColorMap")
  .def(py::init<>())
  .def(py::init<double, double, const std::string&>());
  
  m.def("cppDrawEdge_CppMeshDynTri2D", &PyDrawEdge_CMeshDynTri2D);
  m.def("cppDrawEdge_CppMeshDynTri3D", &PyDrawEdge_CMeshDynTri3D);
  m.def("cppDraw_CppCad2D",            &PyDraw_CCad2D);
  
  
  m.def("draw_mesh_facenorm",  &PyDrawMesh_FaceNorm);
  m.def("draw_mesh_edge",      &PyDrawMesh_Edge);
  
  m.def("drawField_colorMap",   &DrawField_ColorMap);
  m.def("drawField_disp",       &DrawField_Disp);
  m.def("drawField_hedgehog",   &DrawField_Hedgehog);
  
  ///////////////////////////////////
  // gl misc
  m.def("setSomeLighting",  &setSomeLighting, "set some lighting that looks good for me");
  m.def("setup_glsl",       &setUpGLSL, "compile shader program");
  m.def("glew_init",        &glewInit);
  
  m.def("cppDrawSphere",      &DrawSphereAt );
  m.def("cppDrawSphere_Edge", &DrawSphere_Edge);
  m.def("cppDrawTorus_Edge",  &DrawTorus_Edge);
}


