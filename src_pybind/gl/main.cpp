#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>
#include <map>
#include <deque>
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/dtri3_v3dtri.h"
#include "delfem2/mshtopo.h"

// -------------------------
#include "glad/glad.h"
#if defined(__APPLE__) && defined(__MACH__)
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include "delfem2/opengl/gl_funcs.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/cad2dtriv2_glold.h"
#include "delfem2/opengl/caddtri_v3_glold.h"
#include "delfem2/opengl/v2_glold.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/gridcube_glold.h"

namespace py = pybind11;
namespace dfm2 = delfem2;

// TODO put assert for the input python arrays

// --------------------------------------

void init_sampler(py::module &m);
void init_texture(py::module &m);
//void init_rigidasm(py::module &m);


void PyDrawEdge_CMeshDynTri3D(const dfm2::CMeshDynTri3D& dmsh){
  dfm2::opengl::DrawMeshDynTri_Edge(dmsh.aETri,dmsh.aVec3);
}

void PyDrawEdge_CMeshDynTri2D(const dfm2::CMeshDynTri2D& dmsh){
  dfm2::opengl::DrawMeshDynTri_Edge(dmsh.aETri,dmsh.aVec2);
}

void PyDraw_CCad2D(const dfm2::CCad2D& cad){
  dfm2::opengl::Draw_CCad2D(cad);
}

void PyDrawMesh_FaceNorm
(const py::array_t<double>& pos,
 const py::array_t<unsigned int>& elm,
 const dfm2::MESHELEM_TYPE type)
{
  assert(pos.ndim()==2);
  assert(elm.ndim()==2);
  const auto shape_pos = pos.shape();
  const auto shape_elm = elm.shape();
  if( shape_pos[1] == 3 ){ // 3D Mesh
    if( type == dfm2::MESHELEM_TRI  ){  dfm2::opengl::DrawMeshTri3D_FaceNorm( pos.data(), elm.data(), shape_elm[0]); }
    if( type == dfm2::MESHELEM_QUAD ){  dfm2::opengl::DrawMeshQuad3D_FaceNorm(pos.data(), elm.data(), shape_elm[0]); }
    if( type == dfm2::MESHELEM_HEX  ){  dfm2::opengl::DrawMeshHex3D_FaceNorm( pos.data(), elm.data(), shape_elm[0]); }
    if( type == dfm2::MESHELEM_TET  ){  dfm2::opengl::DrawMeshTet3D_FaceNorm( pos.data(), elm.data(), shape_elm[0]); }
  }
}

void PyDrawMesh_Edge
(const py::array_t<double>& pos,
 const py::array_t<unsigned int>& elm,
 const dfm2::MESHELEM_TYPE type)
{
  assert(pos.ndim()==2);
  assert(elm.ndim()==2);
  const auto shape_pos = pos.shape();
  const auto shape_elm = elm.shape();
  if( shape_pos[1] == 3 ){ // 3D Mesh
    if( type == dfm2::MESHELEM_TRI  ){  dfm2::opengl::DrawMeshTri3D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == dfm2::MESHELEM_QUAD ){  dfm2::opengl::DrawMeshQuad3D_Edge(pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == dfm2::MESHELEM_HEX  ){  dfm2::opengl::DrawMeshHex3D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == dfm2::MESHELEM_TET  ){  dfm2::opengl::DrawMeshTet3D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == dfm2::MESHELEM_LINE ){  dfm2::opengl::DrawMeshLine3D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
  }
  if( shape_pos[1] == 2 ){ // 2D Mesh
    if( type == dfm2::MESHELEM_TRI  ){  dfm2::opengl::DrawMeshTri2D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == dfm2::MESHELEM_QUAD ){  dfm2::opengl::DrawMeshQuad2D_Edge(pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
  }
}



void DrawField_ColorMap
(const py::array_t<double>& pos,
 const py::array_t<unsigned int>& elm,
 const py::array_t<double>& val,
 const dfm2::CColorMap& color_map)
{
  //  DrawMeshTri2D_ScalarP1(me.aPos,me.aElem,a.data(),1,0,colorMap);
  const int np = pos.shape()[0];
  const int ndim = pos.shape()[1];
  const int nelm = elm.shape()[0];
  assert( val.shape()[0] == np);
  const int nstride = val.strides()[0] / sizeof(double);
  if( elm.shape()[1] == 3 ){
    if( ndim == 3 ){
      dfm2::opengl::DrawMeshTri3D_ScalarP1(pos.data(), np,
                                     elm.data(), nelm,
                                     val.data(),
                                     color_map.aColor);
    }
    else if( ndim == 2 ){
      dfm2::opengl::DrawMeshTri2D_ScalarP1(pos.data(), np,
                                     elm.data(), nelm,
                                     val.data(), nstride,
                                     color_map.aColor);
    }
  }
  if( elm.shape()[1] == 4 ){
    if( ndim == 3 ){
      dfm2::opengl::DrawMeshTet3D_ScalarP1(pos.data(), np,
                                     elm.data(), nelm,
                                     val.data(),
                                     color_map.aColor);
    }
  }
}

void DrawField_Disp
(const py::array_t<double>& pos,
 const py::array_t<unsigned int>& elm,
 dfm2::MESHELEM_TYPE meshelem_type,
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
    if( meshelem_type == dfm2::MESHELEM_TET ){
      dfm2::opengl::DrawMeshTet3D_FaceNormDisp(pos.data(), np,
                                         elm.data(), nelm,
                                         disp.data());
    }
  }
  else if( ndim == 2 ){
    if( meshelem_type == dfm2::MESHELEM_TRI ){
      dfm2::opengl::DrawMeshTri2D_FaceDisp2D(pos.data(), np,
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
    dfm2::opengl::DrawPoints2D_Vectors(pos.data(), np,
                                 disp.data(), nstride, 0, mag);
  }
}

void PyGLExtensionInit(){
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
}


PYBIND11_MODULE(c_gl, m) {
  m.doc() = "pybind11 delfem2 binding";
  
  // ---------------------
  init_sampler(m);
  init_texture(m);
//  init_rigidasm(m);
   
 // ----------------------

  py::class_<dfm2::CColorMap>(m,"ColorMap")
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

  
  // ----------------
  // gl misc
  m.def("setSomeLighting",  &dfm2::opengl::setSomeLighting, "set some lighting that looks good for me");
  m.def("setup_glsl",       &dfm2::opengl::setUpGLSL, "compile shader program");
  m.def("glew_init",        &PyGLExtensionInit);
  
  m.def("cppDrawSphere",      &dfm2::opengl::DrawSphereAt );
  m.def("cppDrawSphere_Edge", &dfm2::opengl::DrawSphere_Edge);
  m.def("cppDrawTorus_Edge",  &dfm2::opengl::DrawTorus_Edge);
}


