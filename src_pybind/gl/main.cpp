#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vector>
#include <map>
#include <deque>

#include "glad/glad.h"
#if defined(__APPLE__) && defined(__MACH__)
  #include <OpenGL/gl.h>
  #include <OpenGL/glu.h>
#else
  #include <GL/gl.h>
  #include <GL/glu.h>
#endif


#include "delfem2/dtri_v2.h"
#include "delfem2/dtri_v3.h"
#include "delfem2/mshtopo.h"
#include "delfem2/eigen_rigidbody.h"

#include "delfem2/opengl/gl2ew_funcs.h" // have to be included in the beginning
#include "delfem2/opengl/gl2_funcs.h"
#include "delfem2/opengl/gl24_funcs.h"
#include "delfem2/opengl/gl2_color.h"
#include "delfem2/opengl/gl2_v23dtricad.h"
#include "delfem2/opengl/gl2_v23.h"
#include "delfem2/opengl/gl_voxbv.h"

namespace py = pybind11;


// --------------------------------------


void DrawRigidBodyAssemblyStatic(const CRigidBodyAssembly_Static& rba)
{
  const double small_rad = 0.1;
  const double big_rad = 0.1;
  GLUquadricObj *quadricSphere=gluNewQuadric();
  
  ::glDisable(GL_LIGHTING);
  const std::vector<CRigidBody>& aRigidBody = rba.aRigidBody;
  for(unsigned int irb=0;irb<aRigidBody.size();irb++){
    const CRigidBody rb = aRigidBody[irb];
    CVector3 cg = rb.cg;
    if( rba.is_draw_deformed ){
      cg += rb.u;
    }
    
    if( rba.is_draw_skeleton ){
      ::glColor3d(0,1,0);
      ::glPushMatrix();
      ::glTranslated(cg.x,cg.y,cg.z);
        //::glutWireSphere(0.1, 16, 16);
      gluSphere(quadricSphere, big_rad, 16, 16);
      ::glPopMatrix();
    }
    for(unsigned int icp=0;icp<rb.aCP.size();icp++){
      CVector3 p = rb.aCP[icp];
      if( rba.is_draw_deformed ){
        p = rb.R*(p-rb.cg)+rb.cg + rb.u;
      }
      ::glColor3d(1,0,0);
      ::glPushMatrix();
      ::glTranslated(p.x,p.y,p.z);
        //::glutWireSphere(0.02, 16, 16);
      gluSphere(quadricSphere, small_rad, 16, 16);
      ::glPopMatrix();
      
      if( rba.is_draw_skeleton ){
        ::glLineWidth(5);
        ::glColor3d(0,0,0);
        ::glBegin(GL_LINES);
        opengl::myGlVertex(cg);
        opengl::myGlVertex( p);
        ::glEnd();
      }
      
      if( rba.is_draw_force ){
        if (!rb.aCForce.empty()) {
          CVector3 q = p + rba.scale_force*rb.aCForce[icp];
          ::glLineWidth(5);
          ::glColor3d(1,0,0);
          ::glBegin(GL_LINES);
          opengl::myGlVertex( p);
          opengl::myGlVertex( q);
          ::glEnd();
        }
      }
    }
  }
  /////////////////////////////////////////////
  const std::vector<CJoint>& aJoint = rba.aJoint;
  for(unsigned int ij=0;ij<aJoint.size();ij++){
    const CJoint& joint = aJoint[ij];
    int irb0 = joint.irb0;
    int irb1 = joint.irb1;
    const CRigidBody& rb0 = aRigidBody[irb0];
    const CRigidBody& rb1 = aRigidBody[irb1];
    
    CVector3 p0 = joint.p;
    CVector3 p1 = joint.p;
    if( rba.is_draw_deformed ){
      p0 = rb0.R*(p0-rb0.cg)+rb0.cg + rb0.u;
      p1 = rb1.R*(p1-rb1.cg)+rb1.cg + rb1.u;
    }
    
    
      // joint point seen from rb0
    ::glColor3d(0,0,1);
    ::glPushMatrix();
    ::glTranslated(p0.x,p0.y,p0.z);
      //::glutWireSphere(0.02, 16, 16);
    gluSphere(quadricSphere, small_rad, 16, 16);
    ::glPopMatrix();
    
      // joint point seen from rb1
    ::glPushMatrix();
    ::glTranslated(p1.x,p1.y,p1.z);
      //::glutWireSphere(0.02, 16, 16);
    gluSphere(quadricSphere, small_rad, 16, 16);
    ::glPopMatrix();
    
    CVector3 cg0 = rb0.cg;
    CVector3 cg1 = rb1.cg;
    if( rba.is_draw_deformed ){
      cg0 += rb0.u;
      cg1 += rb1.u;
    }
    if( rba.is_draw_skeleton ){
      ::glLineWidth(5);
      ::glColor3d(0,0,0);
      ::glBegin(GL_LINES);
      opengl::myGlVertex(cg0);
      opengl::myGlVertex( p0);
      opengl::myGlVertex(cg1);
      opengl::myGlVertex( p1);
      ::glEnd();
    }
    
    if( rba.is_draw_force ){
      ::glLineWidth(5);
      CVector3 q0 = p0 + rba.scale_force*joint.linear;
      CVector3 q1 = p1 - rba.scale_force*joint.linear;
      ::glColor3d(1,0,1);
      ::glBegin(GL_LINES);
      opengl::myGlVertex( p0);
      opengl::myGlVertex( q0);
      opengl::myGlVertex( p1);
      opengl::myGlVertex( q1);
      ::glEnd();
      CVector3 r0 = p0 + rba.scale_torque*joint.torque;
      CVector3 r1 = p1 - rba.scale_torque*joint.torque;
      ::glColor3d(0,1,1);
      ::glBegin(GL_LINES);
      opengl::myGlVertex( p0);
      opengl::myGlVertex( r0);
      opengl::myGlVertex( p1);
      opengl::myGlVertex( r1);
      ::glEnd();
    }
    
  }
  
  glLineWidth(1.0f);
  gluDeleteQuadric(quadricSphere);
  
}



void DrawFloorGL(const CRigidBodyAssembly_Static& rba) {
  
  if( !rba.is_draw_grid ) { return; }
  
    // draw floor
  ::glLineWidth(1);
  ::glBegin(GL_LINES);
  ::glColor3d(0,0,0);
  double grid_x_min = -10;
  double grid_x_max = +10;
  double grid_z_min = -10;
  double grid_z_max = +10;
  unsigned int ndiv_grid = 30;
  for(unsigned int ix=0;ix<ndiv_grid+1;ix++){
    double x0 = (grid_x_max-grid_x_min) / ndiv_grid * ix + grid_x_min;
    ::glVertex3d(x0,0,grid_z_min);
    ::glVertex3d(x0,0,grid_z_max);
  }
  for(unsigned int iz=0;iz<ndiv_grid+1;iz++){
    double z0 = (grid_z_max-grid_z_min) / ndiv_grid * iz + grid_z_min;
    ::glVertex3d(grid_x_min,0,z0);
    ::glVertex3d(grid_x_max,0,z0);
  }
  ::glEnd();

}


// --------------------------------------

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
    if( type == MESHELEM_TRI  ){  opengl::DrawMeshTri3D_FaceNorm( pos.data(), elm.data(), shape_elm[0]); }
    if( type == MESHELEM_QUAD ){  opengl::DrawMeshQuad3D_FaceNorm(pos.data(), elm.data(), shape_elm[0]); }
    if( type == MESHELEM_HEX  ){  opengl::DrawMeshHex3D_FaceNorm( pos.data(), elm.data(), shape_elm[0]); }
    if( type == MESHELEM_TET  ){  opengl::DrawMeshTet3D_FaceNorm( pos.data(), elm.data(), shape_elm[0]); }
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
    if( type == MESHELEM_TRI  ){  opengl::DrawMeshTri3D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == MESHELEM_QUAD ){  opengl::DrawMeshQuad3D_Edge(pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == MESHELEM_HEX  ){  opengl::DrawMeshHex3D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == MESHELEM_TET  ){  opengl::DrawMeshTet3D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == MESHELEM_LINE ){  opengl::DrawMeshLine3D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
  }
  if( shape_pos[1] == 2 ){ // 2D Mesh
    if( type == MESHELEM_TRI  ){  opengl::DrawMeshTri2D_Edge( pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
    if( type == MESHELEM_QUAD ){  opengl::DrawMeshQuad2D_Edge(pos.data(), shape_pos[0], elm.data(), shape_elm[0]); }
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
      opengl::DrawMeshTri3D_ScalarP1(pos.data(), np,
                                     elm.data(), nelm,
                                     val.data(),
                                     color_map.aColor);
    }
    else if( ndim == 2 ){
      opengl::DrawMeshTri2D_ScalarP1(pos.data(), np,
                                     elm.data(), nelm,
                                     val.data(), nstride,
                                     color_map.aColor);
    }
  }
  if( elm.shape()[1] == 4 ){
    if( ndim == 3 ){
      opengl::DrawMeshTet3D_ScalarP1(pos.data(), np,
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
      opengl::DrawMeshTet3D_FaceNormDisp(pos.data(), np,
                                         elm.data(), nelm,
                                         disp.data());
    }
  }
  else if( ndim == 2 ){
    if( meshelem_type == MESHELEM_TRI ){
      opengl::DrawMeshTri2D_FaceDisp2D(pos.data(), np,
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
    opengl::DrawPoints2D_Vectors(pos.data(), np,
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
  
  // ----------------
  // gl misc
  m.def("setSomeLighting",  &opengl::setSomeLighting, "set some lighting that looks good for me");
  m.def("setup_glsl",       &setUpGLSL, "compile shader program");
  m.def("glew_init",        &PyGLExtensionInit);
  
  m.def("cppDrawSphere",      &opengl::DrawSphereAt );
  m.def("cppDrawSphere_Edge", &opengl::DrawSphere_Edge);
  m.def("cppDrawTorus_Edge",  &opengl::DrawTorus_Edge);
}


