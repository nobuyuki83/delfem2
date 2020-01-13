//
//  py_rigidasm.cpp
//  c_gl
//
//  Created by Nobuyuki Umetani on 2019-11-01.
//

#include <cstdio>
#include <iostream>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "delfem2/eigen/eigen_rigidbody.h"

// -------------
#if defined(__APPLE__) && defined(__MACH__)
  #include <OpenGL/gl.h>
#elif defined(_WIN32) // windows
  #include <windows.h>
  #include <GL/gl.h>
#else
  #include <GL/gl.h>
#endif
#include "delfem2/opengl/glold_v23.h"
#include "delfem2/opengl/glold_funcs.h"


namespace py = pybind11;
namespace dfm2 = delfem2;


void DrawRigidBodyAssemblyStatic(const CRigidBodyAssembly_Static& rba)
{
  const double small_rad = 0.1;
  const double big_rad = 0.1;
  ::glDisable(GL_LIGHTING);
  const std::vector<CRigidBody>& aRigidBody = rba.aRigidBody;
  for(auto rb : aRigidBody){
    CVector3 cg = rb.cg;
    if( rba.is_draw_deformed ){
      cg += rb.u;
    }
    
    if( rba.is_draw_skeleton ){
      ::glColor3d(0,1,0);
      ::glPushMatrix();
      ::glTranslated(cg.x,cg.y,cg.z);
        //::glutWireSphere(0.1, 16, 16);
      //gluSphere(quadricSphere, big_rad, 16, 16);
      dfm2::opengl::DrawSphereAt(16, 16, big_rad, 0.0, 0.0, 0.0);
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
      //gluSphere(quadricSphere, small_rad, 16, 16);
      dfm2::opengl::DrawSphereAt(16, 16, small_rad, 0.0, 0.0, 0.0);
      ::glPopMatrix();
      
      if( rba.is_draw_skeleton ){
        ::glLineWidth(5);
        ::glColor3d(0,0,0);
        ::glBegin(GL_LINES);
        dfm2::opengl::myGlVertex(cg);
        dfm2::opengl::myGlVertex( p);
        ::glEnd();
      }
      
      if( rba.is_draw_force ){
        if (!rb.aCForce.empty()) {
          CVector3 q = p + rba.scale_force*rb.aCForce[icp];
          ::glLineWidth(5);
          ::glColor3d(1,0,0);
          ::glBegin(GL_LINES);
          dfm2::opengl::myGlVertex( p);
          dfm2::opengl::myGlVertex( q);
          ::glEnd();
        }
      }
    }
  }
  // -------------------------------------------
  const std::vector<CJoint>& aJoint = rba.aJoint;
  for(const auto & joint : aJoint){
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
    // gluSphere(quadricSphere, small_rad, 16, 16);
    dfm2::opengl::DrawSphereAt(16, 16, small_rad, 0.0, 0.0, 0.0);
    ::glPopMatrix();
    
      // joint point seen from rb1
    ::glPushMatrix();
    ::glTranslated(p1.x,p1.y,p1.z);
      //::glutWireSphere(0.02, 16, 16);
      //    gluSphere(quadricSphere, small_rad, 16, 16);
    dfm2::opengl::DrawSphereAt(16, 16, small_rad, 0.0, 0.0, 0.0);
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
      dfm2::opengl::myGlVertex(cg0);
      dfm2::opengl::myGlVertex( p0);
      dfm2::opengl::myGlVertex(cg1);
      dfm2::opengl::myGlVertex( p1);
      ::glEnd();
    }
    
    if( rba.is_draw_force ){
      ::glLineWidth(5);
      CVector3 q0 = p0 + rba.scale_force*joint.linear;
      CVector3 q1 = p1 - rba.scale_force*joint.linear;
      ::glColor3d(1,0,1);
      ::glBegin(GL_LINES);
      dfm2::opengl::myGlVertex( p0);
      dfm2::opengl::myGlVertex( q0);
      dfm2::opengl::myGlVertex( p1);
      dfm2::opengl::myGlVertex( q1);
      ::glEnd();
      CVector3 r0 = p0 + rba.scale_torque*joint.torque;
      CVector3 r1 = p1 - rba.scale_torque*joint.torque;
      ::glColor3d(0,1,1);
      ::glBegin(GL_LINES);
      dfm2::opengl::myGlVertex( p0);
      dfm2::opengl::myGlVertex( r0);
      dfm2::opengl::myGlVertex( p1);
      dfm2::opengl::myGlVertex( r1);
      ::glEnd();
    }
  }
  glLineWidth(1.0f);
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




void init_rigidasm(py::module &m)
{
  m.def("drawRigidBodyAssemblyStatic", &DrawRigidBodyAssemblyStatic);
}
