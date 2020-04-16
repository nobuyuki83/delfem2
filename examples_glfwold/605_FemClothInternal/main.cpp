/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include "delfem2/mshmisc.h"
#include "delfem2/mshtopo.h"
#include "delfem2/mats.h"
// -------
#include "delfem2/ilu_mats.h"
#include "delfem2/fem_emats.h"
#include "delfem2/cloth_internal.h"

// ------------------------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

/* ------------------------------------------------------------------------ */


class CInput_ContactNothing: public dfm2::CInput_Contact
{
public:
  double penetrationNormal(double& nx, double &ny, double& nz,
                           double px, double py, double pz) const
  {
    return -100;
  }
};


class CInput_ContactPlane: public dfm2::CInput_Contact
{
  double penetrationNormal(double& nx, double &ny, double& nz,
                                double px, double py, double pz) const
  {
    nx = 0.0;  ny = 0.0;  nz = 1.0; // normal of the plane
    return -0.5 - pz; // penetration depth
  }
};

class CInput_ContactSphere: public dfm2::CInput_Contact
{
  double penetrationNormal(double& nx, double &ny, double& nz,
                           double px, double py, double pz) const
  {
    const double center[3] = { 0.1, 0.5, -0.8 };
    const double radius = 0.3;
    nx = px-center[0];
    ny = py-center[1];
    nz = pz-center[2];
    double len = sqrt(nx*nx+ny*ny+nz*nz);
    nx /= len;
    ny /= len;
    nz /= len;
    return radius-len; // penetration depth
  }
};


/* ------------------------------------------------------------------------ */
// input parameter for simulation

const int ndiv = 25;  // (in) number of division of the square cloth edge
const double cloth_size = 1; // square cloth 1m x 1m
std::vector<double> aXYZ0; // undeformed vertex positions
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aUVW; // deformed vertex velocity
std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)
std::vector<unsigned int> aTri;  // index of triangles
std::vector<unsigned int> aQuad; // index of 4 vertices required for bending
const double lambda = 1.0; // Lame's 1st parameter
const double myu    = 4.0; // Lame's 2nd parameter
const double stiff_bend = 1.0e-3; // bending stiffness
const double areal_density = 1.0; // areal density of a cloth
const double gravity[3] = {0,0,-10}; // gravitatinal accereration
double time_step_size = 0.02; // size of time step
const double stiff_contact = 1.0e+3;
const double contact_clearance = 0.02;
int imode_contact; // mode of contacting object
double mass_point; // mass for a point

// variables for sparse solver
dfm2::CMatrixSparse<double> mat_A; // coefficient matrix
dfm2::CPreconditionerILU<double>  ilu_A; // ilu decomposition of the coefficient matrix

bool is_animation;
int imode_draw = 0;

/* ------------------------------------------------------------------------ */


void StepTime()
{
  CInput_ContactPlane c0;
  CInput_ContactNothing c1;
  CInput_ContactSphere c2;
  dfm2::CInput_Contact* ic = 0;
  if(      imode_contact == 1 ){ ic = &c0; }
  if(      imode_contact == 0 ){ ic = &c1; }
  if(      imode_contact == 2 ){ ic = &c2; }
  // solving lienar system using conjugate gradient method with ILU(0) preconditioner
  StepTime_InternalDynamicsILU(aXYZ, aUVW, mat_A, ilu_A,
                               aXYZ0, aBCFlag,
                               aTri, aQuad,
                               time_step_size,
                               lambda, myu, stiff_bend,
                               gravity, mass_point,
                               stiff_contact,contact_clearance,*ic);
//  MakeNormal(aNormal, aXYZ, aTri);
}

void myGlutDisplay(void)
{
  //	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  
  bool is_lighting = glIsEnabled(GL_LIGHTING);
  
  ::glColor3d(0,0,0);
  delfem2::opengl::DrawMeshTri3D_Edge(aXYZ, aTri);
  {
    ::glEnable(GL_LIGHTING);
    float color[4] = {200.0/256.0, 200.0/256.0, 200.0/256.0,1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
  }
  delfem2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aTri);
  
  { // fixed boundary condition
    ::glDisable(GL_LIGHTING);        
    ::glPointSize(5);
    ::glColor3d(0,0,1);
    ::glBegin(GL_POINTS);
    for(int ip=0;ip<aXYZ.size()/3;ip++){
      if( aBCFlag[ip*3+0] == 0 && aBCFlag[ip*3+1] == 0 && aBCFlag[ip*3+2] == 0 ) continue;
      ::glVertex3d(aXYZ0[ip*3+0],aXYZ0[ip*3+1],aXYZ0[ip*3+2]);
    }
    ::glEnd();
  }
  
  if(      imode_contact == 1 ){     // draw floor
    ::glDisable(GL_LIGHTING);        
    ::glLineWidth(1);
    ::glColor3d(1,0,0);
    ::glBegin(GL_LINES);
    double grid_x_min = -10;
    double grid_x_max = +10;
    double grid_y_min = -10;
    double grid_y_max = +10;
    int ndiv_grid = 30;
    for(unsigned int ix=0;ix<ndiv_grid+1;ix++){
      double x0 = (grid_x_max-grid_x_min) / ndiv_grid * ix + grid_x_min;
      ::glVertex3d(x0,grid_y_min,-0.5);
      ::glVertex3d(x0,grid_y_max,-0.5);
    }
    for(unsigned int iz=0;iz<ndiv_grid+1;iz++){
      double z0 = (grid_y_max-grid_y_min) / ndiv_grid * iz + grid_y_min;
      ::glVertex3d(grid_x_min,z0,-0.5);
      ::glVertex3d(grid_x_max,z0,-0.5);
    }
    ::glEnd();        
  }
  else if( imode_contact == 2 ){
    ::glDisable(GL_LIGHTING);        
    ::glLineWidth(1);
    ::glColor3d(1,0,0);    
    ::glPushMatrix();
    ::glTranslated(0.1, 0.5, -0.8);
//    ::glutWireSphere(0.3, 16, 16);
    ::glPopMatrix();
  }
  
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }
  else{              ::glDisable(GL_LIGHTING); }  
  
}


int main(int argc,char* argv[])
{
  { // initialze data
    SetClothShape_Square(aXYZ0,aBCFlag,aTri,aQuad,
                         ndiv,cloth_size);
    const int np = aXYZ0.size()/3.0;
    double total_area = cloth_size*cloth_size;
    mass_point = total_area*areal_density / (double)np;
    // initialize deformation
    aXYZ = aXYZ0;
    aUVW.assign(np*3,0.0);
//    MakeNormal(aNormal, aXYZ, aTri);
    mat_A.Initialize(np,3,true);
    std::vector<unsigned int> psup_ind,psup;
    dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                                                      aQuad.data(),aQuad.size()/4, 4, np);
    dfm2::JArray_Sort(psup_ind, psup);
    mat_A.SetPattern(psup_ind.data(),psup_ind.size(),
                     psup.data(),psup.size());
    ilu_A.Initialize_ILU0(mat_A);
  }
  
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  delfem2::opengl::setSomeLighting();
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);
  delfem2::opengl::setSomeLighting();
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_ZTOP;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_ZTOP;
  viewer.nav.camera.psi = 3.1415*0.2;
  viewer.nav.camera.theta = 3.1415*0.1;
  viewer.nav.camera.view_height = 2;
  while(!glfwWindowShouldClose(viewer.window)) {
    StepTime();
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


