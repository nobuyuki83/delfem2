#include <iostream>
#include <math.h>
#include "delfem2/mshmisc.h"
#include "delfem2/mshtopo.h"
#include "delfem2/primitive.h"
#include "delfem2/cloth_internal.h"

#if defined(_MSC_VER)
  #include <windows.h>
#endif

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#endif

#include "delfem2/opengl/gl24_funcs.h"
#include "delfem2/opengl/gl4_funcs.h"
#include "delfem2/opengl/gl4_mshcolor.h"
#include "delfem2/opengl/glfw_cam.h"
#include "delfem2/opengl/glfw_viewer.hpp"

namespace dfm2 = delfem2;

// ---------------------------------------


class CInput_ContactNothing: public CInput_Contact
{
public:
  double penetrationNormal(double& nx, double &ny, double& nz,
                           double px, double py, double pz) const
  {
    return -100;
  }
};


class CInput_ContactPlane: public CInput_Contact
{
  double penetrationNormal(double& nx, double &ny, double& nz,
                           double px, double py, double pz) const
  {
    nx = 0.0;  ny = 0.0;  nz = 1.0; // normal of the plane
    return -0.5 - pz; // penetration depth
  }
};

class CInput_ContactSphere: public CInput_Contact
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

// ---------------------------------------------------------

dfm2::opengl::CViewer_GLFW viewer;
CShader_TriMesh shdr_trimsh;

std::vector<double> aXYZ0; // (out) undeformed vertex positions
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aUVW; // deformed vertex velocity
std::vector<int> aBCFlag; // (out) boundary condition flag (0:free 1:fixed)
std::vector<unsigned int> aTri; // (out) index of triangles
std::vector<unsigned int> aQuad; // (out) index of 4 vertices required for bending
dfm2::CMatrixSparse<double> mat_A; // coefficient matrix
dfm2::CPreconditionerILU<double>  ilu_A; // ilu decomposition of the coefficient matrix
double mass_point; // mass for a point

const int ndiv = 25;
const double cloth_size = 1;
const double areal_density = 1.0; // areal density of a cloth
const double lambda = 1.0; // Lame's 1st parameter
const double myu    = 4.0; // Lame's 2nd parameter
const double stiff_bend = 1.0e-3; // bending stiffness
const double gravity[3] = {0,0,-10}; // gravitatinal accereration
double time_step_size = 0.02; // size of time step
const double stiff_contact = 1.0e+3;
const double contact_clearance = 0.02;

// ---------------------------------------------------------

void StepTime()
{
  CInput_ContactPlane c0;
  CInput_ContactNothing c1;
  CInput_ContactSphere c2;
  CInput_Contact* ic = 0;
  int imode_contact = 0;
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

void draw(GLFWwindow* window)
{
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  
  StepTime();
  shdr_trimsh.UpdateVertex(aXYZ, aTri);

  float mP[16], mMV[16];
  viewer.nav.Matrix_MVP(mMV, mP, window);
  shdr_trimsh.Draw(mP,mMV);
  
  viewer.DrawEnd_oldGL();
}

int main(void)
{
  {
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
    std::vector<int> psup_ind,psup;
    JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                                                aQuad.data(),aQuad.size()/4, 4, np);
    dfm2::JArray_Sort(psup_ind, psup);
    mat_A.SetPattern(psup_ind.data(),psup_ind.size(),
                     psup.data(),psup.size());
    ilu_A.Initialize_ILU0(mat_A);
  }
  
  viewer.Init_newGL();
      
  // glad: load all OpenGL function pointers
  // ---------------------------------------
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
  
  shdr_trimsh.Compile();
  shdr_trimsh.Initialize(aXYZ, aTri);  
  
  viewer.nav.camera.view_height = 1.5;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
   
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

