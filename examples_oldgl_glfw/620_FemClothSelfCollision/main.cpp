
#include <iostream>
#include <vector>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/cloth_internal.h"
#include "delfem2/cloth_selfcollision.h"
#include "delfem2/vec3.h"
#include "delfem2/srchbvh.h"
#include "delfem2/mshuni.h"
#include "delfem2/jagarray.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"

namespace dfm2 = delfem2;

/* ------------------------------------------------------------------------ */


/**
 *
 * @param[out] aXYZ0 undeformed vertex positions
 * @param[out] aBCFlag boundary condition flag (0:free 1:fixed)
 * @param[out] aTri index of triangles
 * @param[out] aQuad index of 4 vertices required for bending
 * @param[out] total_area total area of cloth
 * @param[in] elem_length number of division of the square cloth edge
 * @param cloth_size_x
 * @param[in] cloth_size_z  size of square cloth
 */
void SetClothShape_Square(
    std::vector<double> &aXYZ0, //
    std::vector<int> &aBCFlag, //
    std::vector<unsigned int> &aTri, //
    std::vector<unsigned int> &aQuad, //
    double &total_area, //
    //
    double elem_length, //
    double cloth_size_x,
    double cloth_size_z) {
  // make vertex potision array 頂点位置配列を作る
  const int ndiv_x = (int) (cloth_size_x / elem_length); // size of an element
  const int ndiv_z = (int) (cloth_size_z / elem_length); // size of an element
  const int nxyz = (ndiv_x + 1) * (ndiv_z + 1); // number of points
  aXYZ0.reserve(nxyz * 3);
  for (int ix = 0; ix < ndiv_x + 1; ix++) {
    for (int iz = 0; iz < ndiv_z + 1; iz++) {
      aXYZ0.push_back(ix * elem_length);
      aXYZ0.push_back(0.0);
      aXYZ0.push_back(iz * elem_length);
    }
  }

  // make triangle index array, 三角形の頂点配列を作る
  const int ntri = ndiv_x * ndiv_z * 2;
  aTri.reserve(ntri * 3);
  for (int ix = 0; ix < ndiv_x; ix++) {
    for (int iz = 0; iz < ndiv_z; iz++) {
      aTri.push_back(ix * (ndiv_z + 1) + iz);
      aTri.push_back((ix + 1) * (ndiv_z + 1) + iz);
      aTri.push_back(ix * (ndiv_z + 1) + (iz + 1));
      ////
      aTri.push_back((ix + 1) * (ndiv_z + 1) + (iz + 1));
      aTri.push_back(ix * (ndiv_z + 1) + (iz + 1));
      aTri.push_back((ix + 1) * (ndiv_z + 1) + iz);
    }
  }

  aQuad.reserve(ndiv_x * ndiv_z + ndiv_x * (ndiv_z - 1) + (ndiv_x - 1) * ndiv_z);
  for (int ix = 0; ix < ndiv_x; ix++) {
    for (int iz = 0; iz < ndiv_z; iz++) {
      aQuad.push_back((ix + 0) * (ndiv_z + 1) + (iz + 0));
      aQuad.push_back((ix + 1) * (ndiv_z + 1) + (iz + 1));
      aQuad.push_back((ix + 1) * (ndiv_z + 1) + (iz + 0));
      aQuad.push_back((ix + 0) * (ndiv_z + 1) + (iz + 1));
    }
  }
  for (int ix = 0; ix < ndiv_x; ix++) {
    for (int iz = 0; iz < ndiv_z - 1; iz++) {
      aQuad.push_back((ix + 1) * (ndiv_z + 1) + (iz + 0));
      aQuad.push_back((ix + 0) * (ndiv_z + 1) + (iz + 2));
      aQuad.push_back((ix + 1) * (ndiv_z + 1) + (iz + 1));
      aQuad.push_back((ix + 0) * (ndiv_z + 1) + (iz + 1));
    }
  }
  for (int ix = 0; ix < ndiv_x - 1; ix++) {
    for (int iz = 0; iz < ndiv_z; iz++) {
      aQuad.push_back((ix + 0) * (ndiv_z + 1) + (iz + 1));
      aQuad.push_back((ix + 2) * (ndiv_z + 1) + (iz + 0));
      aQuad.push_back((ix + 1) * (ndiv_z + 1) + (iz + 0));
      aQuad.push_back((ix + 1) * (ndiv_z + 1) + (iz + 1));
    }
  }

  // compute total area 全体の面積を計算
  total_area = cloth_size_x * cloth_size_z;

  // set fixed boundary condition 固定境界条件を指定する
  aBCFlag = std::vector<int>(nxyz * 3, 0);
  /*
  for(int ix=0;ix<ndiv+1;ix++){
    aBCFlag[ ix*(ndiv+1) + ndiv ] = 1;
  }
   */
}





/* ------------------------------------------------------------------------ */
// input parameter for simulation
double elem_length = 0.1;
double cloth_size_x = 0.4;
double cloth_size_z = 5.0;
std::vector<double> aXYZ0; // undeformed vertex positions
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aUVW; // deformed vertex velocity
std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)
std::vector<unsigned int> aTri;  // index of triangles
std::vector<unsigned int> aQuad; // index of 4 vertices required for bending
const double lambda = 0.0; // Lame's 1st parameter
const double myu = 30.0; // Lame's 2nd parameter
const double stiff_bend = 1.0e-2; // bending stiffness
const double areal_density = 0.04; // areal density of a cloth,
const double gravity[3] = {0, 0.05, -10}; // gravitatinal accereration
double time_step_size = 0.01; // size of time step
const double stiff_contact = 300;
const double contact_clearance = 0.01;
int imode_contact; // mode of contacting object
double mass_point; // mass for a point

// TODO: make variables non-global
dfm2::CMatrixSparse<double> mat_A;
dfm2::CPreconditionerILU<double> ilu_A;

std::vector<double> aNormal; // deformed vertex noamals，変形中の頂点の法線(可視化用)

// variables for self-collision
int iroot_bvh; // index of the root BVH node
std::vector<delfem2::CNodeBVH2> aNodeBVH; // array of BVH node
std::vector<delfem2::CBV3d_AABB> aBB_BVH; // array of AABB
//CJaggedArray aEdge;
std::vector<unsigned int> psup_ind, psup;

// data for camera
bool is_animation;
int imode_draw = 0;
/* ------------------------------------------------------------------------ */


void MakeNormal() { // make normal
  const size_t np = aXYZ.size() / 3;
  const size_t ntri = aTri.size() / 3;
  aNormal.assign(np * 3, 0);
  for (unsigned int itri = 0; itri < ntri; itri++) {
    const unsigned int ip0 = aTri[itri * 3 + 0];
    const unsigned int ip1 = aTri[itri * 3 + 1];
    const unsigned int ip2 = aTri[itri * 3 + 2];
    double c[3][3] = {
        {aXYZ[ip0 * 3 + 0], aXYZ[ip0 * 3 + 1], aXYZ[ip0 * 3 + 2]},
        {aXYZ[ip1 * 3 + 0], aXYZ[ip1 * 3 + 1], aXYZ[ip1 * 3 + 2]},
        {aXYZ[ip2 * 3 + 0], aXYZ[ip2 * 3 + 1], aXYZ[ip2 * 3 + 2]}};
    double n[3], area;
    dfm2::UnitNormalAreaTri3(n, area, c[0], c[1], c[2]);
    aNormal[ip0 * 3 + 0] += n[0];
    aNormal[ip0 * 3 + 1] += n[1];
    aNormal[ip0 * 3 + 2] += n[2];
    aNormal[ip1 * 3 + 0] += n[0];
    aNormal[ip1 * 3 + 1] += n[1];
    aNormal[ip1 * 3 + 2] += n[2];
    aNormal[ip2 * 3 + 0] += n[0];
    aNormal[ip2 * 3 + 1] += n[1];
    aNormal[ip2 * 3 + 2] += n[2];
  }
  for (unsigned int ip = 0; ip < np; ip++) {
    double sqlen =
        +aNormal[ip * 3 + 0] * aNormal[ip * 3 + 0]
            + aNormal[ip * 3 + 1] * aNormal[ip * 3 + 1]
            + aNormal[ip * 3 + 2] * aNormal[ip * 3 + 2];
    double invlen = 1.0 / sqrt(sqlen);
    aNormal[ip * 3 + 0] *= invlen;
    aNormal[ip * 3 + 1] *= invlen;
    aNormal[ip * 3 + 2] *= invlen;
  }
}

class CInput_ContactPlane : public dfm2::CInput_Contact {
  double penetrationNormal(double &nx, double &ny, double &nz,
                           [[maybe_unused]] double px,
                           [[maybe_unused]] double py, double pz) const override {
    nx = 0.0;
    ny = 0.0;
    nz = 1.0; // normal of the plane
    return -0.5 - pz; // penetration depth
  }
};

class CInput_ContactSphere : public dfm2::CInput_Contact {
  double penetrationNormal(double &nx, double &ny, double &nz,
                           double px, double py, double pz) const override {
    const double center[3] = {0.5, 0.0, +2};
    const double radius = 3.0;
    nx = px - center[0];
    ny = py - center[1];
    nz = pz - center[2];
    double len = sqrt(nx * nx + ny * ny + nz * nz);
    nx /= -len;
    ny /= -len;
    nz /= -len;
    return -(radius - len); // penetration depth
  }
};

void StepTime() {
  CInput_ContactPlane c0;
  CInput_ContactSphere c1;
  dfm2::CInput_Contact *ic = 0;
  if (imode_contact == 0) { ic = &c0; }
  if (imode_contact == 1) { ic = &c1; }
  //
  std::vector<double> aXYZ1 = aXYZ;
  ::StepTime_InternalDynamicsILU(aXYZ, aUVW, mat_A, ilu_A,
                                 aXYZ0, aBCFlag,
                                 aTri, aQuad,
                                 time_step_size,
                                 lambda, myu, stiff_bend,
                                 gravity, mass_point,
                                 stiff_contact, contact_clearance, *ic);
  //
  bool is_impulse_applied;
  GetIntermidiateVelocityContactResolved(aUVW,
                                         is_impulse_applied,
                                         time_step_size,
                                         contact_clearance,
                                         mass_point,
                                         stiff_contact,
                                         aXYZ1,
                                         aTri,
                                         psup_ind, psup,
                                         iroot_bvh, aNodeBVH, aBB_BVH);
  if (is_impulse_applied) {
    std::cout << "update middle velocity" << std::endl;
    for (unsigned int ip = 0; ip < aXYZ.size() / 3; ip++) {
      aXYZ[ip * 3 + 0] = aXYZ1[ip * 3 + 0] + aUVW[ip * 3 + 0] * time_step_size;
      aXYZ[ip * 3 + 1] = aXYZ1[ip * 3 + 1] + aUVW[ip * 3 + 1] * time_step_size;
      aXYZ[ip * 3 + 2] = aXYZ1[ip * 3 + 2] + aUVW[ip * 3 + 2] * time_step_size;
    }
    // now aXYZ is collision free
    ::UpdateIntermidiateVelocity(aUVW,
                                 aXYZ, aXYZ0, aBCFlag,
                                 aTri, aQuad,
                                 time_step_size,
                                 lambda, myu, stiff_bend,
                                 gravity, mass_point);
  }
  MakeNormal();
}

void myGlutDisplay() {

  bool is_lighting = glIsEnabled(GL_LIGHTING);

  { // draw triangle
    ::glDisable(GL_LIGHTING);
    ::glLineWidth(1);
    ::glColor3d(0, 0, 0);
    ::glBegin(GL_LINES);
    for (unsigned int itri = 0; itri < aTri.size() / 3; itri++) {
      const unsigned int ip0 = aTri[itri * 3 + 0];
      const unsigned int ip1 = aTri[itri * 3 + 1];
      const unsigned int ip2 = aTri[itri * 3 + 2];
      ::glVertex3d(aXYZ[ip0 * 3 + 0], aXYZ[ip0 * 3 + 1], aXYZ[ip0 * 3 + 2]);
      ::glVertex3d(aXYZ[ip1 * 3 + 0], aXYZ[ip1 * 3 + 1], aXYZ[ip1 * 3 + 2]);
      ::glVertex3d(aXYZ[ip1 * 3 + 0], aXYZ[ip1 * 3 + 1], aXYZ[ip1 * 3 + 2]);
      ::glVertex3d(aXYZ[ip2 * 3 + 0], aXYZ[ip2 * 3 + 1], aXYZ[ip2 * 3 + 2]);
      ::glVertex3d(aXYZ[ip2 * 3 + 0], aXYZ[ip2 * 3 + 1], aXYZ[ip2 * 3 + 2]);
      ::glVertex3d(aXYZ[ip0 * 3 + 0], aXYZ[ip0 * 3 + 1], aXYZ[ip0 * 3 + 2]);
    }
    ::glEnd();
  }
  { // draw triangle face
    ::glEnable(GL_LIGHTING);
    float color[4] = {200.0 / 256.0, 200.0 / 256.0, 200.0 / 256.0, 1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
    glShadeModel(GL_SMOOTH);
    ::glBegin(GL_TRIANGLES);
    for (unsigned int itri = 0; itri < aTri.size() / 3; itri++) {
      const unsigned int ip0 = aTri[itri * 3 + 0];
      const unsigned int ip1 = aTri[itri * 3 + 1];
      const unsigned int ip2 = aTri[itri * 3 + 2];
      double c[3][3] = {
          {aXYZ[ip0 * 3 + 0], aXYZ[ip0 * 3 + 1], aXYZ[ip0 * 3 + 2]},
          {aXYZ[ip1 * 3 + 0], aXYZ[ip1 * 3 + 1], aXYZ[ip1 * 3 + 2]},
          {aXYZ[ip2 * 3 + 0], aXYZ[ip2 * 3 + 1], aXYZ[ip2 * 3 + 2]}};
      ::glNormal3d(aNormal[ip0 * 3 + 0], aNormal[ip0 * 3 + 1], aNormal[ip0 * 3 + 2]);
      ::glVertex3dv(c[0]);
      ::glNormal3d(aNormal[ip1 * 3 + 0], aNormal[ip1 * 3 + 1], aNormal[ip1 * 3 + 2]);
      ::glVertex3dv(c[1]);
      ::glNormal3d(aNormal[ip2 * 3 + 0], aNormal[ip2 * 3 + 1], aNormal[ip2 * 3 + 2]);
      ::glVertex3dv(c[2]);
    }
    ::glEnd();
  }

  { // fixed boundary condition
    ::glDisable(GL_LIGHTING);
    ::glPointSize(5);
    ::glColor3d(0, 0, 1);
    ::glBegin(GL_POINTS);
    for (unsigned int ip = 0; ip < aXYZ.size() / 3; ip++) {
      if (aBCFlag[ip] == 0) continue;
      ::glVertex3d(aXYZ0[ip * 3 + 0], aXYZ0[ip * 3 + 1], aXYZ0[ip * 3 + 2]);
    }
    ::glEnd();
  }

  if (imode_contact == 0) {     // draw floor
    ::glDisable(GL_LIGHTING);
    ::glLineWidth(1);
    ::glColor3d(1, 0, 0);
    ::glBegin(GL_LINES);
    double grid_x_min = -10;
    double grid_x_max = +10;
    double grid_y_min = -10;
    double grid_y_max = +10;
    const unsigned int ndiv_grid = 30;
    for (unsigned int ix = 0; ix < ndiv_grid + 1; ix++) {
      double x0 = (grid_x_max - grid_x_min) / ndiv_grid * ix + grid_x_min;
      ::glVertex3d(x0, grid_y_min, -0.5);
      ::glVertex3d(x0, grid_y_max, -0.5);
    }
    for (unsigned int iz = 0; iz < ndiv_grid + 1; iz++) {
      double z0 = (grid_y_max - grid_y_min) / ndiv_grid * iz + grid_y_min;
      ::glVertex3d(grid_x_min, z0, -0.5);
      ::glVertex3d(grid_x_max, z0, -0.5);
    }
    ::glEnd();
  } else if (imode_contact == 1) {
    ::glDisable(GL_LIGHTING);
    ::glLineWidth(1);
    ::glColor3d(1, 0, 0);
    ::glPushMatrix();
//    const double radius = 3.0;
    ::glTranslated(0.5, 0.0, +2.0);
//    ::glutWireSphere(radius, 16, 16);
    ::glPopMatrix();
  }

  ::glColor3d(0, 0, 0);
//  ShowFPS();

  if (is_lighting) { ::glEnable(GL_LIGHTING); }
  else { ::glDisable(GL_LIGHTING); }

//  ::glutSwapBuffers();
}

int main() {
  { // initialze data
    double total_area;
    SetClothShape_Square(
        aXYZ0, aBCFlag, aTri, aQuad, total_area,
        elem_length, cloth_size_x, cloth_size_z);
    const size_t np = aXYZ0.size() / 3;
    mass_point = total_area * areal_density / (double) np;
    // initialize deformation
    aXYZ = aXYZ0;
    aUVW.assign(np * 3, 0.0);
    MakeNormal();
    //
    {
      const size_t ntri = aTri.size() / 3;
      std::vector<double> aElemCenter(ntri * 3);
      for (unsigned int itri = 0; itri < ntri; ++itri) {
        dfm2::CVec3d p = dfm2::CG_Tri3(itri, aTri, aXYZ);
        aElemCenter[itri * 3 + 0] = p.x;
        aElemCenter[itri * 3 + 1] = p.y;
        aElemCenter[itri * 3 + 2] = p.z;
      }
      std::vector<unsigned int> aTriSuTri;
      dfm2::ElSuEl_MeshElem(
          aTriSuTri,
          aTri.data(), aTri.size() / 3,
          delfem2::MESHELEM_TRI, aXYZ.size() / 3);
      iroot_bvh = dfm2::BVHTopology_TopDown_MeshElem(
          aNodeBVH,
          3, aTriSuTri, aElemCenter);
    }
    std::cout << "aNodeBVH.size(): " << aNodeBVH.size() << std::endl;
//    aEdge.SetEdgeOfElem(aTri,(int)aTri.size()/3,3, np,false);
    dfm2::JArray_PSuP_MeshElem(
        psup_ind, psup,
        aQuad.data(), aQuad.size() / 4, 4, np);
    dfm2::JArray_Sort(psup_ind, psup);
    mat_A.Initialize(np, 3, true);
    /*
    CJaggedArray crs;
    crs.SetEdgeOfElem(aQuad, (int)aQuad.size()/4, 4, np, false);
    crs.Sort();
     */
    mat_A.SetPattern(
        psup_ind.data(), psup_ind.size(),
        psup.data(), psup.size());
    ilu_A.SetPattern0(mat_A);
  }

  // ---------------------------


  delfem2::glfw::CViewer3 viewer;
  //viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::ZTOP;
  //viewer.camera.psi = 3.1415 * 0.2;
  //viewer.camera.theta = 3.1415 * 0.1;
  viewer.projection.view_height = 2.0;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  delfem2::opengl::setSomeLighting();
  while (!glfwWindowShouldClose(viewer.window)) {
    StepTime();
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


