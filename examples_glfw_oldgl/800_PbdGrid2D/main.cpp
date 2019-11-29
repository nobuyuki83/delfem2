#include <cmath>
#include "delfem2/mshmisc.h"
#include "delfem2/mshtopo.h"
#include "delfem2/primitive.h"
//
#include "delfem2/objfunc_v23.h"

// --------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/gl2_funcs.h"

namespace dfm2 = delfem2;

// -------------------------------------------------

void stepTime
(std::vector<double>& aXY1,
 std::vector<double>& aUV1,
 std::vector<double>& aTmp,
 double dt,
 int nitr,
 const std::vector<unsigned int>& clstr_ind,
 const std::vector<unsigned int>& clstr,
 const std::vector<int>& aBC,
 const std::vector<unsigned int>& aQuad,
 const std::vector<double>& aXY0)
{
  const int ndof = aXY0.size();
  for (int idof=0; idof<ndof; idof++){
    aTmp[idof] = aXY1[idof]+aUV1[idof]*dt;
  }
  for(size_t ip=0;ip<aXY0.size()/2;++ip){
    if( aBC[ip] == 0 ){ continue; }
    aTmp[ip*2+0] = aXY1[ip*2+0];
    aTmp[ip*2+1] = aXY1[ip*2+1];
  }
  { // deform
    for (int itr=0; itr<nitr; itr++){
      PBD_ConstProj_Rigid2D(aTmp.data(),
                            0.5,
                            clstr_ind.data(), clstr_ind.size(),
                            clstr.data(), clstr.size(),
                            aXY0.data(), aXY0.size());
    }
  }
  for(size_t ip=0;ip<aXY0.size()/2;++ip){
    if( aBC[ip] == 0 ){ continue; }
    aTmp[ip*2+0] = aXY1[ip*2+0];
    aTmp[ip*2+1] = aXY1[ip*2+1];
  }
  for (int idof=0; idof<ndof; ++idof){
    aUV1[idof] = (aTmp[idof]-aXY1[idof])*(1.0/dt);
  }
  for (int idof=0; idof<ndof; idof++){
    aXY1[idof] = aTmp[idof];
  }
}

const int nX = 8;
const int nY = 8;
std::vector<double> aXY0;
std::vector<double> aXY1;
std::vector<double> aUV1;
std::vector<double> aXYt;
std::vector<unsigned int> aQuad;
std::vector<unsigned int> clstr_ind, clstr;
std::vector<int> aBC;

void myGlutDisplay()
{
  {
    int viewport[8];
    glGetIntegerv(GL_VIEWPORT, viewport);
    double w = (double)viewport[2];
    double h = (double)viewport[3];
    double asp = w/h;
    const double win_size = 10.0;
    ::glMatrixMode(GL_PROJECTION);
    ::glLoadIdentity();
    ::glOrtho(-asp*win_size, +asp*win_size, -win_size, +win_size, -10, +10);
    ::glMatrixMode(GL_MODELVIEW);
    ::glLoadIdentity();
  }
  
  delfem2::opengl::DrawMeshQuad2D_Edge(aXY1.data(), aXY1.size()/2,
                              aQuad.data(), aQuad.size()/4);
}

void myGlutIdle(){
  double dt = 0.1;
  static double t = 0;
  t += dt;
  for(int ix=0;ix<nX+1;++ix){
    aXY1[ix*2+0] = ix + 2*sin(t*2);
    aXY1[ix*2+1] = 0;
  }
  stepTime(aXY1, aUV1, aXYt,
           dt, 1,
           clstr_ind, clstr,
           aBC,
           aQuad, aXY0);
}

int main(int argc,char* argv[])
{
  // --------------------------
  
  delfem2::MeshQuad2D_Grid(aXY0, aQuad, nX, nY);
  aXY1 = aXY0;
  aXYt = aXY0;
  aUV1.resize(aXY0.size());
  
  std::vector<unsigned int> psup_ind, psup;
  JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                                              aQuad.data(), aQuad.size()/4, 4,
                                              aXY0.size()/2);
    //  Print_IndexedArray(psup_ind, psup);
  dfm2::JArray_AddDiagonal(clstr_ind, clstr,
                           psup_ind.data(), psup_ind.size(),  psup.data(), psup.size());
    //  JArray_Print(clstr_ind, clstr);
  
  aBC.assign(aXY0.size()/2,0);
  for(int ix=0;ix<nX+1;++ix){
    aBC[ix] = 1;
  }
    
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    {
      static int iframe = 0;
      if( iframe % 5 == 0 ){
        myGlutIdle();
      }
      iframe++;
    }
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

