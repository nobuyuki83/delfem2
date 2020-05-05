#include <iostream>
#include <cmath>
#include "delfem2/mshtopo.h"
#include "delfem2/mats.h"
#include "delfem2/vec2.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/ilu_mats.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/fem_emats.h"

// ----------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// ------------------------

void InitializeProblem_PlateBendingMITC3
(dfm2::CMatrixSparse<double>& mat_A,
 dfm2::CPreconditionerILU<double>& ilu_A,
 std::vector<int>& aBCFlag,
 const double lenx,
 const double leny,
 const std::vector<double>& aXY0,
 const std::vector<unsigned int>& aTri)
{
  const std::size_t np = aXY0.size()/2;
  aBCFlag.assign(np*3, 0);
  for(unsigned int ip=0;ip<np;++ip){
    const double px = aXY0[ip*2+0];
//    const double py = aXY0[ip*2+1];
    if( fabs(px-(-lenx*0.5)) < 0.0001 ){
      aBCFlag[ip*3+0] = 1;
      aBCFlag[ip*3+1] = 1;
      aBCFlag[ip*3+2] = 1;
    }
  }
  //
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                             aTri.data(), aTri.size()/3, 3,
                             (int)aXY0.size()/2);
  dfm2::JArray_Sort(psup_ind, psup);
  //
  mat_A.Initialize(np, 3, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(),psup.size());
//  ilu_A.Initialize_ILU0(mat_A);
    ilu_A.Initialize_ILUk(mat_A,0);
}

void SolveProblem_PlateBendingMITC3
(std::vector<double>& aVal,
 std::vector<double>& vec_b,
 std::vector<int>& aMSFlag,
 dfm2::CMatrixSparse<double>& mat_A,
 dfm2::CPreconditionerILU<double>& ilu_A,
 const std::vector<int>& aBCFlag,
 const double thickness,
 const double myu,
 const double lambda,
 const double rho,
 const double gravity_z,
 const std::vector<double>& aXY0,
 const std::vector<unsigned int>& aTri)
{
  const std::size_t np = aXY0.size()/2;
  const std::size_t nDoF = np*3;
  //
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_ShellStaticPlateBendingMITC3_MeshTri2D(mat_A,vec_b.data(),
                                                           thickness,lambda,myu,
                                                           rho,gravity_z,
                                                           aXY0.data(), aXY0.size()/2,
                                                           aTri.data(), aTri.size()/3,
                                                           aVal.data());
  std::cout << dfm2::Dot(vec_b, vec_b) << std::endl;
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag,0);
  //
  std::vector<double> vec_x;
  {
    ilu_A.SetValueILU(mat_A);
    ilu_A.DoILUDecomp();
    vec_x.resize(vec_b.size());
    std::vector<double> conv = Solve_PCG(vec_b.data(), vec_x.data(),
                                         vec_b.size(),
                                         1.0e-5, 1000,
                                         mat_A, ilu_A);
    std::cout << "convergence   nitr:" << conv.size() << "    res:" << conv[conv.size()-1] << std::endl;
  }
  //
  dfm2::XPlusAY(aVal,nDoF,aBCFlag,
                1.0,vec_x);
}

void myGlutDisplay
(const std::vector<double>& aXY0,
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aVal)
{
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  delfem2::opengl::DrawMeshTri2D_Edge(aTri,aXY0);
  {
    assert( aVal.size()/3 == aXY0.size()/2 );
    ::glColor3d(1,0,0);
    ::glBegin(GL_LINES);
    for(size_t itri=0;itri<aTri.size()/3;itri++){
      const unsigned int i0 = aTri[itri*3+0];
      const unsigned int i1 = aTri[itri*3+1];
      const unsigned int i2 = aTri[itri*3+2];
      const double p0[3] = { aXY0[i0*2+0], aXY0[i0*2+1], aVal[i0*3+0] };
      const double p1[3] = { aXY0[i1*2+0], aXY0[i1*2+1], aVal[i1*3+0] };
      const double p2[3] = { aXY0[i2*2+0], aXY0[i2*2+1], aVal[i2*3+0] };
      ::glVertex3dv( p0 );
      ::glVertex3dv( p1 );
      ::glVertex3dv( p1 );
      ::glVertex3dv( p2 );
      ::glVertex3dv( p2 );
      ::glVertex3dv( p0 );
    }
    ::glEnd();
  }
}

int main(int argc,char* argv[])
{
  const double lenx = 1.0;
  const double leny = 0.2;
  std::vector<unsigned int> aTri;
  std::vector<double> aXY0;
  {
    std::vector< std::vector<double> > aaXY;
    {
      aaXY.resize(1);
      aaXY[0].push_back(-lenx*0.5); aaXY[0].push_back(-leny*0.5);
      aaXY[0].push_back(+lenx*0.5); aaXY[0].push_back(-leny*0.5);
      aaXY[0].push_back(+lenx*0.5); aaXY[0].push_back(+leny*0.5);
      aaXY[0].push_back(-lenx*0.5); aaXY[0].push_back(+leny*0.5);
    }
    std::vector<dfm2::CDynPntSur> aPo2D;
    std::vector<dfm2::CDynTri> aETri;
    std::vector<dfm2::CVec2d> aVec2;
    GenMesh(aPo2D, aETri, aVec2,
            aaXY, 0.03, 0.03);
    MeshTri2D_Export(aXY0,aTri,
                     aVec2,aETri);
    std::cout<<"  ntri;"<<aTri.size()/3<<"  nXY:"<<aXY0.size()/2<<std::endl;
  }
  // -------------
  std::vector<double> aVal(aXY0.size()/2*3, 0.0);
  std::vector<int> aBCFlag; // boundary condition flag
  dfm2::CMatrixSparse<double> mat_A;
  dfm2::CPreconditionerILU<double> ilu_A;
  InitializeProblem_PlateBendingMITC3(mat_A,ilu_A,aBCFlag,
                                      lenx,leny,aXY0,aTri);
  // -----------
  const double thickness = 0.05;
  const double myu = 10000.0;
  const double lambda = 0.0;
  const double rho = 1.0;
  const double gravity_z = -10.0;
  std::vector<double> vec_b;
  std::vector<int> aMSFlag; // master slave flag
  SolveProblem_PlateBendingMITC3(aVal,vec_b,aMSFlag,mat_A,ilu_A,aBCFlag,
                                 thickness,myu,lambda,rho,gravity_z,
                                 aXY0,aTri);
 
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 0.8;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_ZTOP;
  viewer.nav.camera.psi = 0.2;
  viewer.nav.camera.theta = 0.2;
  delfem2::opengl::setSomeLighting();

  {
    assert( fabs(lambda)<1.0e-10 );
    const double E = myu*2.0;
    const double I = thickness*thickness*thickness*leny/12.0;
    const double W = thickness*lenx*leny*rho*gravity_z;
    const double w = W/lenx;
    const double disp = w*(lenx*lenx*lenx*lenx)/(8.0*E*I);
    std::cout << "disp:" << disp << std::endl;
    for(size_t ip=0;ip<aXY0.size()/2;++ip){
      const double px = aXY0[ip*2+0];
      if( fabs(px-(+lenx*0.5)) > 0.0001 ){ continue; }
      std::cout << aVal[ip*3+0] << std::endl;
    }
  }
  
  while(!glfwWindowShouldClose(viewer.window)){
    viewer.DrawBegin_oldGL();
    myGlutDisplay(aXY0, aTri, aVal);
    viewer.DrawEnd_oldGL();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);  
}


