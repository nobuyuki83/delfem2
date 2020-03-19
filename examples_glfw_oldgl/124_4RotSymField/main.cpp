/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief demo for view navigation in 3D
 * @details this demo is for showing CViewer_GLFW funtinoalities
 */

#include <cstdlib>
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vec3.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_v23.h"

namespace dfm2 = delfem2;

// ---------------

void FindNearestOrientation(dfm2::CVec3d& d0,
                            dfm2::CVec3d& d1,
                            const dfm2::CVec3d& o0,
                            const dfm2::CVec3d& p0,
                            const dfm2::CVec3d& o1,
                            const dfm2::CVec3d& p1)
{
  const dfm2::CVec3d ad0[4] = { o0, p0, -o0, -p0 };
  const dfm2::CVec3d ad1[2] = { o1, p1 };
  double dot_max = -2;
  for(int j=0;j<2;++j){
    for(int i=0;i<4;++i){
      const double dot = ad0[i]*ad1[j];
      if( dot < dot_max ){ continue; }
      d0 = ad0[i];
      d1 = ad1[j];
      dot_max = dot;
    }
  }
}

void Smooth4RotSym
(std::vector<double>& aOdir,
 const std::vector<double>& aNorm,
 const std::vector<unsigned int>& psup_ind,
 const std::vector<unsigned int>& psup)
{
  for(int iip=0;iip<aOdir.size()/3;++iip){
    const unsigned int ip0 = iip;
    assert( ip0 < psup_ind.size() );
    const unsigned int npj = psup_ind[ip0+1] - psup_ind[ip0+0];
    const dfm2::CVec3d n0 = dfm2::CVec3d(aNorm.data()+ip0*3);
    dfm2::CVec3d o_new = dfm2::CVec3d(aOdir.data()+ip0*3);
    double weight = 0.0;
    for(int jjp=0;jjp<npj;++jjp){
      unsigned int jp1 = psup[psup_ind[ip0]+jjp];
      const dfm2::CVec3d n1 = dfm2::CVec3d(aNorm.data()+jp1*3);
      const dfm2::CVec3d o1 = dfm2::CVec3d(aOdir.data()+jp1*3);
      dfm2::CVec3d d0, d1;
      FindNearestOrientation(d0,d1,
                             o_new,n0^o_new, o1,n1^o1);
      o_new = d0 * weight + d1;
      o_new = (o_new - (o_new*n0)*n0).Normalize();
      weight += 1.0;
    }
    o_new.CopyValueTo(aOdir.data()+ip0*3);
  }
}

void Smooth4RotSym_RandomPermutation
 (std::vector<double>& aOdir,
  std::vector<unsigned int>& permutation0,
  const std::vector<double>& aNorm,
  const std::vector<unsigned int>& psup_ind,
  const std::vector<unsigned int>& psup)
{
  std::vector<unsigned int> permutation1(10);
  std::random_shuffle( permutation0.begin(), permutation0.end() );
  for(int iip=0;iip<aOdir.size()/3;++iip){
    const unsigned int ip0 = permutation0[iip];
    assert( ip0 < psup_ind.size() );
    const unsigned int npj = psup_ind[ip0+1] - psup_ind[ip0+0];
    const dfm2::CVec3d n0 = dfm2::CVec3d(aNorm.data()+ip0*3);
    dfm2::CVec3d o_new = dfm2::CVec3d(aOdir.data()+ip0*3);
    double weight = 0.0;
    {
      permutation1.resize(npj);
      for(int jjp=0;jjp<npj;++jjp){ permutation1[jjp] = jjp; }
      std::random_shuffle( permutation1.begin(), permutation1.end() );
    }
    for(int jjp=0;jjp<npj;++jjp){
      unsigned int jp1 = psup[psup_ind[ip0]+permutation1[jjp]];
      const dfm2::CVec3d n1 = dfm2::CVec3d(aNorm.data()+jp1*3);
      const dfm2::CVec3d o1 = dfm2::CVec3d(aOdir.data()+jp1*3);
      dfm2::CVec3d d0, d1;
      FindNearestOrientation(d0,d1,
                             o_new,n0^o_new, o1,n1^o1);
      o_new = d0 * weight + d1;
      o_new = (o_new - (o_new*n0)*n0).Normalize();
      weight += 1.0;
    }
    o_new.CopyValueTo(aOdir.data()+ip0*3);
  }
}


void InitializeTangentField
(std::vector<double>& aOdir,
 const std::vector<double>& aNorm)
{
  unsigned int np = aNorm.size()/3;
  aOdir.resize(np*3);
  for(int ip=0;ip<np;++ip){
    dfm2::CVec3d o = dfm2::CVec3d::Random();
    dfm2::CVec3d n = dfm2::CVec3d(aNorm.data()+ip*3).Normalize();
    o = (o - (o*n)*n).Normalize();
    o.CopyValueTo(aOdir.data()+ip*3);
  }
}


// ---------------

int main()
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  delfem2::Read_Ply(std::string(PATH_INPUT_DIR)+"/bunny_1k.ply",
                    aXYZ,aTri);
  delfem2::Normalize_Points3(aXYZ);
  std::vector<double> aNorm(aXYZ.size());
  dfm2::Normal_MeshTri3D(aNorm.data(),
                         aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
  
  std::vector<double> aOdir;
  InitializeTangentField(aOdir,aNorm);
  
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                             aTri.data(), aTri.size()/3, 3, aXYZ.size()/3);
  
  std::vector<unsigned int> permutation0(aXYZ.size()/3);
  for(unsigned int i=0;i<aXYZ.size()/3;++i){ permutation0[i] = i; }
  
  // ------------------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.camera_rot_mode = dfm2::CAMERA_ROT_TBALL;
  dfm2::opengl::setSomeLighting();
  unsigned int iframe = 0;
  while (true)
  {
    if( iframe == 0 ){
        InitializeTangentField(aOdir,aNorm);
    }
    if( iframe > 30 ){
      Smooth4RotSym(aOdir,
                    aNorm, psup_ind, psup);
    }
    // --------------------
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ.data(), aTri.data(), aTri.size()/3);
    {
      ::glDisable(GL_LIGHTING);
      double len = 0.03;
      ::glLineWidth(3);
      unsigned int np = aXYZ.size()/3;
      for(int ip=0;ip<np;++ip){
        const dfm2::CVec3d p = dfm2::CVec3d(aXYZ.data()+ip*3);
        const dfm2::CVec3d n = dfm2::CVec3d(aNorm.data()+ip*3).Normalize();
        const dfm2::CVec3d o = dfm2::CVec3d(aOdir.data()+ip*3).Normalize();
        const dfm2::CVec3d q = dfm2::Cross(n,o);
        ::glBegin(GL_LINES);
        ::glColor3d(0,0,0);
        dfm2::opengl::myGlVertex(p);
        dfm2::opengl::myGlVertex(p+len*n);
        ::glColor3d(1,0,0);
        dfm2::opengl::myGlVertex(p-len*o);
        dfm2::opengl::myGlVertex(p+len*o);
        dfm2::opengl::myGlVertex(p-len*q);
        dfm2::opengl::myGlVertex(p+len*q);
        ::glEnd();
      }
    }
    iframe = (iframe+1)%60;
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
