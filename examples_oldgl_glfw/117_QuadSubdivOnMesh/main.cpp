/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/srchuni_v3.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshsubdiv.h"
#include "delfem2/vec3.h"
#include "delfem2/mshprimitive.h"
#include <GLFW/glfw3.h>
#include <vector>
#include <string>
#include <cstdlib>

namespace dfm2 = delfem2;

// ---------------------------------------------------

std::vector<double> aXYZ;
std::vector<unsigned int> aTri;

const unsigned int nsubdiv = 3;
std::vector<double> aXYZ_Quad;
std::vector<double> aNorm_Quad;
std::vector< std::vector<unsigned int>> aaQuad;


// ---------------------------------------------------

void InitializeProblem() {
    
  dfm2::MeshTri3D_Sphere(
      aXYZ, aTri,
      1.0, 12, 34);
  {
    const double bbmin[3] = {-1,-1,-1};
    const double bbmax[3] = {+1,+1,+1};
    aaQuad.resize(1);
    dfm2::MeshQuad3_CubeVox(aXYZ_Quad, aaQuad[0],
        bbmin, bbmax);
  }
  
  for(unsigned int ip=0;ip<aXYZ_Quad.size()/3;++ip){
    dfm2::CVec3d p0(aXYZ_Quad[ip*3+0], aXYZ_Quad[ip*3+1], aXYZ_Quad[ip*3+2]);
    dfm2::CPtElm2<double> pes0 = Nearest_Point_MeshTri3D(p0,
        aXYZ, aTri);
    dfm2::CVec3d q0 = pes0.Pos_Tri(aXYZ, aTri);
    aXYZ_Quad[ip*3+0] = q0.x;
    aXYZ_Quad[ip*3+1] = q0.y;
    aXYZ_Quad[ip*3+2] = q0.z;
  }
  
  dfm2::Normal_MeshQuad3(aNorm_Quad,
                         aXYZ_Quad, aaQuad[0]);
  
  aaQuad.resize(nsubdiv+1);
  for(unsigned int isubdiv=0;isubdiv<nsubdiv;++isubdiv){
    std::vector<unsigned int> psup_ind, psup;
    std::vector<unsigned int> aEdgeFace0;
    dfm2::SubdivTopo_MeshQuad(
        aaQuad[isubdiv+1], psup_ind, psup, aEdgeFace0,
        aaQuad[isubdiv].data(), aaQuad[isubdiv].size()/4,
        aXYZ_Quad.size()/3);
    const unsigned int nv0 = aXYZ_Quad.size()/3;
    const unsigned int ne0 = aEdgeFace0.size()/4;
    const unsigned int nq0 = aaQuad[isubdiv].size()/4;
    aXYZ_Quad.resize((nv0+ne0+nq0)*3);
    aNorm_Quad.resize((nv0+ne0+nq0)*3);
    { // enter dangerous zone
      const std::vector<unsigned int>& aQuad0 = aaQuad[isubdiv];
      const std::vector<double>& aXYZ0 = aXYZ_Quad;
      std::vector<double>& aXYZ1 = aXYZ_Quad;
      for(unsigned int ie=0;ie<ne0;++ie){
        const unsigned int iv0 = aEdgeFace0[ie*4+0];
        const unsigned int iv1 = aEdgeFace0[ie*4+1];
        dfm2::AverageTwo3(aXYZ1.data()+(nv0+ie)*3,
                          aXYZ0.data()+iv0*3, aXYZ0.data()+iv1*3);
        dfm2::AverageTwo3(aNorm_Quad.data()+(nv0+ie)*3,
                          aNorm_Quad.data()+iv0*3, aNorm_Quad.data()+iv1*3);
        dfm2::Normalize3(aNorm_Quad.data()+(nv0+ie)*3);
      }
      for(unsigned int iq=0;iq<nq0;++iq){
        const unsigned int iv0 = aQuad0[iq*4+0];
        const unsigned int iv1 = aQuad0[iq*4+1];
        const unsigned int iv2 = aQuad0[iq*4+2];
        const unsigned int iv3 = aQuad0[iq*4+3];
        dfm2::AverageFour3(aXYZ1.data()+(nv0+ne0+iq)*3,
                           aXYZ0.data()+iv0*3, aXYZ0.data()+iv1*3, aXYZ0.data()+iv2*3, aXYZ0.data()+iv3*3);
        dfm2::AverageFour3(aNorm_Quad.data()+(nv0+ne0+iq)*3,
                           aNorm_Quad.data()+iv0*3, aNorm_Quad.data()+iv1*3, aNorm_Quad.data()+iv2*3, aNorm_Quad.data()+iv3*3);
        dfm2::Normalize3(aNorm_Quad.data()+(nv0+ne0+iq)*3);
      }
    }
    for(unsigned int ip=nv0;ip<nv0+ne0+nq0;++ip){
      const dfm2::CVec3d p0( aXYZ_Quad[ip*3+0], aXYZ_Quad[ip*3+1], aXYZ_Quad[ip*3+2] );
      const dfm2::CVec3d n0( aNorm_Quad[ip*3+0], aNorm_Quad[ip*3+1], aNorm_Quad[ip*3+2] );
      std::vector<dfm2::CPtElm2<double>> aPES = IntersectionLine_MeshTri3(
          p0, n0,
          aTri, aXYZ,
          1.0e-3);
      std::map<double, dfm2::CPtElm2<double>> mapPES;
      for(auto & pes : aPES){
        dfm2::CVec3d q0 = pes.Pos_Tri(aXYZ, aTri);
        double h = (q0-p0)*n0;
        mapPES.insert( std::make_pair(-h, pes) );
      }
      if( !aPES.empty() ){
        dfm2::CPtElm2<double>& pes0 = mapPES.begin()->second;
        dfm2::CVec3d q0 = pes0.Pos_Tri(aXYZ, aTri);
        aXYZ_Quad[ip*3+0] = q0.x;
        aXYZ_Quad[ip*3+1] = q0.y;
        aXYZ_Quad[ip*3+2] = q0.z;
      }
    }
    { // make normal for new vertices
      aNorm_Quad.resize((nv0+ne0+nq0)*3);
      std::vector<double> aNorm1;
      dfm2::Normal_MeshQuad3(aNorm1,
                             aXYZ_Quad, aaQuad[isubdiv+1]);
      for(unsigned int ip=nv0;ip<nv0+ne0+nq0;++ip){
        aNorm_Quad[ip*3+0] = aNorm1[ip*3+0];
        aNorm_Quad[ip*3+1] = aNorm1[ip*3+1];
        aNorm_Quad[ip*3+2] = aNorm1[ip*3+2];
      }
    }
    
  }
}

void UpdateProblem() {
}

// -----------------------------------------------------

void myGlutDisplay()
{
  ::glEnable(GL_LIGHTING);
  { //
    ::glDisable(GL_TEXTURE_2D);
    float gray[4] = {1.0f,0.0f,0.0f,1.f};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, gray);
    float shine[4] = {0,0,0,0};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
    ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 127.0);
    dfm2::opengl::DrawMeshQuad3D_FaceNorm(aXYZ_Quad,aaQuad[nsubdiv]);
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0.0, 0.0, 0.0);
    dfm2::opengl::DrawMeshQuad3D_Edge(aXYZ_Quad,aaQuad[nsubdiv]);
    /*
    ::glDisable(GL_LIGHTING);
    ::glColor3d(1.0, 1.0, 1.0);
    dfm2::opengl::DrawPoints3d_NormVtx(aXYZ_Quad, aNorm_Quad, 0.1);
     */
  }
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1.0, 1.0, 1.0);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ, aTri);
  
  ::glEnd();
  ::glEnable(GL_DEPTH_TEST);
}

int main(int argc,char* argv[])
{
  // ----------------

  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.view_height = 1.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  dfm2::opengl::setSomeLighting();

  InitializeProblem();
  UpdateProblem();

  int iframe = 0;
  while (!glfwWindowShouldClose(viewer.window))
  {
    iframe += 1;
    UpdateProblem();
    // -------
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.SwapBuffers();
    glfwPollEvents();
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
