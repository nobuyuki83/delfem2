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
#include "delfem2/points.h"
#include "delfem2/mshio.h"
#include "delfem2/mshuni.h"
#include "delfem2/vec3.h"
#include "delfem2/slice.h"
#include <GLFW/glfw3.h>
#include <vector>
#include <string>
#include <cassert>
#include <cstdlib>
#include <set>

namespace dfm2 = delfem2;

// -------------------------

std::vector<double> aXYZ;
std::vector<unsigned int> aTri;
std::vector<delfem2::CSliceTriMesh> aCS;
std::vector< std::set<unsigned int> > ReebGraphCS;
std::vector<dfm2::CVec3d> aCG_CS;

// ---------------------------

void myGlutDisplay()
{
  ::glEnable(GL_LIGHTING);
  delfem2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aTri);
  
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,0,0);
  ::glLineWidth(5);
  for(auto & cs : aCS){
    ::glBegin(GL_LINE_LOOP);
    for(const auto & seg : cs.aTriInfo){
      double pA[3],pB[3]; seg.Pos3D(pA,pB,
                                    aXYZ,aTri);
      ::glVertex3d(pA[0],pA[1],pA[2]);
    }
    ::glEnd();
  }
  
  ::glDisable(GL_DEPTH_TEST);
  
  for(size_t ics=0;ics<ReebGraphCS.size();++ics){
    ::glColor3d(0,0,0);
    ::glPointSize(10);
    ::glBegin(GL_POINTS);
    ::glVertex3d(aCG_CS[ics].x, aCG_CS[ics].y, aCG_CS[ics].z);
    ::glEnd();
  }
  for(size_t ics=0;ics<ReebGraphCS.size();++ics){
    for(auto itr = ReebGraphCS[ics].begin();itr!=ReebGraphCS[ics].end();++itr){
      const unsigned int jcs = *itr;
      assert( jcs < aCS.size());
      assert( abs(aCS[ics].IndHeight() - aCS[jcs].IndHeight()) == 1 );
      ::glColor3d(0,0,0);
      ::glLineWidth(3);
      ::glBegin(GL_LINES);
      ::glVertex3d(aCG_CS[ics].x, aCG_CS[ics].y, aCG_CS[ics].z);
      ::glVertex3d(aCG_CS[jcs].x, aCG_CS[jcs].y, aCG_CS[jcs].z);
      ::glEnd();
    }
  }
  ::glEnable(GL_DEPTH_TEST);
}

void Hoge(){
  delfem2::Read_Ply(std::string(PATH_INPUT_DIR)+"/bunny_1k.ply",
           aXYZ,aTri);
  delfem2::Normalize_Points3(aXYZ);
  std::vector<unsigned int> aTriSuTri;
  ElSuEl_MeshElem(aTriSuTri,
                  aTri.data(), aTri.size()/3, delfem2::MESHELEM_TRI, aXYZ.size()/3);
  
  
  std::vector<double> aHeight;
  aHeight.push_back(-0.3);
  aHeight.push_back(-0.2);
  aHeight.push_back(-0.1);
  aHeight.push_back(-0.0);
  aHeight.push_back(+0.1);
  aHeight.push_back(+0.2);
  aHeight.push_back(+0.3);
  const double nrm[3] = {0,1,0};
  const double org[3] = {0,0,0};
  {
    std::vector<double> aHeightVtx(aXYZ.size()/3);
    for(size_t ip=0;ip<aXYZ.size()/3;++ip){
      double x0 = aXYZ[ip*3+0] - org[0];
      double y0 = aXYZ[ip*3+1] - org[1];
      double z0 = aXYZ[ip*3+2] - org[2];
      aHeightVtx[ip] = nrm[0]*x0 + nrm[1]*y0 + nrm[2]*z0;
    }
    Slice_MeshTri3D_Heights(aCS,
                            aHeight,
                            aHeightVtx,
                            aTri,aTriSuTri);
  }
  MakeReebGraph(ReebGraphCS,
                aCS, aTri, aTriSuTri);
  assert( aCS.size() == ReebGraphCS.size() );
  
  aCG_CS.resize(aCS.size());
  for(size_t ics=0;ics<aCS.size();++ics){
    const double h0 = aHeight[aCS[ics].IndHeight()];
    const double po[3] = {org[0]+nrm[0]*h0,  org[1]+nrm[1]*h0,  org[2]+nrm[2]*h0 };
    double sum_area = 0.0;
    dfm2::CVec3d cg(0,0,0);
    for(auto & iseg : aCS[ics].aTriInfo){
      double pA[3],pB[3];
      iseg.Pos3D(pA,pB,
                 aXYZ,aTri);
      double n0[3]; dfm2::NormalTri3(n0,
                                      pA,pB,po);
      const double area0 = n0[0]*nrm[0] + n0[1]*nrm[1] + n0[2]*nrm[2];
      sum_area += area0;
      cg.p[0] += area0*(po[0]+pA[0]+pB[0])/3.0;
      cg.p[1] += area0*(po[1]+pA[1]+pB[1])/3.0;
      cg.p[2] += area0*(po[2]+pA[2]+pB[2])/3.0;
    }
    cg /= sum_area;
    aCG_CS[ics] = cg;
  }
}

int main(int argc,char* argv[])
{
  delfem2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  
  delfem2::opengl::setSomeLighting();
  
  viewer.camera.view_height = 0.5;
  viewer.camera.camera_rot_mode  = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  
  Hoge();
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
