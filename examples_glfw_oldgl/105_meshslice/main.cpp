#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <set>
#include <stack>
#include "delfem2/msh.h"
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vec3.h"
#include "delfem2/slice.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/gl2_funcs.h"
#include "delfem2/opengl/gl2_color.h"

// -------------------------

std::vector<double> aXYZ;
std::vector<unsigned int> aTri;
std::vector<CSliceTriMesh> aCS;
std::vector< std::set<unsigned int> > ReebGraphCS;
std::vector<CVector3> aCG_CS;

// ---------------------------

void myGlutDisplay(void)
{
  ::glEnable(GL_LIGHTING);
  delfem2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aTri);
  
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,0,0);
  ::glLineWidth(5);
  for(int iloop=0;iloop<aCS.size();++iloop){
    ::glBegin(GL_LINE_LOOP);
    for(int iseg=0;iseg<aCS[iloop].aTriInfo.size();++iseg){
      const CVector3& pA = aCS[iloop].aTriInfo[iseg].pA;
      ::glVertex3d(pA.x,pA.y,pA.z);
    }
    ::glEnd();
  }
  
  ::glDisable(GL_DEPTH_TEST);
  
  for(int ics=0;ics<ReebGraphCS.size();++ics){
    ::glColor3d(0,0,0);
    ::glPointSize(10);
    ::glBegin(GL_POINTS);
    ::glVertex3d(aCG_CS[ics].x, aCG_CS[ics].y, aCG_CS[ics].z);
    ::glEnd();
  }
  for(int ics=0;ics<ReebGraphCS.size();++ics){
    for(std::set<unsigned int>::iterator itr = ReebGraphCS[ics].begin();itr!=ReebGraphCS[ics].end();++itr){
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
  Read_Ply(std::string(PATH_INPUT_DIR)+"/bunny_1k.ply",
           aXYZ,aTri);
  Normalize(aXYZ);
  std::vector<int> aTriSurRel;
  makeSurroundingRelationship(aTriSurRel,
                              aTri.data(), aTri.size()/3, MESHELEM_TRI, aXYZ.size()/3);
  
  
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
  Slice_MeshTri3D_Heights(aCS,
                          aHeight,
                          nrm, org,
                          aXYZ,aTri,aTriSurRel);
  MakeReebGraph(ReebGraphCS,
                aCS, aTri, aTriSurRel);
  assert( aCS.size() == ReebGraphCS.size() );
  
  aCG_CS.resize(aCS.size());
  for(int ics=0;ics<aCS.size();++ics){
    const double h0 = aHeight[aCS[ics].IndHeight()];
    const double po[3] = {org[0]+nrm[0]*h0,  org[1]+nrm[1]*h0,  org[2]+nrm[2]*h0 };
    double sum_area = 0.0;
    CVector3 cg(0,0,0);
    for(int iseg=0;iseg<aCS[ics].aTriInfo.size();++iseg){
      double n0[3]; NormalTri3D(n0,
                                aCS[ics].aTriInfo[iseg].pA,
                                aCS[ics].aTriInfo[iseg].pB,
                                po);
      double area0 = n0[0]*nrm[0] + n0[1]*nrm[1] + n0[2]*nrm[2];
      sum_area += area0;
      cg.x += area0*(po[0]+aCS[ics].aTriInfo[iseg].pA[0]+aCS[ics].aTriInfo[iseg].pB[0])/3.0;
      cg.y += area0*(po[1]+aCS[ics].aTriInfo[iseg].pA[1]+aCS[ics].aTriInfo[iseg].pB[1])/3.0;
      cg.z += area0*(po[2]+aCS[ics].aTriInfo[iseg].pA[2]+aCS[ics].aTriInfo[iseg].pB[2])/3.0;
    }
    cg /= sum_area;
    aCG_CS[ics] = cg;
  }
}

int main(int argc,char* argv[])
{
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_GLold();
  
  delfem2::opengl::setSomeLighting();
  
  viewer.nav.camera.view_height = 0.5;
  viewer.nav.camera.camera_rot_mode  = delfem2::CAMERA_ROT_TBALL;
  
  Hoge();
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    viewer.DrawBegin_Glold();
    myGlutDisplay();
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
