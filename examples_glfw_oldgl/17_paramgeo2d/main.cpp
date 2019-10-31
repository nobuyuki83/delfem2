/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <math.h>
#include "delfem2/vec2.h"
#include "delfem2/camera.h"
#include "delfem2/paramgeo_v23.h"

#include <GLFW/glfw3.h>
#include "delfem2/glfw_viewer.hpp"
#include "delfem2/gl2_v23.h"
#include "delfem2/gl2_funcs.h"

// ---------------------------------------------

int ipoint_picked;
double pos_down[2];
double pos_cur[2];

int ndegree = 2;
std::vector<int> aKnotMulti;
std::vector<double> aKnot;
std::vector<double> aKnotFlat;
std::vector<CVector2> aCtrlPoint;

const int nsmpl = 100;
std::vector<CVector2> polyline0; // current test

// -----------------------------------------------

void SetExample(int ndeg, int ncp)
{
  ndegree = ndeg;
//  const int nk = ncp+ndeg+1;
  const int ndiv = ncp-ndeg;
  aKnot.assign(ndiv+1,0);
  for(int idiv=0;idiv<ndiv+1;++idiv){
    aKnot[idiv] = (double)idiv/ndiv;
  }
  aKnotMulti.assign(ndiv+1,1);
  aKnotMulti[   0] = ndeg+1;
  aKnotMulti[ndiv] = ndeg+1;
  FlatKnot(aKnotFlat, aKnotMulti, aKnot);
  for(unsigned int ik=0;ik<aKnotFlat.size();++ik){
    std::cout << "knot" << ik << " " << aKnotFlat[ik] << std::endl;
  }
  //
  aCtrlPoint.assign(ncp, CVector2(0,0));  //7+2+1 = 10
  for(unsigned int i=0;i<aCtrlPoint.size();++i){
    aCtrlPoint[i].x = i*2.0/(aCtrlPoint.size()-1)-1.0;
  }
}

void myGlutDisplay(void)
{
  ::glDisable(GL_LIGHTING);
  
  ::glLineWidth(2);
  ::glPointSize(5);
  ::glColor3d(1,0,0);
  
  ::glColor3d(0,0.5,0);
  drawPolyLine2D(aCtrlPoint);
      
  ::glPointSize(2);
  ::glColor3d(0,0,0);
  drawPolyLine2D(polyline0);
}

int main(int argc,char* argv[])
{
  CViewer_GLFW viewer;
  viewer.Init_GLold();
  
  // -----------
  
  SetExample(3,6);
  SampleBSpline(polyline0, nsmpl, ndegree, aKnotFlat, aCtrlPoint);
  
  // -----------
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    {
      static int iframe = 0;
      if( iframe == 0 ){
        for(int icp=0;icp<aCtrlPoint.size();++icp){
          aCtrlPoint[icp].x += ((double)rand()/RAND_MAX-0.5)*0.1;
          aCtrlPoint[icp].y += ((double)rand()/RAND_MAX-0.5)*0.1;
        }
        SampleBSpline(polyline0, nsmpl, ndegree, aKnotFlat, aCtrlPoint);
      }
      iframe = (iframe+1)%50;
    }
    viewer.DrawBegin_Glold();
    
    myGlutDisplay();
        
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


