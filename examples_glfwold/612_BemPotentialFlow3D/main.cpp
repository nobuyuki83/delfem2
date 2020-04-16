/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <iostream>
#include "delfem2/bem.h"
#include "delfem2/primitive.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/v3q_glold.h"

namespace dfm2 = delfem2;

// ------------------------------------

dfm2::CVec3d posGrid(int ig, int jg, double lenGrid, int ngrid){
  dfm2::CVec3d p;
  p.p[0] = ig*lenGrid/ngrid-lenGrid*0.5;
  p.p[1] = jg*lenGrid/ngrid-lenGrid*0.5;
  p.p[2] = 0;
  return p;
}

// ------------------------------------

dfm2::CVec3d Velo;

std::vector<unsigned int> aTri;
std::vector<double> aXYZ;

std::vector<double> aSol;
double min_sol,max_sol;

double lenGrid = 1.0;
int ngrid = 45;
std::vector<dfm2::CVec3d> aValGrid;

int imode_draw = 0;

// ------------------------------------------

void SetProblem()
{
  dfm2::MeshTri3D_Cube(aXYZ, aTri, 10);
  for(int i=0;i<aXYZ.size();++i){ aXYZ[i] *= 0.3; }
  Velo.SetZero();
  Velo[0] = 1.0;
  {
    std::vector<double> A, f;
    makeLinearSystem_PotentialFlow_Order0th(A,f,
                                            //
                                            Velo,
                                            1,
                                            aXYZ,
                                            aTri);
    aSol.assign(aTri.size()/3,0.0);
    double conv_ratio = 1.0e-8;
    int iteration = 1000;
    dfm2::Solve_BiCGSTAB(conv_ratio, iteration, aSol,
                   A, f);
    std::cout << conv_ratio << " " << iteration << std::endl;
    min_sol = max_sol = aSol[0];
    for (int itri = 0; itri<aTri.size()/3; itri++){
      double v = aSol[itri];
      min_sol = (v<min_sol) ? v : min_sol;
      max_sol = (v>max_sol) ? v : max_sol;
    }
    std::cout << min_sol << " " << max_sol << std::endl;
  }
  // ----------
  aValGrid.resize( (ngrid+1)*(ngrid+1) );
  for(int ig=0;ig<ngrid+1;ig++){
    for(int jg=0;jg<ngrid+1;jg++){
      dfm2::CVec3d p = posGrid(ig,jg,lenGrid,ngrid);
      dfm2::CVec3d gradphi_pos;
      double phi_pos;
      evaluateField_PotentialFlow_Order0th(phi_pos,gradphi_pos,
                                           p, Velo,
                                           1,
                                           aSol,aXYZ,aTri);
      aValGrid[ ig*(ngrid+1)+jg ] = gradphi_pos;
    }
  }
}

// -------------------------------------------

// -------------------------------------------


void myGlutDisplay(void)
{
  {
//    Draw_SurfaceMeshFaceEdge(aXYZ, aTri);
    {
      const int nTri = (int)aTri.size()/3;
      //
      ::glBegin(GL_TRIANGLES);
      for (int itri = 0; itri<nTri; itri++){
        {
//          double v = (aSol[itri]-min_sol)/(max_sol-min_sol);
          double v = aSol[itri]+0.5;
//          double c[4]; dfm2::opengl::heatmap(v, c);
          double c[4] = {1,0,0,0};
          ::glColor3dv(c);
        }
        const int i1 = aTri[itri*3+0];
        const int i2 = aTri[itri*3+1];
        const int i3 = aTri[itri*3+2];
        dfm2::opengl::myGlVertex3(i1, aXYZ);
        dfm2::opengl::myGlVertex3(i2, aXYZ);
        dfm2::opengl::myGlVertex3(i3, aXYZ);
      }
      ::glEnd();
      // ---------
      ::glDisable(GL_LIGHTING);
      ::glColor3d(0, 0, 0);
      ::glBegin(GL_LINES);
      for (int itri = 0; itri<nTri; itri++){
        const unsigned int i1 = aTri[itri*3+0];
        const unsigned int i2 = aTri[itri*3+1];
        const unsigned int i3 = aTri[itri*3+2];
        dfm2::opengl::myGlVertex3(i1, aXYZ);
        dfm2::opengl::myGlVertex3(i2, aXYZ);
        dfm2::opengl::myGlVertex3(i2, aXYZ);
        dfm2::opengl::myGlVertex3(i3, aXYZ);
        dfm2::opengl::myGlVertex3(i3, aXYZ);
        dfm2::opengl::myGlVertex3(i1, aXYZ);
      }
      ::glEnd();
    }
  }
  
  {
    for(int ig=0;ig<ngrid+1;ig++){
    for(int jg=0;jg<ngrid+1;jg++){
      dfm2::CVec3d p = posGrid(ig, jg, lenGrid, ngrid);
      dfm2::CVec3d val0 = aValGrid[ig*(ngrid+1)+jg];
      ::glColor3d(1, 0, 0);
      dfm2::opengl::DrawArrow(p, val0*0.03);
      //
      /*
      CVector3 val1 = Velo;
      ::glColor3d(0, 0, 1);
      dfm2::opengl::DrawArrow(p, val1*0.03);
       */
    }
    }
  }
}

/*
void myGlutKeyboard(unsigned char Key, int x, int y)
{
	switch(Key)
	{
    case 'q':
    case 'Q':
    case '\033':
      exit(0);
    case 'a':
      is_animation = !is_animation;
      break;
    case 'd': // change draw mode
      imode_draw = (imode_draw+1) % 3;
      if( imode_draw != 0 ){ setValueField(); }
      break;
    case 'k':
    {
      for (int i = 0; i<6;i++){
        Velo.SetZero();
        Omg.SetZero();
        if (i==0){ Velo[0] = 1.0; }
        if (i==1){ Velo[1] = 1.0; }
        if (i==2){ Velo[2] = 1.0; }
        if (i==3){ Omg[0] = 1.0; }
        if (i==4){ Omg[1] = 1.0; }
        if (i==5){ Omg[2] = 1.0; }
        Solve();
        CVector3 Vj(0, 0, 0), Oj(0, 0, 0);
        for (int itri = 0; itri<aTri.size()/3; itri++){
          CVector3 pi = MidPoint(itri, aTri, aXYZ);
          CVector3 ni = NormalTri(itri, aTri, aXYZ);
          double areai = ni.Length()*0.5;
          ni.SetNormalizedVector();
          double phii = 0;
          for (int jtri = 0; jtri<aTri.size()/3; jtri++){
            CVector3 pj = MidPoint(jtri, aTri, aXYZ);
            CVector3 nj = NormalTri(jtri, aTri, aXYZ);
            double areaj = nj.Length()*0.5;
            if (jtri==itri){
              phii += aSol[jtri]*0.5;
            }
            else{
              CVector3 v = pi-pj;
              double vlen = v.Length();
              phii += areaj*aSol[jtri]/vlen;
            }
          }
//          std::cout<<itri<<" "<<aSol[itri]<<" "<<phii<<std::endl;
          Vj += areai*phii*ni;
          Oj += areai*phii*(pi^ni);
        }
        std::cout<<i<<"  "<<Vj<<"  " << Oj<<std::endl;
      }
      Velo.SetZero();
      Omg.SetZero();
      Solve();
      break;
    }
    case ' ':
    {
      static int iprob = 0;
      iprob++;
      if( iprob >= 6 ){ iprob=0; }
      std::cout<<iprob<<std::endl;
      Velo.SetZero();
      Omg.SetZero();
      if (iprob==0){ Velo[0] = 1.0; }
      if (iprob==1){ Velo[1] = 1.0; }
      if (iprob==2){ Velo[2] = 1.0; }
      if (iprob==3){ Omg[0] = 1.0; }
      if (iprob==4){ Omg[1] = 1.0; }
      if (iprob==5){ Omg[2] = 1.0; }
      Solve();
      setValueField();
    }
  }
	::glutPostRedisplay();
}
        */

int main(int argc,char* argv[])
{
  SetProblem();
  
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  delfem2::opengl::setSomeLighting();
  while (!glfwWindowShouldClose(viewer.window))
  {
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


