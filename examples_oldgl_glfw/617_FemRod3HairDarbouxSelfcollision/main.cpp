/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/ls_solver_block_sparse.h"
#include "delfem2/ls_block_sparse.h"
#include "delfem2/hair_darboux_util.h"
#include "delfem2/hair_darboux_collision.h"
#include "delfem2/geo_ccd.h"
#include "delfem2/geo_line.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/srch_bv3_sphere.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/funcs.h"

namespace dfm2 = delfem2;

// -------------------------------------

void myDraw(
	const std::vector<dfm2::CVec3d>& aP,
	const std::vector<dfm2::CVec3d>& aS,
	std::vector<unsigned int>& aIP_HairRoot)
{
  assert(!aIP_HairRoot.empty());
  const unsigned int nhair = static_cast<unsigned int>(aIP_HairRoot.size())-1;
  for(unsigned int ihair=0;ihair<nhair;++ihair){
    const unsigned int ips = aIP_HairRoot[ihair];
    const unsigned int ipe = aIP_HairRoot[ihair+1];
    assert( aP.size() == aS.size() );
    ::glDisable(GL_LIGHTING);
    ::glColor3d(1,0,0);
    ::glPointSize(3);
    ::glBegin(GL_POINTS);
    for(unsigned int ip=ips;ip<ipe;++ip){
      ::glVertex3d(aP[ip].x, aP[ip].y, aP[ip].z);
    }
    ::glEnd();
    // ------------
    ::glColor3d(0,0,0);
    ::glLineWidth(3);
    ::glBegin(GL_LINES);
    unsigned int ns = ipe-ips-1;
    for(unsigned int is=0;is<ns;++is){
      const unsigned int ip0 = ips+is+0; assert( ip0 < aP.size() );
      const unsigned int ip1 = ips+is+1; assert( ip1 < aP.size() );
      ::glVertex3d(aP[ip0].x, aP[ip0].y, aP[ip0].z);
      ::glVertex3d(aP[ip1].x, aP[ip1].y, aP[ip1].z);
    }
    ::glEnd();
    // --------------
    ::glBegin(GL_LINES);
    for(unsigned int is=0;is<ns;++is){
      const unsigned int ip0 = ips+is+0; assert( ip0 < aP.size() );
      const unsigned int ip1 = ips+is+1; assert( ip1 < aP.size() );
      dfm2::CVec3d p01 = 0.5*(aP[ip0]+aP[ip1]);
      double l01 = (aP[ip0]-aP[ip1]).norm();
      dfm2::opengl::myGlVertex(p01);
      dfm2::opengl::myGlVertex(p01+(l01*0.5)*aS[is]);
    }
    ::glEnd();
  }

  { // draw rectangle wall
    ::glEnable(GL_LIGHTING);
    ::glBegin(GL_QUADS);
    ::glNormal3d(1.0, 1.2, 0);
    ::glVertex3d(-0.8,-1,-2);
    ::glVertex3d(-1.4,+2,-2);
    ::glVertex3d(-1.4,+2,+2);
    ::glVertex3d(-0.8,-1,+2);
    ::glEnd();
  }
}

void FindRodHairContactCCD(
    std::vector<dfm2::CContactHair>& aCollision,
    const double clearance,
    const std::vector<dfm2::CVec3d>& aP,
    const std::vector<unsigned int>& aIP_HairRoot,
    const std::vector<dfm2::CVec3d>& aPt) // indexes of root point
{
  const auto nr = static_cast<unsigned int>(aP.size());
  assert( aPt.size() == nr );
  std::vector<bool> aIsRod(nr,true);
  for(unsigned int ihair=0;ihair<aIP_HairRoot.size()-1;++ihair){
    unsigned int ip0 = aIP_HairRoot[ihair];
    unsigned int ip1 = aIP_HairRoot[ihair+1]-1;
    aIsRod[ip0] = false;
    aIsRod[ip1] = false;
  }
  std::vector<dfm2::CBV3_Sphere<double>> aBV(nr);
  for(unsigned int ir=0;ir<nr;++ir) {
    if (!aIsRod[ir]) {
      aBV[ir].Set_Inactive();
      continue;
    }
    const dfm2::CVec3d p0s = aP[ir];
    const dfm2::CVec3d p1s = aP[ir + 1];
    const dfm2::CVec3d p0e = aPt[ir];
    const dfm2::CVec3d p1e = aPt[ir + 1];
    aBV[ir].SetPoints4(p0s.p,p1s.p,p0e.p,p1e.p, clearance*0.5);
  }

  for(unsigned int ir=0;ir<nr;++ir){
    if( !aIsRod[ir] ){ continue; }
    const dfm2::CVec3d p0s = aP[ir+0];
    const dfm2::CVec3d p1s = aP[ir+1];
    const dfm2::CVec3d p0e = aPt[ir+0];
    const dfm2::CVec3d p1e = aPt[ir+1];
    for(unsigned int jr=ir+2;jr<nr;++jr) {
//    for(unsigned int jr=0;jr<nr;++jr) {
      if( !aIsRod[jr] ){ continue; }
//      if( jr == ir || jr == ir+1 || jr == jr-1 ){ continue; } // neighbouring rod is always intersecting
      if( !aBV[ir].IsIntersect(aBV[jr]) ) continue;
      // -------
      const dfm2::CVec3d q0s = aP[jr+0];
      const dfm2::CVec3d q1s = aP[jr+1];
      const dfm2::CVec3d q0e = aPt[jr+0];
      const dfm2::CVec3d q1e = aPt[jr+1];
      double p[3] = {0.5, 0.5, 0.5};
      {
        double D;
        dfm2::CVec3d Da,Db;
        dfm2::Nearest_Line3_Line3(
          D, Da, Db,
          p0s, p1s-p0s, q0s, q1s-q0s);
        Da /= D;
        Db /= D;
        p[0] = (Da-p0s).dot(p1s-p0s)/(p1s-p0s).squaredNorm();
        p[1] = (Db-q0s).dot(q1s-q0s)/(q1s-q0s).squaredNorm();
//        std::cout << ir << " " << jr << " --> " << p[0] << " " << p[1] << " " << (Da-Db)*(p1s-p0s) << " " << (Da-Db)*(q1s-q0s) << " " << (Da-Db).Length() << std::endl;
        if( p[0] > 1 ){ p[0] = 1; } else if( p[0] < 0 ){ p[0] = 0; }
        if( p[1] > 1 ){ p[1] = 1; } else if( p[1] < 0 ){ p[1] = 0; }
        p[2] = 0.0;
      }
      double ps0 = p[0];
      double pt0 = p[1];
      // space time collision
      {
        double len0 = dfm2::Nearest_LineSeg_LineSeg_CCD_Iteration(
            p,
            p0s, p0e, p1s, p1e, q0s, q0e, q1s, q1e, 10);
        bool is_near = false;
        if( len0 < clearance ){ is_near = true; }
        if( !is_near ){ continue; }
      }
//      double s1 = p[0], t1 = p[1];
      double s1 = ps0, t1 = pt0;
      dfm2::CVec3d vs = (1-s1)*p0s + s1*p1s - (1-t1)*q0s - t1*q1s; // difference of positions at the begining of a time step
      /*
      double u1 = p[2];
      dfm2::CVec3d vm =
          + ((1 - s1) * (1 - u1)) * p0s + ((1 - s1) * u1) * p0e + (s1 * (1 - u1)) * p1s + (s1 * u1) * p1e
          - ((1 - t1) * (1 - u1)) * q0s - ((1 - t1) * u1) * q0e - (t1 * (1 - u1)) * q1s - (t1 * u1) * q1e;
      std::cout << vm.Length () << " " << clearance << std::endl;
       */
      dfm2::CContactHair ch{ir + 0, ir + 1, s1,
                            jr + 0, jr + 1, t1,
                            vs.normalized()};
      aCollision.push_back(ch);
    } // jr
  } // ir
}

int main()
{
  double dt = 0.02;
  double mass = 1.0e-1;
  dfm2::CVec3d gravity(0,-10,0);
  const double stiff_stretch = 1000;
  const double stiff_bendtwist[3] = { 300, 300, 300 };
  const double stiff_contact = 1.0e+4;
  const double clearance = 0.02;
  //
  class CViewerDemo : public dfm2::glfw::CViewer3 {
  public:
    CViewerDemo() : CViewer3(1.5) {}
    void  key_press(int key, 
		[[maybe_unused]] int mods) override {
      if( key == GLFW_KEY_A ){
        is_animation = !is_animation;
      }
      if( key == GLFW_KEY_S ){
        is_animation = true;
        is_step = true;
      }
    }
  public:
    bool is_animation = true;
    bool is_step = false;
  } viewer;
  //
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  // -------
  std::vector<dfm2::CVec3d> aP0, aS0; // initial director vector
  std::vector<unsigned int> aIP_HairRoot; // indexes of root point
  { // make the un-deformed shape of hair
    dfm2::HairDarbouxShape hs{
      30, 1, 0.2, 0.1};
    std::vector<dfm2::CVec3d> ap, as;
    hs.MakeConfigDaboux(ap, as);
    for(int ihair=0;ihair<2;++ihair){
      std::vector<dfm2::CVec3d> ap1 = ap;
      for(auto& p : ap1){
        p.y += ihair;
        p.x += -0.2*ihair-1;
      }
      aIP_HairRoot.push_back(aP0.size());
      aP0.insert(aP0.end(), ap1.begin(), ap1.end());
      aS0.insert(aS0.end(), as.begin(), as.end());
    }
    aIP_HairRoot.push_back(aP0.size());
    assert( aS0.size() == aP0.size() );
  }
  for(int itr=0;itr<10;++itr){ // relax director vectors
    ParallelTransport_RodHair(aP0, aS0, aIP_HairRoot);
  }
  dfm2::LinearSystemSolver_BlockSparse ls_solver; // sparse matrix
  {
    std::vector<unsigned int> psup_ind, psup;
    delfem2::JArray_PSuP_Hair(
        psup_ind,psup,
        aIP_HairRoot);
    ls_solver.Initialize(
        aP0.size(),4,
        psup_ind, psup);
    dfm2::MakeBCFlag_RodHair( // set fixed boundary condition
        ls_solver.dof_bcflag,
        aIP_HairRoot);
    assert(ls_solver.dof_bcflag.size() == aP0.size()*4 );
  }
  // -----------------
  std::vector<dfm2::CVec3d> aP = aP0, aS = aS0;
  std::vector<dfm2::CVec3d> aPV (aP0.size(), dfm2::CVec3d(0,0,0)); // velocity
  std::vector<dfm2::CVec3d> aPt = aP; // temporally positions
  std::vector<dfm2::CContactHair> aContact; // collision in the previous time-step
  double time_cur = 0.0;
  //
//  dfm2::CContactHair ch0 = {2,3,0.5, 40,41,0.5, dfm2::CVec3d(0,0,1)};
  while (true)
  {
    if( viewer.is_animation ){
      time_cur += dt;
      { // set fixed boundary condition
        unsigned int ip0 = aIP_HairRoot[1];
        double z0 = 0.5*sin(0.5*time_cur+0.5);
        aP[ip0].p[2] = aP[ip0+1].p[2] = z0;
        aPt[ip0].p[2] = aPt[ip0+1].p[2] = z0;
      }
      for(unsigned int ip=0;ip<aP.size();++ip){
        if( ls_solver.dof_bcflag[ip*4+0] !=0 ) { continue; } // this is not fixed boundary
        aPt[ip] = aP[ip] + dt * aPV[ip] + dt * dt * gravity;
      }
      const std::vector<dfm2::CVec3d> aPt0 = aPt; // initial temporal positions
      // -----------
      aContact.clear();
      for(int itr=0;itr<1;++itr){
        FindRodHairContactCCD(
            aContact,
            clearance,
            aP, aIP_HairRoot,aPt);
        dfm2::MakeDirectorOrthogonal_RodHair(aS,aPt);
        Solve_RodHairContact(
            aPt, aS, ls_solver,
            stiff_stretch, stiff_bendtwist, mass/(dt*dt),
            aPt0, aP0, aS0, aIP_HairRoot,
            clearance, stiff_contact, aContact);
      }
      // --------------
      for(unsigned int ip=0;ip<aP.size();++ip){
        if( ls_solver.dof_bcflag[ip*4+0] != 0 ){ continue; }
        aPV[ip] = (aPt[ip] - aP[ip])/dt;
        aP[ip] = aPt[ip];
      }
      if( viewer.is_step ){
        viewer.is_step = false;
        viewer.is_animation = false;
      }
    }
    // -------------------------------------
    viewer.DrawBegin_oldGL();
    myDraw(aP, aS, aIP_HairRoot);
    for(const auto& ch : aContact ){
      double s = ch.s;
      double t = ch.t;
      const dfm2::CVec3d p0 = aP[ch.ip0];
      const dfm2::CVec3d p1 = aP[ch.ip1];
      const dfm2::CVec3d q0 = aP[ch.iq0];
      const dfm2::CVec3d q1 = aP[ch.iq1];
      const dfm2::CVec3d p = (1-s)*p0 + s*p1;
      const dfm2::CVec3d q = (1-t)*q0 + t*q1;
      ::glDisable(GL_LIGHTING);
      ::glColor3d(1,0,1);
      ::glLineWidth(4);
      ::glBegin(GL_LINES);
      ::glVertex3dv(p.p);
      ::glVertex3dv(q.p);
      ::glEnd();
    }
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
