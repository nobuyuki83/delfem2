/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/geoplygn2_v2.h"
#include "delfem2/vec2.h"
#include "delfem2/mat3.h"
#include <GLFW/glfw3.h>

namespace dfm2 = delfem2;

// ------------------------------------------------------

dfm2::CMat3d Mat3_Affine(
    const dfm2::CVec2d& posl,
    const double theta,
    const dfm2::CVec2d& posg)
{
  dfm2::CMat3d mT0; dfm2::Mat3_AffineTranslation(mT0.mat,(-posl).p);
  dfm2::CMat3d mR; dfm2::Mat3_AffineRotation(mR.mat,theta);
  dfm2::CMat3d mT1; dfm2::Mat3_AffineTranslation(mT1.mat,posg.p);
  return mT1*mR*mT0;
}

class CRigidState2
{
public:
  bool is_fix;
  std::vector<dfm2::CVec2d> shape;
  // below: derived from shape
  dfm2::CVec2d posl;
  double mass;
  double I;
  //
  // below: change over time
  dfm2::CVec2d posg;
  double theta;
  dfm2::CVec2d velo;
  double omega;
  // below: temp data for PBD
  dfm2::CVec2d posg_tmp;
  double theta_tmp;
};

void Draw(
    const CRigidState2& rs)
{
  const dfm2::CMat3d mT1RT0 = Mat3_Affine(rs.posl,rs.theta,rs.posg);
  ::glBegin(GL_LINES);
  ::glColor3d(0,0,0);
  for(unsigned int i0=0;i0<rs.shape.size();++i0){
    unsigned int i1 = (i0+1)%rs.shape.size();
    double p0[2]; dfm2::Vec2_Mat3Vec2_AffineProjection(p0,mT1RT0.mat,rs.shape[i0].p);
    double p1[2]; dfm2::Vec2_Mat3Vec2_AffineProjection(p1,mT1RT0.mat,rs.shape[i1].p);
    ::glVertex2dv(p0);
    ::glVertex2dv(p1);
  }
  ::glEnd();
}

bool Nearest_Polygon2(
    dfm2::CVec2d& p_near,
    dfm2::CVec2d& n_near,
    const std::vector<dfm2::CVec2d>& shape,
    const dfm2::CVec2d& q)
{
  const unsigned int nej = shape.size();
  double min_dist = -1.0;
  double sum_theta = 0.0;
  for(unsigned int iej=0;iej<nej;++iej) {
    const dfm2::CVec2d p0 = shape[(iej + 0) % nej];
    const dfm2::CVec2d p1 = shape[(iej + 1) % nej];
    double theta = atan2((p0-q)^(p1-q),(p0-q)*(p1-q));
    sum_theta += theta;
    const dfm2::CVec2d pe = dfm2::GetNearest_LineSeg_Point(q, p0, p1);
    const dfm2::CVec2d ne = (p1 - p0).Rotate(-M_PI * 0.5).Normalize();
    const double dist0 = (pe-q).Length();
    if( min_dist > 0 && dist0 > min_dist ){ continue; }
    //
    min_dist = dist0;
    p_near = pe;
    if( dist0 < 1.0e-5 ){ n_near = ne; }
    else{ n_near = (pe-q).Normalize(); }
  }
//  std::cout << "sum_t" << sum_theta << std::endl;
  return fabs(sum_theta - M_PI*2)<1.0e-3;
}

class CContact {
public:
  unsigned int ir, jr;
  dfm2::CVec2d Pi, Pjn, Njn;
  double lambda;
};

void Steptime_Rgd2(
    std::vector<CRigidState2>& aRS,
    double dt,
    const dfm2::CVec2d& gravity)
{
  for(CRigidState2& rs : aRS){
    if( rs.is_fix ) {
      rs.posg_tmp = rs.posg;
      rs.theta_tmp = rs.theta;
      continue;
    }
    rs.velo += gravity * dt;
    rs.posg_tmp = rs.posg + dt * rs.velo;
    rs.theta_tmp = rs.theta + dt * rs.omega;
  }
  std::vector<CContact> aContact;
  for(unsigned int ir=0;ir<aRS.size();++ir) {
    CRigidState2& ri = aRS[ir];
    const dfm2::CMat3d Ai = Mat3_Affine(ri.posl,ri.theta_tmp,ri.posg_tmp);
    for(dfm2::CVec2d& Pi : ri.shape) {
      const dfm2::CVec2d pi = Pi.Mat3Vec2_AffineProjection(Ai.mat);
      for(unsigned int jr=0;jr<aRS.size();++jr) {
        if( ir == jr ){ continue; }
        CRigidState2& rj = aRS[jr];
        const dfm2::CMat3d Aj = Mat3_Affine(rj.posl,rj.theta_tmp,rj.posg_tmp);
        const dfm2::CVec2d Pj = Pi.Mat3Vec2_AffineProjection((Aj.Inverse() * Ai).mat);
        dfm2::CVec2d Pjn, Njn;
        const bool is_inside = Nearest_Polygon2(Pjn, Njn, rj.shape, Pj);
        if( !is_inside ){ continue; }
//          std::cout << " a " << (Pj-Pjn).Length() << std::endl;
//          std::cout << ir << " " << jr << "  ( " << Pi << "  )( " << Pjn << "  )( " << Njn << ")" << std::endl;
        const double penetration = (Pjn-Pj)*Njn;
        const dfm2::CVec2d pj = Pjn.Mat3Vec2_AffineProjection(Aj.mat);
        const dfm2::CVec2d nj = Njn.Mat3Vec2_AffineDirection(Aj.mat);
        double deno = 0.0;
        if( !ri.is_fix ){
          deno += 1/ri.mass;
          const double t0 = (pi-ri.posg) ^ nj;
          deno += t0*t0/ri.I;
        }
        if( !rj.is_fix ){
          deno += 1/rj.mass;
          const double t0 = (pj-rj.posg) ^ (-nj);
          deno += t0*t0/rj.I;
        }
        const double lambda = penetration/deno; // force
        if( !ri.is_fix ){
          ri.posg_tmp += (lambda/ri.mass)*nj;
          ri.theta_tmp += ((pi-ri.posg) ^ nj) * lambda / ri.I;
        }
        if( !rj.is_fix ){
          rj.posg_tmp += (lambda/rj.mass)*(-nj);
          rj.theta_tmp += ((pj-rj.posg) ^ (-nj)) * lambda / rj.I;
        }
        aContact.push_back({ir,jr,Pi,Pjn,Njn,lambda});
      }
    }
  }
  for(CRigidState2& rs : aRS) {
    if( rs.is_fix ) {
      rs.velo = dfm2::CVec2d(0,0);
      rs.omega = 0.0;
      continue;
    }
    rs.velo = (rs.posg_tmp - rs.posg)/dt;
    rs.omega = (rs.theta_tmp - rs.theta)/dt;
    rs.posg = rs.posg_tmp;
    rs.theta = rs.theta_tmp;
  }
  for(CContact& c : aContact ){
    CRigidState2& ri = aRS[c.ir];
    CRigidState2& rj = aRS[c.jr];
    dfm2::CVec2d vi = ri.omega*(c.Pi - ri.posl).Rotate(ri.theta + M_PI*0.5) + ri.velo;
    dfm2::CVec2d vj = rj.omega*(c.Pjn - rj.posl).Rotate(rj.theta + M_PI*0.5) + rj.velo;
    dfm2::CVec2d pi = (c.Pi - ri.posl).Rotate(ri.theta) + ri.posg;
    dfm2::CVec2d pj = (c.Pjn - rj.posl).Rotate(rj.theta) + rj.posg;
    dfm2::CVec2d tj = c.Njn.Rotate(rj.theta - M_PI*0.5);
    double vslip = (vi-vj)*tj;
    double force = c.lambda/(dt*dt);
    double deno = 0.0;
    if( !ri.is_fix ){
      deno += 1/ri.mass;
      const double t0 = (pi-ri.posg) ^ tj;
      deno += t0*t0/ri.I;
    }
    if( !rj.is_fix ){
      deno += 1/rj.mass;
      const double t0 = (pj-rj.posg) ^ tj;
      deno += t0*t0/rj.I;
    }
    double mass = 1.0/deno;
    double dvelo = force/mass*dt; // maximum velocity difference caused by friction
    if( dvelo < fabs(vslip) ){
      double sign = (vslip < 0 ) ? -1 : +1;
      vslip = dvelo * sign;
    }
    if( !ri.is_fix ) {
      ri.velo += vslip * (-tj) * (1 / ri.mass) / deno;
      const double t0 = (pi-ri.posg) ^ (-tj);
      ri.omega += vslip*(t0/ri.I)/deno;
    }
    if( !rj.is_fix ) {
      rj.velo += vslip * tj  * ( 1 / rj.mass) / deno;
      const double t0 = (pj-rj.posg) ^ tj;
      rj.omega += vslip*(t0/rj.I)/deno;
    }
  }
}

int main(int argc,char* argv[])
{
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.camera.view_height = 1.5;

  std::vector<CRigidState2> aRS(4);
  {
    CRigidState2& rs = aRS[0];
    rs.is_fix = true;
    rs.shape = {
        {-1.0, 0.0},
        {+1.0, 0.0},
        {+1.0, 0.8},
        {+0.9, 0.8},
        {+0.9, 0.2},
        {-0.9, 0.2},
        {-0.9, 0.8},
        {-1.0, 0.8},
    };
  }
  {
    CRigidState2& rs = aRS[1];
    rs.is_fix = false;
    rs.shape = { {0.0, 0.0}, {0.4, 0.0}, {0.4, 0.2}, {0.0, 0.2} };
  }
  {
    CRigidState2& rs = aRS[2];
    rs.is_fix = false;
    rs.shape = { {0.0, 0.0}, {0.4, 0.0}, {0.4, 0.2}, {0.0, 0.2} };
  }
  {
    CRigidState2& rs = aRS[3];
    rs.is_fix = false;
    rs.shape = { {0.0, 0.0}, {0.4, 0.0}, {0.4, 0.2}, {0.0, 0.2} };
  }

  for(CRigidState2& rs : aRS){
    dfm2::CgArea_Polygon(rs.posl,rs.mass, rs.shape);
    rs.I = dfm2::RotationalMomentPolar_Polygon2(rs.shape,rs.posl);
    double rho = 1.0;
    rs.mass *= rho;
    rs.I *= rho;
    rs.omega = 0.0;
    rs.theta = 0;
    rs.velo = dfm2::CVec2d(0, 0);
    rs.posg = dfm2::CVec2d(0, 0.0);
  }

  aRS[1].posg.p[1] = 0.5;
  aRS[1].theta = M_PI*0.1;

  aRS[2].posg.p[0] = 0.1;
  aRS[2].posg.p[1] = 1.0;

  aRS[3].posg.p[0] = -0.2;
  aRS[3].posg.p[1] = 1.3;

  const dfm2::CVec2d gravity(0,-10);
  double dt = 0.005;

  while(true){
    Steptime_Rgd2(aRS, dt, gravity);
    //
    viewer.DrawBegin_oldGL();
    for(const CRigidState2& rs : aRS) {
      Draw(rs);
    }
    viewer.SwapBuffers();
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ break; }
  }
  viewer.ExitIfClosed();
  return 0;
}


