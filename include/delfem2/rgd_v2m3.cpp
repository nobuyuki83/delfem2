#include "delfem2/rgd_v2m3.h"

namespace delfem2{
namespace rgd_v2m3{

CMat3d Mat3_Affine(
    const CVec2d& posl,
    const double theta,
    const CVec2d& posg)
{
  CMat3d mT0; Mat3_AffineTranslation(mT0.mat,(-posl).p);
  CMat3d mR; Mat3_AffineRotation(mR.mat,theta);
  CMat3d mT1; Mat3_AffineTranslation(mT1.mat,posg.p);
  return mT1*mR*mT0;
}

bool Nearest_Polygon2(
    CVec2d& p_near,
    CVec2d& n_near,
    const std::vector<CVec2d>& shape,
    const CVec2d& q)
{
  const unsigned int nej = shape.size();
  double min_dist = -1.0;
  double sum_theta = 0.0;
  for(unsigned int iej=0;iej<nej;++iej) {
    const CVec2d p0 = shape[(iej + 0) % nej];
    const CVec2d p1 = shape[(iej + 1) % nej];
    double theta = atan2((p0-q)^(p1-q),(p0-q)*(p1-q));
    sum_theta += theta;
    const CVec2d pe = GetNearest_LineSeg_Point(q, p0, p1);
    const CVec2d ne = (p1 - p0).Rotate(-M_PI * 0.5).Normalize();
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

}
}

void delfem2::Steptime_Rgd2(
    std::vector<delfem2::CRigidState2>& aRS,
    double dt,
    const delfem2::CVec2d& gravity)
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
  // find & resolve contact
  std::vector<CContact> aContact;
  for(unsigned int ir=0;ir<aRS.size();++ir) {
    CRigidState2& ri = aRS[ir];
    const CMat3d Ai = rgd_v2m3::Mat3_Affine(ri.posl,ri.theta_tmp,ri.posg_tmp);
    for(CVec2d& Pi : ri.shape) {
      const CVec2d pi = Pi.Mat3Vec2_AffineProjection(Ai.mat);
      for(unsigned int jr=0;jr<aRS.size();++jr) {
        if( ir == jr ){ continue; }
        CRigidState2& rj = aRS[jr];
        const CMat3d Aj = rgd_v2m3::Mat3_Affine(rj.posl,rj.theta_tmp,rj.posg_tmp);
        const CVec2d Pj = Pi.Mat3Vec2_AffineProjection((Aj.Inverse() * Ai).mat);
        CVec2d Pjn, Njn;
        const bool is_inside = rgd_v2m3::Nearest_Polygon2(Pjn, Njn, rj.shape, Pj);
        if( !is_inside ){ continue; }
//          std::cout << " a " << (Pj-Pjn).Length() << std::endl;
//          std::cout << ir << " " << jr << "  ( " << Pi << "  )( " << Pjn << "  )( " << Njn << ")" << std::endl;
        const double penetration = (Pjn-Pj)*Njn;
        const CVec2d pj = Pjn.Mat3Vec2_AffineProjection(Aj.mat);
        const CVec2d nj = Njn.Mat3Vec2_AffineDirection(Aj.mat);
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
  // finalize position, set velocity
  for(CRigidState2& rs : aRS) {
    if( rs.is_fix ) {
      rs.velo = CVec2d(0,0);
      rs.omega = 0.0;
      continue;
    }
    rs.velo = (rs.posg_tmp - rs.posg)/dt;
    rs.omega = (rs.theta_tmp - rs.theta)/dt;
    rs.posg = rs.posg_tmp;
    rs.theta = rs.theta_tmp;
  }
  // add frictional velocity change
  for(CContact& c : aContact ){
    CRigidState2& ri = aRS[c.ir];
    CRigidState2& rj = aRS[c.jr];
    const CVec2d vi = ri.omega*(c.Pi - ri.posl).Rotate(ri.theta + M_PI*0.5) + ri.velo;
    const CVec2d vj = rj.omega*(c.Pjn - rj.posl).Rotate(rj.theta + M_PI*0.5) + rj.velo;
    const CVec2d pi = (c.Pi - ri.posl).Rotate(ri.theta) + ri.posg;
    const CVec2d pj = (c.Pjn - rj.posl).Rotate(rj.theta) + rj.posg;
    const CVec2d tj = c.Njn.Rotate(rj.theta - M_PI*0.5);
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