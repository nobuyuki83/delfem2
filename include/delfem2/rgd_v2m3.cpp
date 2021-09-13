#include "delfem2/rgd_v2m3.h"

#ifndef M_PI
  #define M_PI 3.141592653589793
#endif  


namespace delfem2{
namespace rgd_v2m3{

CMat3d Mat3_Affine(
    const CVec2d& posl,
    const double theta,
    const CVec2d& posg)
{
  CMat3d mT0; Mat3_AffineTranslation(mT0.p_,(-posl).p);
  CMat3d mR; Mat3_AffineRotation(mR.p_,theta);
  CMat3d mT1; Mat3_AffineTranslation(mT1.p_,posg.p);
  return mT1*mR*mT0;
}

/**
 *
 * @param[out] ie_near
 * @param[out] re_near
 * @param[out] n_near normal
 * @param[in] shape of the polygon
 * @param[in] q input point
 * @return true: collision, false: no-collision
 */
bool Nearest_Polygon2(
    unsigned int& ie_near,
    double& re_near,
    CVec2d& n_near,
    const std::vector<CVec2d>& shape,
    const CVec2d& q)
{
  const unsigned int nej = static_cast<unsigned int>(shape.size());
  double min_dist = -1.0;
  double winding_nubmer = 0.0;
  for(unsigned int iej=0;iej<nej;++iej) {
    const CVec2d p0 = shape[(iej + 0) % nej];
    const CVec2d p1 = shape[(iej + 1) % nej];
    double theta = atan2((p0-q)^(p1-q),(p0-q).dot(p1-q));
    winding_nubmer += theta;
    const CVec2d pe = GetNearest_LineSeg_Point(q, p0, p1);
    double re0 = (pe-p0).norm()/(p1-p0).norm();
    const CVec2d ne = (p1 - p0).Rotate(-M_PI * 0.5).normalized();
    const double dist0 = (pe-q).norm();
    if( min_dist > 0 && dist0 > min_dist ){ continue; }
    //
    min_dist = dist0;
    ie_near = iej;
    re_near = re0;
    if( dist0 < 1.0e-5 ){ n_near = ne; } // if distance is small use edge's normal
    else{ n_near = (pe-q).normalized(); } // if distance is not so small, use the contacting direciton
  }
  return fabs(winding_nubmer - M_PI*2)<1.0e-3; // winding_number==M_PI: inside else: outside
}

/**
 *
 * @param[in,out] ri rigid body A
 * @param[in,out] rj rigid body B
 * @param[in] penetration penetration depth
 * @param[in] pi contacting point on A in global coordinate
 * @param[in] pj contacting point on B in global coordinate
 * @param[in] nj contact normal on B in global coordinate
 * @return force*dt*dt to resolve contact
 */
double ResolveContact(
    CRigidState2& rbA,
    CRigidState2& rbB,
    double penetration,
    const CVec2d& pA,
    const CVec2d& pB,
    const CVec2d& nB)
{
  double deno = 0.0;
  if( !rbA.is_fix ){
    deno += 1/rbA.mass;
    const double t0 = (pA-rbA.posg) ^ nB;
    deno += t0*t0/rbA.I;
  }
  if( !rbB.is_fix ){
    deno += 1/rbB.mass;
    const double t0 = (pB-rbB.posg) ^ (-nB);
    deno += t0*t0/rbB.I;
  }
  const double lambda = penetration/deno; // force*dt*dt
  if( !rbA.is_fix ){
    rbA.posg_tmp += (lambda/rbA.mass)*nB;
    rbA.theta_tmp += ((pA-rbA.posg) ^ nB) * lambda / rbA.I;
  }
  if( !rbB.is_fix ){
    rbB.posg_tmp += (lambda/rbB.mass)*(-nB);
    rbB.theta_tmp += ((pB-rbB.posg) ^ (-nB)) * lambda / rbB.I;
  }
  return lambda;
}

/**
 *
 * @param[in,out] rbA rigid body A
 * @param[in,out] rbB rigid body B
 * @param[in] pA contact point on A in global coordinate
 * @param[in] pB contact point on B in global coordinate
 * @param[in] dir_tangent  tangent direction in global coordinate
 * @param[in] vslip // slipping velocity
 * @param impulse force*dt: perpendiculer force times time step.
 */
void AddFriction(
    CRigidState2& rbA,
    CRigidState2& rbB,
    const CVec2d& pA,
    const CVec2d& pB,
    const CVec2d& dir_tangent,
    double vslip,
    double impulse)
{
//  double contact_force = lambda/(dt*dt); // force=implulse/(dt*dt)
  double deno = 0.0;
  if( !rbA.is_fix ){
    deno += 1/rbA.mass;
    const double t0 = (pA-rbA.posg) ^ dir_tangent;
    deno += t0*t0/rbA.I;
  }
  if( !rbB.is_fix ){
    deno += 1/rbB.mass;
    const double t0 = (pB-rbB.posg) ^ dir_tangent;
    deno += t0*t0/rbB.I;
  }
  const double mass = 1.0/deno;
  const double dvelo = impulse/mass; // maximum velocity difference caused by friction
  if( dvelo < fabs(vslip) ){
    double sign = (vslip < 0 ) ? -1 : +1;
    vslip = dvelo * sign;
  }
  if( !rbA.is_fix ) {
    rbA.velo += vslip * (-dir_tangent) * (1 / rbA.mass) / deno;
    const double t0 = (pA-rbA.posg) ^ (-dir_tangent);
    rbA.omega += vslip*(t0/rbA.I)/deno;
  }
  if( !rbB.is_fix ) {
    rbB.velo += vslip * dir_tangent  * ( 1 / rbB.mass) / deno;
    const double t0 = (pB-rbB.posg) ^ dir_tangent;
    rbB.omega += vslip*(t0/rbB.I)/deno;
  }
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
  std::vector<CContactInfo2> aContact;
  for(unsigned int irbA=0;irbA<aRS.size();++irbA) { // loop over rigid bodies
    for(unsigned int ip=0;ip<aRS[irbA].shape.size();++ip) {
      for(unsigned int irbB=0;irbB<aRS.size();++irbB) { // loop over rigid bodies
        if( irbA == irbB ){ continue; }
        const CRigidState2& rbA = aRS[irbA];
        const CMat3d& matAffineA = rgd_v2m3::Mat3_Affine(rbA.posl,rbA.theta_tmp,rbA.posg_tmp);
        const CVec2d& pA = rbA.shape[ip].Mat3Vec2_AffineProjection(matAffineA.p_);
        const CRigidState2& rbB = aRS[irbB];
        const CMat3d& matAffineB = rgd_v2m3::Mat3_Affine(rbB.posl,rbB.theta_tmp,rbB.posg_tmp); // j-Rigidbody's affine matrix
        const CVec2d& pAonB = rbA.shape[ip].Mat3Vec2_AffineProjection((matAffineB.Inverse() * matAffineA).p_); // Pi in j-Rigidbody's coordinate
        unsigned int ieB;
        double reB;
        CVec2d NrmB;
        const bool is_inside = rgd_v2m3::Nearest_Polygon2(ieB, reB, NrmB, rbB.shape, pAonB);
        if( !is_inside ){ continue; } // not penetrating
        //
        CVec2d PB = (1-reB)*rbB.shape[(ieB+0)%rbB.shape.size()] + reB*rbB.shape[(ieB+1)%rbB.shape.size()];
        const double penetration = (PB-pAonB).dot(NrmB);
        const CVec2d pB = PB.Mat3Vec2_AffineProjection(matAffineB.p_);
        const CVec2d nrmB = NrmB.Mat3Vec2_AffineDirection(matAffineB.p_);
        const double lambda = rgd_v2m3::ResolveContact(aRS[irbA],aRS[irbB],
                                                       penetration, pA, pB, nrmB);
        aContact.push_back({irbA,irbB,ip,ieB,reB,NrmB,lambda});
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
  for(CContactInfo2& c : aContact ){
    CRigidState2& rbA = aRS[c.irbA];
    CRigidState2& rbB = aRS[c.irbB];
    assert( rbA.shape.size() == rbA.shape_velo.size() );
    assert( rbB.shape.size() == rbB.shape_velo.size() );
    const CVec2d& PA = aRS[c.irbA].shape[c.ipA];
    const CVec2d& VA = aRS[c.irbA].shape_velo[c.ipA];
    const unsigned int ieB = c.ieB;
    const double reB = c.reB;
    const unsigned int npB = static_cast<unsigned int>(rbB.shape.size());
    const CVec2d PB = (1-reB)*rbB.shape[(ieB+0)%npB] + reB*rbB.shape[(ieB+1)%npB];
    const CVec2d VB = (1-reB)*rbB.shape_velo[(ieB+0)%npB] + reB*rbB.shape_velo[(ieB+1)%npB];
    const CVec2d vA = rbA.omega*(PA - rbA.posl).Rotate(rbA.theta + M_PI*0.5) + rbA.velo + VA.Rotate(rbA.theta);
    const CVec2d vB = rbB.omega*(PB - rbB.posl).Rotate(rbB.theta + M_PI*0.5) + rbB.velo + VB.Rotate(rbB.theta);
    const CVec2d pA = (PA - rbA.posl).Rotate(rbA.theta) + rbA.posg;
    const CVec2d pB = (PB - rbB.posl).Rotate(rbB.theta) + rbB.posg;
    const CVec2d tangent_dir = c.Njn.Rotate(rbB.theta - M_PI*0.5); // tangent direction
    double velo_slip = (vA-vB).dot(tangent_dir); // slipping velocity
    rgd_v2m3::AddFriction(rbA,rbB,
                          pA,pB,tangent_dir,velo_slip,c.lambda/dt);
  }
}
