/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/eigen/eigen_rigidbody.h"

#if defined(__APPLE__) && defined(__MACH__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/glu.h>
#elif defined(_WIN32)
#include <windows.h>
#include <GL/glu.h>
#else
#include <GL/glu.h>
#endif

#if defined(__APPLE__) && defined(__MACH__)
#include "Eigen/Dense"
#else
#include <Eigen/Dense>
#endif

#include "delfem2/mat3_funcs.h"

void myUpdateMinMaxXYZ(
    double bb[6],
    const delfem2::CVec3d p)
{
  const double x = p.x;
  const double y = p.y;
  const double z = p.z;
  if( bb[0] > bb[1] ){
    bb[0] = bb[1] = x;
    bb[2] = bb[3] = y;
    bb[4] = bb[5] = z;
    return;
  }
  bb[0] = (bb[0] < x) ? bb[0] : x;
  bb[1] = (bb[1] > x) ? bb[1] : x;
  bb[2] = (bb[2] < y) ? bb[2] : y;
  bb[3] = (bb[3] > y) ? bb[3] : y;
  bb[4] = (bb[4] < z) ? bb[4] : z;
  bb[5] = (bb[5] > z) ? bb[5] : z;
}

// -----------------------------------------------

void EdEd_Potential(
    double& energy,
    delfem2::CVec3d& dEdu,
    delfem2::CVec3d& dEdw,
    const delfem2::CVec3d& cg,
    const delfem2::CVec3d& u,
    const double mass,
    const delfem2::CVec3d& g       // floor normal
)
{
  energy = mass*u.dot(g);
  dEdu = mass*g;
  dEdw = delfem2::CVec3d(0,0,0);
}

void EdEddE_Exforce(
    double& energy,
    delfem2::CVec3d& dEdu,
    delfem2::CVec3d& dEdw,
    delfem2::CMat3d&  ddEddu,
    delfem2::CMat3d&  ddEddw,
    delfem2::CMat3d&  ddEdudw,
    const delfem2::CVec3d& cg, // the center of gravity
    const delfem2::CVec3d& pex, // external force position
    const delfem2::CVec3d& fex, // external force
    const delfem2::CVec3d& u, // displacement
    const delfem2::CMat3d& R // rigid rotation
)
{
  delfem2::CVec3d Rv = R * (pex-cg);
  delfem2::CVec3d qex = Rv + cg + u; // current external force position
  energy = +qex.dot(fex);
  dEdu = +fex;
  dEdw = +Rv.cross(fex);
  ddEddu  = delfem2::CMat3d(0.0);
  ddEddw  = +delfem2::CMat3d(Mat3_Spin(fex))*delfem2::CMat3d(Mat3_Spin(Rv));
  ddEdudw = delfem2::CMat3d(0.0);
}

void EdEddE_Contact(
    double& energy,
    delfem2::CVec3d& dEdu,
    delfem2::CVec3d& dEdw,
    delfem2::CMat3d&  ddEddu,
    delfem2::CMat3d&  ddEddw,
    delfem2::CMat3d&  ddEdudw,
    const delfem2::CVec3d& cg, // the center of gravity
    const delfem2::CVec3d& cp, // contact position
    const delfem2::CVec3d& u, // displacement
    const delfem2::CMat3d& R, // rigid rotation
    const double cont_stiff,
    const delfem2::CVec3d& n       // floor normal
)
{
  delfem2::CVec3d Rv = R * (cp-cg);
  delfem2::CVec3d cq = Rv + cg + u;
  energy = 0.5*cq.dot(n)*cq.dot(n)*cont_stiff;
  dEdu = (cq.dot(n)*cont_stiff)*n;
  dEdw = (cq.dot(n)*cont_stiff)*Rv.cross(n);
  const delfem2::CMat3d m0 = delfem2::Mat3_OuterProduct(n,n);
  const delfem2::CMat3d m1 = delfem2::Mat3_OuterProduct(Rv.cross(n),Rv.cross(n));
  const delfem2::CMat3d m2 = delfem2::CMat3d(Mat3_Spin(n))*delfem2::CMat3d(Mat3_Spin(Rv));
  const delfem2::CMat3d m3 = delfem2::Mat3_OuterProduct(Rv.cross(n),n);
  ddEddu  = cont_stiff*m0;
  ddEddw  = cont_stiff*m1 + (cq.dot(n)*cont_stiff)*m2;
  ddEdudw = cont_stiff*m3;
}

void EdEddE_ContactFriction(
    double& energy,
    delfem2::CVec3d& dEdu,
    delfem2::CVec3d& dEdw,
    delfem2::CMat3d&  ddEddu,
    delfem2::CMat3d&  ddEddw,
    delfem2::CMat3d&  ddEdudw,
    const delfem2::CVec3d& cg, // the center of gravity
    const delfem2::CVec3d& cp, // contact position
    const delfem2::CVec3d& u, // displacement
    const delfem2::CMat3d& R, // rigid rotation
    const double cont_stiff)
{
  delfem2::CVec3d Rv = R * (cp-cg);
  delfem2::CVec3d cq = Rv + cg + u;
  energy = 0.5*(cq-cp).dot(cq-cp)*cont_stiff;
  dEdu = cont_stiff*(cq-cp);
  dEdw = cont_stiff*(Rv.cross(cq-cp));
  ddEddu  = cont_stiff*delfem2::CMat3d::Identity();
  //    ddEddw  = -cont_stiff*delfem2::CMat3d::OuterProduct(Rv,Rv) + cont_stiff*delfem2::CMat3d::Spin(Rv)*delfem2::CMat3d::Spin(cq-cp);
  delfem2::CMat3d m0 = delfem2::CMat3d(Mat3_Spin(Rv))*delfem2::CMat3d(Mat3_Spin(Rv));
  delfem2::CMat3d m1 = delfem2::CMat3d(Mat3_Spin(cq-cp))*delfem2::CMat3d(Mat3_Spin(Rv));
  delfem2::CMat3d m2 = Mat3_Spin(Rv);
  ddEddw  = -cont_stiff*m0 + cont_stiff*m1;
  ddEdudw = cont_stiff*m2;
}


void EdEddE_Joint(
    double& energy,
    delfem2::CVec3d& dEdu0,
    delfem2::CVec3d& dEdw0,
    delfem2::CVec3d& dEdu1,
    delfem2::CVec3d& dEdw1,
    delfem2::CMat3d& ddEdu0du0,
    delfem2::CMat3d& ddEdu0dw0,
    delfem2::CMat3d& ddEdu0du1,
    delfem2::CMat3d& ddEdu0dw1,
    delfem2::CMat3d& ddEdw0dw0,
    delfem2::CMat3d& ddEdw0du1,
    delfem2::CMat3d& ddEdw0dw1,
    delfem2::CMat3d& ddEdu1du1,
    delfem2::CMat3d& ddEdu1dw1,
    delfem2::CMat3d& ddEdw1dw1,
    const double trans_stiff,
    const double rot_stiff,
    const delfem2::CVec3d& pj,
    const delfem2::CVec3d& cg0,
    const delfem2::CVec3d& u0,
    const delfem2::CMat3d& R0,
    const delfem2::CVec3d& cg1,
    const delfem2::CVec3d& u1,
    const delfem2::CMat3d& R1)
{
  delfem2::CVec3d Rv0 = R0 * (pj-cg0);
  delfem2::CVec3d qj0 = Rv0 + cg0 + u0; // after deformation joint pos relative to rigid body 0
  
  delfem2::CVec3d Rv1 = R1 * (pj-cg1);
  delfem2::CVec3d qj1 = Rv1 + cg1 + u1; // after deformation joint pos relative to rigid body 1
  
  energy = 0.5*trans_stiff*(qj0 - qj1).squaredNorm();
  
  dEdu0 = trans_stiff*(qj0-qj1);
  dEdw0 = trans_stiff*(Rv0.cross(qj0-qj1));
  
  dEdu1 = trans_stiff*(qj1-qj0);
  dEdw1 = trans_stiff*(Rv1.cross(qj1-qj0));
  
  ddEdu0du0 =  trans_stiff*delfem2::CMat3d::Identity();
  ddEdu0du1 = -trans_stiff*delfem2::CMat3d::Identity();
  ddEdu0dw0 = +trans_stiff*delfem2::CMat3d(Mat3_Spin(Rv0));
  ddEdu0dw1 = -trans_stiff*delfem2::CMat3d(Mat3_Spin(Rv1));

  delfem2::CMat3d m0 = delfem2::CMat3d(Mat3_Spin(Rv0))*delfem2::CMat3d(Mat3_Spin(Rv0));
  delfem2::CMat3d m1 = delfem2::CMat3d(Mat3_Spin(qj0-qj1))*delfem2::CMat3d(Mat3_Spin(Rv0));
  delfem2::CMat3d m2 = delfem2::CMat3d(Mat3_Spin(Rv1))*delfem2::CMat3d(Mat3_Spin(Rv0));
  ddEdw0dw0 = -trans_stiff*m0 + trans_stiff*m1;
  ddEdw0du1 = +trans_stiff*delfem2::CMat3d(Mat3_Spin(Rv0));
  ddEdw0dw1 = +trans_stiff*m2;
  
  ddEdu1du1 =  trans_stiff*delfem2::CMat3d::Identity();
  ddEdu1dw1 = +trans_stiff*delfem2::CMat3d(Mat3_Spin(Rv1));

  delfem2::CMat3d m3 = delfem2::CMat3d(Mat3_Spin(Rv1))*delfem2::CMat3d(Mat3_Spin(Rv1));
  delfem2::CMat3d m4 = delfem2::CMat3d(Mat3_Spin(qj1-qj0))*delfem2::CMat3d(Mat3_Spin(Rv1));
  ddEdw1dw1 = -trans_stiff*m3 + trans_stiff*m4;
  
  delfem2::CVec3d av(0,0,0);
  delfem2::CMat3d davdw0, davdw1;
  {
    for(unsigned int i=0;i<3;i++){
      delfem2::CVec3d r0( R0.p_[0*3+i], R0.p_[1*3+i], R0.p_[2*3+i] );
      delfem2::CVec3d r1( R1.p_[0*3+i], R1.p_[1*3+i], R1.p_[2*3+i] );
      av += r0.cross(r1);
      davdw0 += delfem2::CMat3d(Mat3_Spin(r1))*delfem2::CMat3d(Mat3_Spin(r0));
      davdw1 -= delfem2::CMat3d(Mat3_Spin(r0))*delfem2::CMat3d(Mat3_Spin(r1));
    }
    av *= 0.5;
    davdw0 *= 0.5;
    davdw1 *= 0.5;
  }
  energy += 0.5*rot_stiff*av.squaredNorm();
  
  // delfem2::CMat3d m0,m1,m2;
  for(unsigned int i=0;i<3;i++){
    delfem2::CVec3d r0( R0.p_[0*3+i], R0.p_[1*3+i], R0.p_[2*3+i] );
    delfem2::CVec3d r1( R1.p_[0*3+i], R1.p_[1*3+i], R1.p_[2*3+i] );
    dEdw0 += 0.5*rot_stiff*r0.cross(r1.cross(av));
    dEdw1 -= 0.5*rot_stiff*r1.cross(r0.cross(av));
    delfem2::CMat3d l0 = delfem2::CMat3d(Mat3_Spin(r1.cross(av)))*delfem2::CMat3d(Mat3_Spin(r0));
    delfem2::CMat3d l1 = delfem2::CMat3d(Mat3_Spin(r0))*delfem2::CMat3d(Mat3_Spin(r1));
    delfem2::CMat3d l2 = delfem2::CMat3d(Mat3_Spin(r0.cross(av)))*delfem2::CMat3d(Mat3_Spin(r1));
    delfem2::CMat3d l3 = delfem2::CMat3d(Mat3_Spin(r1))*delfem2::CMat3d(Mat3_Spin(r0));
    delfem2::CMat3d l4 = delfem2::CMat3d(Mat3_Spin(r1))*delfem2::CMat3d(Mat3_Spin(av))*delfem2::CMat3d(Mat3_Spin(r0));
    delfem2::CMat3d l5 = delfem2::CMat3d(Mat3_Spin(r1))*delfem2::CMat3d(Mat3_Spin(r0));
    ddEdw0dw0 += 0.5*rot_stiff*( l0 + l1*davdw0 );
    ddEdw1dw1 -= 0.5*rot_stiff*( l2 + l3*davdw1 );
    ddEdw0dw1 -= 0.5*rot_stiff*( l4 + l5*davdw0 );
  }
}


delfem2::CVec3d rand_vec(double s)
{
  delfem2::CVec3d v;
  v.x = s*rand()/(RAND_MAX+1.0);
  v.y = s*rand()/(RAND_MAX+1.0);
  v.z = s*rand()/(RAND_MAX+1.0);
  return v;
}

delfem2::CMat3d rand_rot()
{
  double s = 3.5;
  delfem2::CVec3d v;
  v.x = s*rand()/(RAND_MAX+1.0);
  v.y = s*rand()/(RAND_MAX+1.0);
  v.z = s*rand()/(RAND_MAX+1.0);
  delfem2::CMat3d R = delfem2::Mat3_RotMatFromAxisAngleVec(v);
  return R;
}

void CheckDiff_Contact()
{
  double energy = 0.0;
  delfem2::CVec3d dEdu,dEdw;
  delfem2::CMat3d ddEddu, ddEddw, ddEdudw;
  double epsilon = 1.0e-5;
  const double cont_stiff = 1.0e+3;
  delfem2::CVec3d cg = rand_vec(1.0);
  delfem2::CVec3d cp = rand_vec(1.0);
  delfem2::CVec3d u  = rand_vec(1.0);
  delfem2::CMat3d  R  = rand_rot();
  delfem2::CVec3d n(0,0,1);
  EdEddE_Contact(energy,
                 dEdu,dEdw,  ddEddu,ddEddw,ddEdudw,
                 cg,cp,u,R,
                 cont_stiff,n);
  for(unsigned int idim=0;idim<3;idim++){
    delfem2::CVec3d dEdu_, dEdw_;
    delfem2::CMat3d ddEddu_, ddEddw_, ddEdudw_;
    double energy_ = 0.0;
    delfem2::CVec3d u_ = u;
    u_[idim] += epsilon;
    EdEddE_Contact(energy_,  dEdu_,dEdw_,   ddEddu_,ddEddw_,ddEdudw_,
                   cg,cp,u_,R,
                   cont_stiff,n);
    std::cout << "dEdu " << idim << " --> "  << (energy_-energy)/epsilon << " " << dEdu[idim] << std::endl;
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEddu  " << idim << " " << jdim << " --> " << (dEdu_[jdim]-dEdu[jdim])/epsilon << " " << ddEddu.p_[jdim*3+idim] << std::endl;
    }
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEdudw " << idim << " " << jdim << " --> " << (dEdw_[jdim]-dEdw[jdim])/epsilon << " " << ddEdudw.p_[jdim*3+idim] << std::endl;
    }
  }
  for(unsigned int idim=0;idim<3;idim++){
    delfem2::CVec3d dEdu_, dEdw_;
    delfem2::CMat3d ddEddu_, ddEddw_, ddEdudw_;
    double energy_ = 0.0;
    delfem2::CVec3d w(0,0,0);
    w[idim] = epsilon;
    delfem2::CMat3d dR = delfem2::Mat3_RotMatFromAxisAngleVec(w);
    delfem2::CMat3d R_ = dR*R;
    EdEddE_Contact(energy_,  dEdu_,dEdw_,   ddEddu_,ddEddw_,ddEdudw_,
                   cg,cp,u,R_,
                   cont_stiff,n);
    std::cout << "dEdw " << idim << " --> " << (energy_-energy)/epsilon << " " << dEdw[idim] << std::endl;
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEdudw " << idim << " " << jdim << " --> " << (dEdu_[jdim]-dEdu[jdim])/epsilon << " " << ddEdudw.p_[idim*3+jdim] << std::endl;
    }
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEddw  " << idim << " " << jdim << " --> " << (dEdw_[jdim]-dEdw[jdim])/epsilon << " " << ddEddw.p_[jdim*3+idim] << std::endl;
    }
  }
}


void CheckDiff_ContactFriction()
{
  double energy = 0.0;
  delfem2::CVec3d dEdu,dEdw;
  delfem2::CMat3d ddEddu, ddEddw, ddEdudw;
  double epsilon = 1.0e-5;
  const double cont_stiff = 1.0e+3;
  delfem2::CVec3d cg = rand_vec(1.0);
  delfem2::CVec3d cp = rand_vec(1.0);
  delfem2::CVec3d u  = rand_vec(1.0);
  delfem2::CMat3d  R  = rand_rot();
  EdEddE_ContactFriction(energy,
                         dEdu,dEdw,  ddEddu,ddEddw,ddEdudw,
                         cg,cp,u,R,
                         cont_stiff);
  for(unsigned int idim=0;idim<3;idim++){
    delfem2::CVec3d dEdu_, dEdw_;
    delfem2::CMat3d ddEddu_, ddEddw_, ddEdudw_;
    double energy_ = 0.0;
    delfem2::CVec3d u_ = u;
    u_[idim] += epsilon;
    EdEddE_ContactFriction(energy_,  dEdu_,dEdw_,   ddEddu_,ddEddw_,ddEdudw_,
                           cg,cp,u_,R,
                           cont_stiff);
    std::cout << "dEdu " << idim << " --> "  << (energy_-energy)/epsilon << " " << dEdu[idim] << std::endl;
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEddu  " << idim << " " << jdim << " --> " << (dEdu_[jdim]-dEdu[jdim])/epsilon << " " << ddEddu.p_[jdim*3+idim] << std::endl;
    }
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEdudw " << idim << " " << jdim << " --> " << (dEdw_[jdim]-dEdw[jdim])/epsilon << " " << ddEdudw.p_[jdim*3+idim] << std::endl;
    }
  }
  for(unsigned int idim=0;idim<3;idim++){
    delfem2::CVec3d dEdu_, dEdw_;
    delfem2::CMat3d ddEddu_, ddEddw_, ddEdudw_;
    double energy_ = 0.0;
    delfem2::CVec3d w(0,0,0);
    w[idim] = epsilon;
    delfem2::CMat3d dR = delfem2::Mat3_RotMatFromAxisAngleVec(w);
    delfem2::CMat3d R_ = dR*R;
    EdEddE_ContactFriction(energy_,  dEdu_,dEdw_,   ddEddu_,ddEddw_,ddEdudw_,
                           cg,cp,u,R_,
                           cont_stiff);
    std::cout << "dEdw " << idim << " --> " << (energy_-energy)/epsilon << " " << dEdw[idim] << std::endl;
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEdudw " << idim << " " << jdim << " --> " << (dEdu_[jdim]-dEdu[jdim])/epsilon << " " << ddEdudw.p_[idim*3+jdim] << std::endl;
    }
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEddw  " << idim << " " << jdim << " --> " << (dEdw_[jdim]-dEdw[jdim])/epsilon << " " << ddEddw.p_[jdim*3+idim] << std::endl;
    }
    
  }
}



void CheckDiff_Exforce()
{
  double energy = 0.0;
  delfem2::CVec3d dEdu,dEdw;
  delfem2::CMat3d ddEddu, ddEddw, ddEdudw;
  double epsilon = 1.0e-5;
  delfem2::CVec3d cg = rand_vec(1.0);
  delfem2::CVec3d u  = rand_vec(1.0);
  delfem2::CVec3d fex = rand_vec(1.0);
  delfem2::CVec3d pex = rand_vec(1.0);
  delfem2::CMat3d  R  = rand_rot();
  delfem2::CVec3d n(0,0,1);
  EdEddE_Exforce(energy,
                 dEdu,dEdw,  ddEddu,ddEddw,ddEdudw,
                 cg,pex,fex,u,R);
  for(unsigned int idim=0;idim<3;idim++){
    delfem2::CVec3d dEdu_, dEdw_;
    delfem2::CMat3d ddEddu_, ddEddw_, ddEdudw_;
    double energy_ = 0.0;
    delfem2::CVec3d u_ = u;
    u_[idim] += epsilon;
    EdEddE_Exforce(energy_,  dEdu_,dEdw_,   ddEddu_,ddEddw_,ddEdudw_,
                   cg,pex,fex,u_,R);
    std::cout << "dEdu " << idim << " --> "  << (energy_-energy)/epsilon << " " << dEdu[idim] << std::endl;
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEddu  " << idim << " " << jdim << " --> " << (dEdu_[jdim]-dEdu[jdim])/epsilon << " " << ddEddu.p_[jdim*3+idim] << std::endl;
    }
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEdudw " << idim << " " << jdim << " --> " << (dEdw_[jdim]-dEdw[jdim])/epsilon << " " << ddEdudw.p_[jdim*3+idim] << std::endl;
    }
  }
  for(unsigned int idim=0;idim<3;idim++){
    delfem2::CVec3d dEdu_, dEdw_;
    delfem2::CMat3d ddEddu_, ddEddw_, ddEdudw_;
    double energy_ = 0.0;
    delfem2::CVec3d w(0,0,0);
    w[idim] = epsilon;
    delfem2::CMat3d dR = delfem2::Mat3_RotMatFromAxisAngleVec(w);
    delfem2::CMat3d R_ = dR*R;
    EdEddE_Exforce(energy_,  dEdu_,dEdw_,   ddEddu_,ddEddw_,ddEdudw_,
                   cg,pex,fex,u,R_);
    std::cout << "dEdw " << idim << " --> " << (energy_-energy)/epsilon << " " << dEdw[idim] << std::endl;
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEdudw " << idim << " " << jdim << " --> " << (dEdu_[jdim]-dEdu[jdim])/epsilon << " " << ddEdudw.p_[idim*3+jdim] << std::endl;
    }
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEddw  " << idim << " " << jdim << " --> " << (dEdw_[jdim]-dEdw[jdim])/epsilon << " " << ddEddw.p_[jdim*3+idim] << std::endl;
    }
  }
}




void CheckDiff_Joint()
{
  double energy = 0.0;
  delfem2::CVec3d dEdu0(0,0,0);
  delfem2::CVec3d dEdw0(0,0,0);
  delfem2::CVec3d dEdu1(0,0,0);
  delfem2::CVec3d dEdw1(0,0,0);
  delfem2::CMat3d ddEdu0du0;
  delfem2::CMat3d ddEdu0dw0;
  delfem2::CMat3d ddEdu0du1;
  delfem2::CMat3d ddEdu0dw1;
  delfem2::CMat3d ddEdw0dw0;
  delfem2::CMat3d ddEdw0du1;
  delfem2::CMat3d ddEdw0dw1;
  delfem2::CMat3d ddEdu1du1;
  delfem2::CMat3d ddEdu1dw1;
  delfem2::CMat3d ddEdw1dw1;
  double epsilon = 1.0e-5;
  const double trans_stiff = 1.0e+3;
  const double rot_stiff = 1.0e+3;
  delfem2::CVec3d pj  = rand_vec(1.0);
  ////
  delfem2::CVec3d cg0 = rand_vec(1.0);
  delfem2::CVec3d cp0 = rand_vec(1.0);
  delfem2::CVec3d u0  = rand_vec(1.0);
  delfem2::CMat3d  R0  = rand_rot();
  ////
  delfem2::CVec3d cg1 = rand_vec(1.0);
  delfem2::CVec3d cp1 = rand_vec(1.0);
  delfem2::CVec3d u1  = rand_vec(1.0);
  delfem2::CMat3d  R1  = rand_rot();
  ////
  EdEddE_Joint(energy,
               dEdu0, dEdw0, dEdu1,dEdw1,
               ddEdu0du0, ddEdu0dw0, ddEdu0du1, ddEdu0dw1, ddEdw0dw0, ddEdw0du1, ddEdw0dw1, ddEdu1du1, ddEdu1dw1, ddEdw1dw1,
               trans_stiff, rot_stiff, pj,
               cg0, u0, R0,
               cg1, u1, R1);
  for(unsigned int kdim=0;kdim<3;kdim++){
    delfem2::CVec3d dEdu0_, dEdw0_, dEdu1_, dEdw1_;
    delfem2::CMat3d ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_;
    double energy_ = 0.0;
    delfem2::CVec3d u0_ = u0;
    u0_[kdim] += epsilon;
    EdEddE_Joint(energy_, dEdu0_,dEdw0_, dEdu1_,dEdw1_,
                 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_,
                 trans_stiff, rot_stiff, pj,
                 cg0, u0_, R0,
                 cg1, u1,  R1);
    std::cout << "dEdu0: " << kdim << " -->  " << (energy_-energy)/epsilon << " " << dEdu0[kdim] << std::endl;
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdu0du0 " << kdim << " " << idim << " --> " << (dEdu0_[idim]-dEdu0[idim])/epsilon << " " << ddEdu0du0.p_[idim*3+kdim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdw0du0 " << kdim << " " << idim << " --> " << (dEdw0_[idim]-dEdw0[idim])/epsilon << " " << ddEdu0dw0.p_[idim*3+kdim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << "*ddEdu1du0 " << kdim << " " << idim << " --> " << (dEdu1_[idim]-dEdu1[idim])/epsilon << " " << ddEdu0du1.p_[idim*3+kdim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << "*ddEdw1du0 " << kdim << " " << idim << " --> " << (dEdw1_[idim]-dEdw1[idim])/epsilon << " " << ddEdu0dw1.p_[idim*3+kdim] << std::endl;
    }
  }
  for(unsigned int jdim=0;jdim<3;jdim++){
    delfem2::CVec3d dEdu0_, dEdw0_, dEdu1_, dEdw1_;
    delfem2::CMat3d ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_;
    double energy_ = 0.0;
    delfem2::CVec3d w(0,0,0);
    w[jdim] = epsilon;
    delfem2::CMat3d dR = delfem2::Mat3_RotMatFromAxisAngleVec(w);
    delfem2::CMat3d R0_ = dR*R0;
    EdEddE_Joint(energy_, dEdu0_,dEdw0_, dEdu1_,dEdw1_,
                 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_,
                 trans_stiff, rot_stiff, pj,
                 cg0, u0, R0_,
                 cg1, u1, R1);
    std::cout << "dEdw0: " << jdim << " -->  " << (energy_-energy)/epsilon << " " << dEdw0[jdim] << std::endl;
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdu0dw0 " << jdim << " " << idim << " --> " << (dEdu0_[idim]-dEdu0[idim])/epsilon << " " << ddEdu0dw0.p_[jdim*3+idim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdw0dw0 " << jdim << " " << idim << " --> " << (dEdw0_[idim]-dEdw0[idim])/epsilon << " " << ddEdw0dw0.p_[idim*3+jdim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdu1dw0 " << jdim << " " << idim << " --> " << (dEdu1_[idim]-dEdu1[idim])/epsilon << " " << ddEdw0du1.p_[idim*3+jdim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdw1dw0 " << jdim << " " << idim << " --> " << (dEdw1_[idim]-dEdw1[idim])/epsilon << " " << ddEdw0dw1.p_[idim*3+jdim] << std::endl;
    }
  }
  
  for(unsigned int jdim=0;jdim<3;jdim++){
    delfem2::CVec3d dEdu0_, dEdw0_, dEdu1_, dEdw1_;
    delfem2::CMat3d ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_;
    double energy_ = 0.0;
    delfem2::CVec3d u1_ = u1;
    u1_[jdim] += epsilon;
    EdEddE_Joint(energy_, dEdu0_,dEdw0_, dEdu1_,dEdw1_,
                 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_,
                 trans_stiff, rot_stiff, pj,
                 cg0, u0, R0,
                 cg1, u1_,R1);
    std::cout << "dEdu1: " << jdim << " -->  " << (energy_-energy)/epsilon << " " << dEdu1[jdim] << std::endl;
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdu0du1 " << jdim << " " << idim << " --> " << (dEdu0_[idim]-dEdu0[idim])/epsilon << " " << ddEdu0du1.p_[jdim*3+idim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdw0du1 " << jdim << " " << idim << " --> " << (dEdw0_[idim]-dEdw0[idim])/epsilon << " " << ddEdw0du1.p_[jdim*3+idim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdu1du1 " << jdim << " " << idim << " --> " << (dEdu1_[idim]-dEdu1[idim])/epsilon << " " << ddEdu1du1.p_[idim*3+jdim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << "*ddEdw1du1 " << jdim << " " << idim << " --> " << (dEdw1_[idim]-dEdw1[idim])/epsilon << " " << ddEdu1dw1.p_[idim*3+jdim] << std::endl;
    }
  }
  
  for(unsigned int jdim=0;jdim<3;jdim++){
    delfem2::CVec3d dEdu0_, dEdw0_, dEdu1_, dEdw1_;
    delfem2::CMat3d ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_;
    double energy_ = 0.0;
    delfem2::CVec3d w(0,0,0);
    w[jdim] = epsilon;
    delfem2::CMat3d dR = delfem2::Mat3_RotMatFromAxisAngleVec(w);
    delfem2::CMat3d R1_ = dR*R1;
    EdEddE_Joint(energy_, dEdu0_,dEdw0_, dEdu1_,dEdw1_,
                 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_,
                 ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_,
                 ddEdu1du1_, ddEdu1dw1_,
                 ddEdw1dw1_,
                 trans_stiff, rot_stiff, pj,
                 cg0, u0, R0,
                 cg1, u1, R1_);
    std::cout << "dEdw1: " << jdim << " -->  " << (energy_-energy)/epsilon << " " << dEdw1[jdim] << std::endl;
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdu0dw1 " << jdim << " " << idim << " --> " << (dEdu0_[idim]-dEdu0[idim])/epsilon << " " << ddEdu0dw1.p_[jdim*3+idim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdw0dw1 " << jdim << " " << idim << " --> " << (dEdw0_[idim]-dEdw0[idim])/epsilon << " " << ddEdw0dw1.p_[jdim*3+idim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdu1dw1 " << jdim << " " << idim << " --> " << (dEdu1_[idim]-dEdu1[idim])/epsilon << " " << ddEdu1dw1.p_[jdim*3+idim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdw1dw1 " << jdim << " " << idim << " --> " << (dEdw1_[idim]-dEdw1[idim])/epsilon << " " << ddEdw1dw1.p_[idim*3+jdim] << std::endl;
    }
  }
  
}




///////////////////////////////////////////////////////////////////////////


CRigidBodyAssembly_Static::CRigidBodyAssembly_Static()
{
  nitr = 30;
  damping_ratio = 0.01;
  ////
  n = delfem2::CVec3d(0,1,0); // normal direction of floor (should be an unit vector)
  gravity = delfem2::CVec3d(0,-10,0); // gravity
  ////
  cont_stiff = 1.0e+4; // contact_stiffness (insensitive)
  trans_stiff = 1.0e+4; // joint_translation_stiffness (insensitive)
  rot_stiff = 1.0e+9; // joint_rotation_stiffness (insensitive)

  scale_force = 0.10; //0.5;
  scale_torque = 0.05; //0.5;

  is_draw_deformed = false; //false
  is_draw_skeleton = false; //false;
  is_draw_force = false; //false;
  is_draw_section = true; //true;
  is_draw_grid = false; //true;
  is_draw_section_moment = false; //false;
  
//  SetExample();
//  this->Solve();
  this->is_draw_skeleton = true;
  this->is_draw_force = true;
  this->is_draw_grid = true;
}

CRigidBodyAssembly_Static::CRigidBodyAssembly_Static
(const std::vector<CRigidBody>& aRB,
 const std::vector<CJoint>& aJ)
{
  
  nitr = 30;
  damping_ratio = 0.01;
  ////
  n = delfem2::CVec3d(0,1,0); // normal direction of floor (should be an unit vector)
  gravity = delfem2::CVec3d(0,-10,0); // gravity
  ////
  cont_stiff = 1.0e+4; // contact_stiffness (insensitive)
  trans_stiff = 1.0e+4; // joint_translation_stiffness (insensitive)
  rot_stiff = 1.0e+9; // joint_rotation_stiffness (insensitive)
  
  scale_force = 0.10; //0.5;
  scale_torque = 0.05; //0.5;
  
  this->aRigidBody = aRB;
  this->aJoint = aJ;
//  this->Solve();
  this->is_draw_skeleton = true;
  this->is_draw_force = true;
  this->is_draw_grid = true;
}

void CRigidBodyAssembly_Static::AddRigidBody(const double centre_of_mass[3],
                  const double mass,
                  const std::vector<double>& contact_points)
{
  const std::vector<double> pcg(centre_of_mass,centre_of_mass+3);
  CRigidBody rb(mass, pcg);
  const int ncp = (int)contact_points.size()/3;
  for(int icp=0;icp<ncp;icp++) {
    delfem2::CVec3d pc(contact_points[icp*3+0],contact_points[icp*3+1],contact_points[icp*3+2]);
    rb.aCP.push_back(pc);
  }
  aRigidBody.push_back(rb);
}

void CRigidBodyAssembly_Static::AddJoint(const double position[3],
              const int body_index1,
              const int body_index2){
  CJoint j(body_index1, body_index2,
           std::vector<double>(position,position+3) );
//  j.p = delfem2::CVec3d(position[0], position[1], position[2]);
//  j.irb0 = body_index1;
//  j.irb1 = body_index2;
  aJoint.push_back(j);
}

void CRigidBodyAssembly_Static::ClearProblem(){
  aRigidBody.clear();
  aJoint.clear();
//  aPlate.clear();
}


void AddMatrix(Eigen::MatrixXd& M,
               unsigned int i0, unsigned int j0,
               const delfem2::CMat3d& m,
               bool isnt_inverse)
{
  if( isnt_inverse ){
    for(unsigned int idim=0;idim<3;idim++){
      for(unsigned int jdim=0;jdim<3;jdim++){
        M(i0+idim,j0+jdim) += m.p_[idim*3+jdim];
      }
    }
  }
  else{
    for(unsigned int idim=0;idim<3;idim++){
      for(unsigned int jdim=0;jdim<3;jdim++){
        M(i0+idim,j0+jdim) += m.p_[jdim*3+idim];
      }
    }
  }
}


void EdEddE_Total
(double& E,
 Eigen::VectorXd& dE,
 Eigen::MatrixXd& ddE,
 ////
 const std::vector<CRigidBody>& aRigidBody,
 const std::vector<CJoint>& aJoint,
 ////
 const delfem2::CVec3d n,
 const delfem2::CVec3d gravity,
 const double cont_stiff,
 const double trans_stiff,
 const double rot_stiff,
 bool is_friction
 )
{
  E = 0;
  ddE.setConstant(0);
  dE.setConstant(0); // u0,w0, u1,w1, ....
  
  ////
  for(unsigned int irb=0;irb<aRigidBody.size();irb++){
    const CRigidBody& rb = aRigidBody[irb];
    delfem2::CVec3d cg = rb.cg;
    {
      double e = 0;
      delfem2::CVec3d du, dw;
      EdEd_Potential(e, du, dw,
                       cg, rb.u, rb.m,
                       gravity);
      E += e;
      for(unsigned int idim=0;idim<3;idim++){
        dE[irb*6+0+idim] += du[idim];
        dE[irb*6+3+idim] += dw[idim];
      }
    }
    for(std::size_t icp=0;icp<rb.aCP.size();icp++){
      delfem2::CVec3d cp = rb.aCP[icp];
      double e = 0;
      delfem2::CVec3d du, dw;
      delfem2::CMat3d ddu,ddw,dudw;
      if( is_friction ){
        EdEddE_ContactFriction(e,
                                 du,dw,  ddu,ddw,dudw,
                                 cg,cp,rb.u,rb.R,
                                 cont_stiff);
      }
      else{
        EdEddE_Contact(e,
                         du,dw,  ddu,ddw,dudw,
                         cg,cp,rb.u,rb.R,
                         cont_stiff,n);
      }
      E += e;
      for(unsigned int idim=0;idim<3;idim++){
        dE[irb*6+0+idim] += du[idim];
        dE[irb*6+3+idim] += dw[idim];
      }
      for(unsigned int idim=0;idim<3;idim++){
        for(unsigned int jdim=0;jdim<3;jdim++){
          ddE(irb*6+0+idim,irb*6+0+jdim) += ddu.p_[idim*3+jdim];
        }
      }
      AddMatrix(ddE, irb*6+0, irb*6+0, ddu,  true );
      AddMatrix(ddE, irb*6+0, irb*6+3, dudw, false);
      //
      AddMatrix(ddE, irb*6+3, irb*6+0, dudw, true );
      AddMatrix(ddE, irb*6+3, irb*6+3, ddw,  true );
    }
    for(std::size_t iexf=0;iexf<rb.aExForce.size();iexf++){
      delfem2::CVec3d pex = rb.aExForce[iexf].first;
      delfem2::CVec3d fex = rb.aExForce[iexf].second;
      double e = 0;
      delfem2::CVec3d du, dw;
      delfem2::CMat3d ddu,ddw,dudw;
      EdEddE_Exforce(e,
                       du,dw,  ddu,ddw,dudw,
                       cg,pex,fex, rb.u,rb.R);
      E += e;
      for(int idim=0;idim<3;idim++){
        dE[irb*6+0+idim] += du[idim];
        dE[irb*6+3+idim] += dw[idim];
      }
      for(unsigned int idim=0;idim<3;idim++){
        for(unsigned int jdim=0;jdim<3;jdim++){
          ddE(irb*6+0+idim,irb*6+0+jdim) += ddu.p_[idim*3+jdim];
        }
      }
      AddMatrix(ddE, irb*6+0, irb*6+0, ddu,  true );
      AddMatrix(ddE, irb*6+0, irb*6+3, dudw, false);
      //
      AddMatrix(ddE, irb*6+3, irb*6+0, dudw, true );
      AddMatrix(ddE, irb*6+3, irb*6+3, ddw,  true );
    }
  }
  for(const auto & joint : aJoint){
    delfem2::CVec3d pj = joint.p;
    int irb0 = joint.irb0;
    int irb1 = joint.irb1;
    const CRigidBody& rb0 = aRigidBody[irb0];
    const CRigidBody& rb1 = aRigidBody[irb1];
    {
      double e = 0;
      delfem2::CVec3d du0, dw0, du1, dw1;
      delfem2::CMat3d du0du0, du0dw0, du0du1, du0dw1, dw0dw0, dw0du1, dw0dw1, du1du1, du1dw1, dw1dw1;
      EdEddE_Joint(e,
                     du0, dw0, du1, dw1,
                     du0du0, du0dw0, du0du1, du0dw1, dw0dw0, dw0du1, dw0dw1, du1du1, du1dw1, dw1dw1,
                     trans_stiff, rot_stiff, pj,
                     rb0.cg, rb0.u, rb0.R,
                     rb1.cg, rb1.u, rb1.R);
      E += e;
      for(int idim=0;idim<3;idim++){
        dE[irb0*6+0+idim] += du0[idim];
        dE[irb0*6+3+idim] += dw0[idim];
      }
      for(int idim=0;idim<3;idim++){
        dE[irb1*6+0+idim] += du1[idim];
        dE[irb1*6+3+idim] += dw1[idim];
      }
      AddMatrix(ddE, irb0*6+0, irb0*6+0, du0du0, true );
      AddMatrix(ddE, irb0*6+0, irb0*6+3, du0dw0, false);
      AddMatrix(ddE, irb0*6+0, irb1*6+0, du0du1, false);
      AddMatrix(ddE, irb0*6+0, irb1*6+3, du0dw1, false);
      ////
      AddMatrix(ddE, irb0*6+3, irb0*6+0, du0dw0, true );
      AddMatrix(ddE, irb0*6+3, irb0*6+3, dw0dw0, true );
      AddMatrix(ddE, irb0*6+3, irb1*6+0, dw0du1, false);
      AddMatrix(ddE, irb0*6+3, irb1*6+3, dw0dw1, false);
      ////
      AddMatrix(ddE, irb1*6+0, irb0*6+0, du0du1, true );
      AddMatrix(ddE, irb1*6+0, irb0*6+3, dw0du1, true );
      AddMatrix(ddE, irb1*6+0, irb1*6+0, du1du1, true );
      AddMatrix(ddE, irb1*6+0, irb1*6+3, du1dw1, false);
      ////
      AddMatrix(ddE, irb1*6+3, irb0*6+0, du0dw1, true );
      AddMatrix(ddE, irb1*6+3, irb0*6+3, dw0dw1, true );
      AddMatrix(ddE, irb1*6+3, irb1*6+0, du1dw1, true );
      AddMatrix(ddE, irb1*6+3, irb1*6+3, dw1dw1, true );
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

// solve one iteration
void CRigidBodyAssembly_Static::SolveOneIteration()
/*
(
 std::vector<CRigidBody>& aRigidBody,
 std::vector<CJoint>& aJoint,
 ////
 double damping_ratio,
 ////
 const delfem2::CVec3d n,
 const delfem2::CVec3d gravity,
 const double cont_stiff,
 const double trans_stiff,
 const double rot_stiff)*/
{
  int nDof = (int)aRigidBody.size()*6;
  double E = 0;
	Eigen::VectorXd dE(nDof);
	Eigen::MatrixXd ddE(nDof,nDof);
  EdEddE_Total(E,dE,ddE,
               aRigidBody,aJoint,
               n,gravity,cont_stiff,trans_stiff,rot_stiff,true);
  std::cout << "energy : " << E << std::endl;
  /////
	Eigen::VectorXd b(nDof);
	Eigen::MatrixXd A(nDof,nDof);
  b = -dE;
  A = ddE;
  for(int i=0;i<nDof;i++){
    A(i,i) += damping_ratio;
  }
	Eigen::VectorXd x = A.partialPivLu().solve(b);
  for(std::size_t irb=0;irb<aRigidBody.size();irb++){
    CRigidBody& rb = aRigidBody[irb];
    for(int idim=0;idim<3;idim++){
      rb.u[idim] += x(irb*6+0+idim);
    }
    {
      delfem2::CVec3d w;
      w.x = x(irb*6+3+0);
      w.y = x(irb*6+3+1);
      w.z = x(irb*6+3+2);
      delfem2::CMat3d dR = delfem2::Mat3_RotMatFromAxisAngleVec(w);
      rb.R = dR*rb.R;
    }
  }
}

void CRigidBodyAssembly_Static::Solve_InterPlane()
/*
(std::vector<CRigidBody>& aRigidBody,
 std::vector<CJoint>& aJoint,
 ////
 double damping_ratio,
 int nitr,
 ////
 const delfem2::CVec3d n,
 const delfem2::CVec3d gravity,
 const double cont_stiff,
 const double trans_stiff,
 const double rot_stiff)
 */
{
  for(int itr=0;itr<nitr;itr++){
    //SolveOneIteration(aRigidBody, aJoint, damping_ratio, n,gravity,cont_stiff,trans_stiff,rot_stiff);
      SolveOneIteration();
  }
  ComputeForces();
}

void CRigidBodyAssembly_Static::ComputeForces()
/*
(std::vector<CRigidBody>& aRigidBody,
 std::vector<CJoint>& aJoint,
 ////
 const delfem2::CVec3d n,
 const double cont_stiff, 
 const double trans_stiff,
 const double rot_stiff)
 */
{
  
  for(auto & irb : aRigidBody){
    CRigidBody& rb = irb;
    const int ncp = (int)irb.aCP.size();
    rb.aCForce.resize(ncp);
    for(int icp=0;icp<ncp;icp++){
      const delfem2::CVec3d& cp = irb.aCP[icp];
      delfem2::CVec3d Rv = rb.R * (cp-rb.cg);
      delfem2::CVec3d cq = Rv + rb.cg + rb.u;
      rb.aCForce[icp] = (cq.dot(n)*cont_stiff)*n;
    }
  }  
  for(auto & joint : aJoint){
    delfem2::CVec3d pj = joint.p;
    int irb0 = joint.irb0;
    int irb1 = joint.irb1;
    const CRigidBody& rb0 = aRigidBody[irb0];
    const CRigidBody& rb1 = aRigidBody[irb1];
    
    delfem2::CVec3d cg0 = rb0.cg;
    delfem2::CVec3d u0  = rb0.u;
    delfem2::CMat3d  R0  = rb0.R;
    
    delfem2::CVec3d cg1 = rb1.cg;
    delfem2::CVec3d u1  = rb1.u;
    delfem2::CMat3d  R1  = rb1.R;
    
    delfem2::CVec3d trans_f; // translation_force
    {
      delfem2::CVec3d Rv0 = R0 * (pj-cg0);
      delfem2::CVec3d qj0 = Rv0 + cg0 + u0; // after deformation joint pos relative to rigid body 0
      delfem2::CVec3d Rv1 = R1 * (pj-cg1);
      delfem2::CVec3d qj1 = Rv1 + cg1 + u1; // after deformation joint pos relative to rigid body 1
      trans_f = trans_stiff*(qj0 - qj1);
    }
    
    delfem2::CVec3d torque_f; // rotation_force
    {
      delfem2::CVec3d av(0,0,0);
      for(int i=0;i<3;i++){
        delfem2::CVec3d r0( R0.p_[0*3+i], R0.p_[1*3+i], R0.p_[2*3+i] );
        delfem2::CVec3d r1( R1.p_[0*3+i], R1.p_[1*3+i], R1.p_[2*3+i] );
        av += r0.cross(r1);
      }
      av *= 0.5;
      torque_f = rot_stiff*av;
    }
    joint.linear = trans_f;
    joint.torque = torque_f;
  }
}

//void CRigidBodyAssembly_Static::PrintJointForce()
/*
(const std::vector<CRigidBody>& aRigidBody,
 const std::vector<CJoint>& aJoint,
 ////
 const double trans_stiff,
 const double rot_stiff)
 */
/*
{
  std::cout << "force on joint" << std::endl;
    
  for(unsigned int ij=0;ij<aJoint.size();ij++){
    const CJoint& joint = aJoint[ij];
    delfem2::CVec3d pj = joint.p;
    int irb0 = joint.irb0;
    int irb1 = joint.irb1;
    const CRigidBody& rb0 = aRigidBody[irb0];
    const CRigidBody& rb1 = aRigidBody[irb1];
    
    delfem2::CVec3d cg0 = rb0.cg;
    delfem2::CVec3d u0  = rb0.u;
    delfem2::CMat3d  R0  = rb0.R;
    
    delfem2::CVec3d cg1 = rb1.cg;
    delfem2::CVec3d u1  = rb1.u;
    delfem2::CMat3d  R1  = rb1.R;
    
    delfem2::CVec3d trans_f; // translation_force
    {
      delfem2::CVec3d Rv0 = R0 * (pj-cg0);
      delfem2::CVec3d qj0 = Rv0 + cg0 + u0; // after deformation joint pos relative to rigid body 0
      delfem2::CVec3d Rv1 = R1 * (pj-cg1);
      delfem2::CVec3d qj1 = Rv1 + cg1 + u1; // after deformation joint pos relative to rigid body 1
      trans_f = trans_stiff*(qj0 - qj1);
    }
    
    delfem2::CVec3d torque_f; // rotation_force
    {
      delfem2::CVec3d av(0,0,0);
      for(unsigned int i=0;i<3;i++){
        delfem2::CVec3d r0( R0.mat[0*3+i], R0.mat[1*3+i], R0.mat[2*3+i] );
        delfem2::CVec3d r1( R1.mat[0*3+i], R1.mat[1*3+i], R1.mat[2*3+i] );
        av += (r0^r1);
      }
      av *= 0.5;
      torque_f = rot_stiff*av;
    }
    
    std::cout << "force of joint: " << ij << std::endl;
    std::cout << "  trans_force:  " << trans_f.x << " " << trans_f.y << " " << trans_f.z << std::endl;
    std::cout << "  torque_force: " << torque_f.x << " " << torque_f.y << " " << torque_f.z << std::endl;
  }
}
 */



// Setting problem here
void CRigidBodyAssembly_Static::SetExample()
{
  aRigidBody.clear();
  aJoint.clear();
  
  { // making plane 0
    CRigidBody rb(1.0, delfem2::CVec3d(0,1,0).stlvec() );
    aRigidBody.push_back(rb);
  }
  { // making plane 1
    CRigidBody rb(0.1, delfem2::CVec3d(-1, 0.5, -1).stlvec());
    rb.addCP(delfem2::CVec3d(-1.1,0,-1.1).stlvec());
    aRigidBody.push_back(rb);
  }
  { // making plane 2
    CRigidBody rb(0.1, delfem2::CVec3d(-1,0.5,+1).stlvec() );
    rb.addCP( delfem2::CVec3d(-1.1,0,+1.1).stlvec() );
    aRigidBody.push_back(rb);
  }
  { // making plane 3
    CRigidBody rb(0.1, delfem2::CVec3d(+1,0.5,-1).stlvec() );
    rb.addCP( delfem2::CVec3d(+1.1,0,-1.1).stlvec() );
    aRigidBody.push_back(rb);
  }
  { // making plane 4
    CRigidBody rb(0.1, delfem2::CVec3d(+1,0.5,+1).stlvec() );
    rb.addCP( delfem2::CVec3d(+1.1,0,+1.1).stlvec() );
    aRigidBody.push_back(rb);
  }
  
  
  { // joint 0
    CJoint jt(0,1, delfem2::CVec3d(-1,+1,-1).stlvec() );
    aJoint.push_back(jt);
  }
  { // joint 1
    CJoint jt(0,2, delfem2::CVec3d(-1,+1,+1).stlvec() );
    aJoint.push_back(jt);
  }
  { // joint 2
    CJoint jt(0,3, delfem2::CVec3d(+1,+1,-1).stlvec() );
    aJoint.push_back(jt);
  }
  { // joint 3
    CJoint jt(0,4, delfem2::CVec3d(+1,+1,+1).stlvec() );
    aJoint.push_back(jt);
  }
}


static void myGlVertex3d(const delfem2::CVec3d& v){
  ::glVertex3d(v.x,v.y,v.z);
}

std::vector<double> CRigidBodyAssembly_Static::MinMaxXYZ() const
{
  double bb[6] = {1,-1, 0,0, 0,0};
  for(const auto & rb : aRigidBody){
    myUpdateMinMaxXYZ(bb, rb.cg);
  }
  for(const auto & j : aJoint){
    myUpdateMinMaxXYZ(bb, j.p);
  }
  std::vector<double> res(bb,bb+6);
  return res;
}



