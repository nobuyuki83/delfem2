/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/v23m3q.h"


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

#include "delfem2/eigen/eigen_rigidbody.h"

void myUpdateMinMaxXYZ
(double bb[6],
 const CVector3 p)
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

void EdEd_Potential
(double& energy,
 CVector3& dEdu,
 CVector3& dEdw,
 const CVector3& cg,
 const CVector3& u,
 const double mass,
 const CVector3& g       // floor normal
)
{
  energy = mass*(u*g);
  dEdu = mass*g;
  dEdw = CVector3(0,0,0);
}

void EdEddE_Exforce
(double& energy,
 CVector3& dEdu, CVector3& dEdw,
 CMatrix3&  ddEddu,  CMatrix3&  ddEddw, CMatrix3&  ddEdudw,
 const CVector3& cg, // the center of gravity
 const CVector3& pex, // external force position
 const CVector3& fex, // external force
 const CVector3& u, // displacement
 const CMatrix3& R // rigid rotation
)
{
  CVector3 Rv = R * (pex-cg);
  CVector3 qex = Rv + cg + u; // current external force position
  energy = +qex*fex;
  dEdu = +fex;
  dEdw = +Rv^fex;
  ddEddu  = CMatrix3(0.0);
  ddEddw  = +Mat3_Spin(fex)*Mat3_Spin(Rv);
  ddEdudw = CMatrix3(0.0);
}

void EdEddE_Contact
(double& energy,
 CVector3& dEdu, CVector3& dEdw,
 CMatrix3&  ddEddu,  CMatrix3&  ddEddw, CMatrix3&  ddEdudw,
 const CVector3& cg, // the center of gravity
 const CVector3& cp, // contact position
 const CVector3& u, // displacement
 const CMatrix3& R, // rigid rotation
 const double cont_stiff,
 const CVector3& n       // floor normal
)
{
  CVector3 Rv = R * (cp-cg);
  CVector3 cq = Rv + cg + u;
  energy = 0.5*(cq*n)*(cq*n)*cont_stiff;
  dEdu = ((cq*n)*cont_stiff)*n;
  dEdw = ((cq*n)*cont_stiff)*(Rv^n);
  ddEddu  = cont_stiff*Mat3_OuterProduct(n,n);
  ddEddw  = cont_stiff*Mat3_OuterProduct(Rv^n,Rv^n) + ((cq*n)*cont_stiff)*Mat3_Spin(n)*Mat3_Spin(Rv);
  ddEdudw = cont_stiff*Mat3_OuterProduct(Rv^n,n);
}

void EdEddE_ContactFriction
(double& energy,
 CVector3& dEdu, CVector3& dEdw,
 CMatrix3&  ddEddu,  CMatrix3&  ddEddw, CMatrix3&  ddEdudw,
 const CVector3& cg, // the center of gravity
 const CVector3& cp, // contact position
 const CVector3& u, // displacement
 const CMatrix3& R, // rigid rotation
 const double cont_stiff
 )
{
  CVector3 Rv = R * (cp-cg);
  CVector3 cq = Rv + cg + u;
  energy = 0.5*(cq-cp)*(cq-cp)*cont_stiff;
  dEdu = cont_stiff*(cq-cp);
  dEdw = cont_stiff*(Rv^(cq-cp));
  ddEddu  = cont_stiff*CMatrix3::Identity();
  //    ddEddw  = -cont_stiff*CMatrix3::OuterProduct(Rv,Rv) + cont_stiff*CMatrix3::Spin(Rv)*CMatrix3::Spin(cq-cp);
  ddEddw  = -cont_stiff*Mat3_Spin(Rv)*Mat3_Spin(Rv) + cont_stiff*Mat3_Spin(cq-cp)*Mat3_Spin(Rv);
  ddEdudw = cont_stiff*Mat3_Spin(Rv);
}


void EdEddE_Joint
(double& energy,
 CVector3& dEdu0, CVector3& dEdw0, CVector3& dEdu1, CVector3& dEdw1,
 ////
 CMatrix3& ddEdu0du0,  CMatrix3& ddEdu0dw0, CMatrix3& ddEdu0du1, CMatrix3& ddEdu0dw1, CMatrix3& ddEdw0dw0,
 CMatrix3& ddEdw0du1,  CMatrix3& ddEdw0dw1, CMatrix3& ddEdu1du1, CMatrix3& ddEdu1dw1, CMatrix3& ddEdw1dw1,
 ////
 const double trans_stiff,
 const double rot_stiff,
 const CVector3& pj,
 ////
 const CVector3& cg0,  const CVector3& u0,  const CMatrix3& R0,
 const CVector3& cg1,  const CVector3& u1,  const CMatrix3& R1
 )
{
  CVector3 Rv0 = R0 * (pj-cg0);
  CVector3 qj0 = Rv0 + cg0 + u0; // after deformation joint pos relative to rigid body 0
  
  CVector3 Rv1 = R1 * (pj-cg1);
  CVector3 qj1 = Rv1 + cg1 + u1; // after deformation joint pos relative to rigid body 1
  
  energy = 0.5*trans_stiff*(qj0 - qj1).DLength();
  
  dEdu0 = trans_stiff*(qj0-qj1);
  dEdw0 = trans_stiff*(Rv0^(qj0-qj1));
  
  dEdu1 = trans_stiff*(qj1-qj0);
  dEdw1 = trans_stiff*(Rv1^(qj1-qj0));
  
  ddEdu0du0 =  trans_stiff*CMatrix3::Identity();
  ddEdu0du1 = -trans_stiff*CMatrix3::Identity();
  ddEdu0dw0 = +trans_stiff*Mat3_Spin(Rv0);
  ddEdu0dw1 = -trans_stiff*Mat3_Spin(Rv1);
  
  ddEdw0dw0 = -trans_stiff*Mat3_Spin(Rv0)*Mat3_Spin(Rv0) + trans_stiff*Mat3_Spin(qj0-qj1)*Mat3_Spin(Rv0);
  ddEdw0du1 = +trans_stiff*Mat3_Spin(Rv0);
  ddEdw0dw1 = +trans_stiff*Mat3_Spin(Rv1)*Mat3_Spin(Rv0);
  
  ddEdu1du1 =  trans_stiff*CMatrix3::Identity();
  ddEdu1dw1 = +trans_stiff*Mat3_Spin(Rv1);
  
  ddEdw1dw1 = -trans_stiff*Mat3_Spin(Rv1)*Mat3_Spin(Rv1) + trans_stiff*Mat3_Spin(qj1-qj0)*Mat3_Spin(Rv1);
  
  CVector3 av(0,0,0);
  CMatrix3 davdw0, davdw1;
  {
    for(unsigned int i=0;i<3;i++){
      CVector3 r0( R0.mat[0*3+i], R0.mat[1*3+i], R0.mat[2*3+i] );
      CVector3 r1( R1.mat[0*3+i], R1.mat[1*3+i], R1.mat[2*3+i] );
      av += (r0^r1);
      davdw0 += Mat3_Spin(r1)*Mat3_Spin(r0);
      davdw1 -= Mat3_Spin(r0)*Mat3_Spin(r1);
    }
    av *= 0.5;
    davdw0 *= 0.5;
    davdw1 *= 0.5;
  }
  energy += 0.5*rot_stiff*av.DLength();
  
  CMatrix3 m0,m1,m2;
  for(unsigned int i=0;i<3;i++){
    CVector3 r0( R0.mat[0*3+i], R0.mat[1*3+i], R0.mat[2*3+i] );
    CVector3 r1( R1.mat[0*3+i], R1.mat[1*3+i], R1.mat[2*3+i] );
    dEdw0 += 0.5*rot_stiff*r0^(r1^av);
    dEdw1 -= 0.5*rot_stiff*r1^(r0^av);
    ddEdw0dw0 += 0.5*rot_stiff*( Mat3_Spin(r1^av)*Mat3_Spin(r0) + Mat3_Spin(r0)*Mat3_Spin(r1)*davdw0 );
    ddEdw1dw1 -= 0.5*rot_stiff*( Mat3_Spin(r0^av)*Mat3_Spin(r1) + Mat3_Spin(r1)*Mat3_Spin(r0)*davdw1 );
    ddEdw0dw1 -= 0.5*rot_stiff*(  Mat3_Spin(r1)*Mat3_Spin(av)*Mat3_Spin(r0)
                                + Mat3_Spin(r1)*Mat3_Spin(r0)*davdw0 );
  }
}


CVector3 rand_vec(double s)
{
  CVector3 v;
  v.x = s*rand()/(RAND_MAX+1.0);
  v.y = s*rand()/(RAND_MAX+1.0);
  v.z = s*rand()/(RAND_MAX+1.0);
  return v;
}

CMatrix3 rand_rot()
{
  double s = 3.5;
  CVector3 v;
  v.x = s*rand()/(RAND_MAX+1.0);
  v.y = s*rand()/(RAND_MAX+1.0);
  v.z = s*rand()/(RAND_MAX+1.0);
  CMatrix3 R;  R.SetRotMatrix_Cartesian(v.x,v.y,v.z);
  return R;
}

void CheckDiff_Contact()
{
  double energy = 0.0;
  CVector3 dEdu,dEdw;
  CMatrix3 ddEddu, ddEddw, ddEdudw;
  double epsilon = 1.0e-5;
  const double cont_stiff = 1.0e+3;
  CVector3 cg = rand_vec(1.0);
  CVector3 cp = rand_vec(1.0);
  CVector3 u  = rand_vec(1.0);
  CMatrix3  R  = rand_rot();
  CVector3 n(0,0,1);
  EdEddE_Contact(energy,
                 dEdu,dEdw,  ddEddu,ddEddw,ddEdudw,
                 cg,cp,u,R,
                 cont_stiff,n);
  for(unsigned int idim=0;idim<3;idim++){
    CVector3 dEdu_, dEdw_;
    CMatrix3 ddEddu_, ddEddw_, ddEdudw_;
    double energy_ = 0.0;
    CVector3 u_ = u;
    u_[idim] += epsilon;
    EdEddE_Contact(energy_,  dEdu_,dEdw_,   ddEddu_,ddEddw_,ddEdudw_,
                   cg,cp,u_,R,
                   cont_stiff,n);
    std::cout << "dEdu " << idim << " --> "  << (energy_-energy)/epsilon << " " << dEdu[idim] << std::endl;
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEddu  " << idim << " " << jdim << " --> " << (dEdu_[jdim]-dEdu[jdim])/epsilon << " " << ddEddu.mat[jdim*3+idim] << std::endl;
    }
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEdudw " << idim << " " << jdim << " --> " << (dEdw_[jdim]-dEdw[jdim])/epsilon << " " << ddEdudw.mat[jdim*3+idim] << std::endl;
    }
  }
  for(unsigned int idim=0;idim<3;idim++){
    CVector3 dEdu_, dEdw_;
    CMatrix3 ddEddu_, ddEddw_, ddEdudw_;
    double energy_ = 0.0;
    CVector3 w(0,0,0);
    w[idim] = epsilon;
    CMatrix3 dR; dR.SetRotMatrix_Cartesian(w.x,w.y,w.z);
    CMatrix3 R_ = dR*R;
    EdEddE_Contact(energy_,  dEdu_,dEdw_,   ddEddu_,ddEddw_,ddEdudw_,
                   cg,cp,u,R_,
                   cont_stiff,n);
    std::cout << "dEdw " << idim << " --> " << (energy_-energy)/epsilon << " " << dEdw[idim] << std::endl;
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEdudw " << idim << " " << jdim << " --> " << (dEdu_[jdim]-dEdu[jdim])/epsilon << " " << ddEdudw.mat[idim*3+jdim] << std::endl;
    }
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEddw  " << idim << " " << jdim << " --> " << (dEdw_[jdim]-dEdw[jdim])/epsilon << " " << ddEddw.mat[jdim*3+idim] << std::endl;
    }
  }
}


void CheckDiff_ContactFriction()
{
  double energy = 0.0;
  CVector3 dEdu,dEdw;
  CMatrix3 ddEddu, ddEddw, ddEdudw;
  double epsilon = 1.0e-5;
  const double cont_stiff = 1.0e+3;
  CVector3 cg = rand_vec(1.0);
  CVector3 cp = rand_vec(1.0);
  CVector3 u  = rand_vec(1.0);
  CMatrix3  R  = rand_rot();
  EdEddE_ContactFriction(energy,
                         dEdu,dEdw,  ddEddu,ddEddw,ddEdudw,
                         cg,cp,u,R,
                         cont_stiff);
  for(unsigned int idim=0;idim<3;idim++){
    CVector3 dEdu_, dEdw_;
    CMatrix3 ddEddu_, ddEddw_, ddEdudw_;
    double energy_ = 0.0;
    CVector3 u_ = u;
    u_[idim] += epsilon;
    EdEddE_ContactFriction(energy_,  dEdu_,dEdw_,   ddEddu_,ddEddw_,ddEdudw_,
                           cg,cp,u_,R,
                           cont_stiff);
    std::cout << "dEdu " << idim << " --> "  << (energy_-energy)/epsilon << " " << dEdu[idim] << std::endl;
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEddu  " << idim << " " << jdim << " --> " << (dEdu_[jdim]-dEdu[jdim])/epsilon << " " << ddEddu.mat[jdim*3+idim] << std::endl;
    }
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEdudw " << idim << " " << jdim << " --> " << (dEdw_[jdim]-dEdw[jdim])/epsilon << " " << ddEdudw.mat[jdim*3+idim] << std::endl;
    }
  }
  for(unsigned int idim=0;idim<3;idim++){
    CVector3 dEdu_, dEdw_;
    CMatrix3 ddEddu_, ddEddw_, ddEdudw_;
    double energy_ = 0.0;
    CVector3 w(0,0,0);
    w[idim] = epsilon;
    CMatrix3 dR; dR.SetRotMatrix_Cartesian(w.x,w.y,w.z);
    CMatrix3 R_ = dR*R;
    EdEddE_ContactFriction(energy_,  dEdu_,dEdw_,   ddEddu_,ddEddw_,ddEdudw_,
                           cg,cp,u,R_,
                           cont_stiff);
    std::cout << "dEdw " << idim << " --> " << (energy_-energy)/epsilon << " " << dEdw[idim] << std::endl;
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEdudw " << idim << " " << jdim << " --> " << (dEdu_[jdim]-dEdu[jdim])/epsilon << " " << ddEdudw.mat[idim*3+jdim] << std::endl;
    }
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEddw  " << idim << " " << jdim << " --> " << (dEdw_[jdim]-dEdw[jdim])/epsilon << " " << ddEddw.mat[jdim*3+idim] << std::endl;
    }
    
  }
}



void CheckDiff_Exforce()
{
  double energy = 0.0;
  CVector3 dEdu,dEdw;
  CMatrix3 ddEddu, ddEddw, ddEdudw;
  double epsilon = 1.0e-5;
  CVector3 cg = rand_vec(1.0);
  CVector3 u  = rand_vec(1.0);
  CVector3 fex = rand_vec(1.0);
  CVector3 pex = rand_vec(1.0);
  CMatrix3  R  = rand_rot();
  CVector3 n(0,0,1);
  EdEddE_Exforce(energy,
                 dEdu,dEdw,  ddEddu,ddEddw,ddEdudw,
                 cg,pex,fex,u,R);
  for(unsigned int idim=0;idim<3;idim++){
    CVector3 dEdu_, dEdw_;
    CMatrix3 ddEddu_, ddEddw_, ddEdudw_;
    double energy_ = 0.0;
    CVector3 u_ = u;
    u_[idim] += epsilon;
    EdEddE_Exforce(energy_,  dEdu_,dEdw_,   ddEddu_,ddEddw_,ddEdudw_,
                   cg,pex,fex,u_,R);
    std::cout << "dEdu " << idim << " --> "  << (energy_-energy)/epsilon << " " << dEdu[idim] << std::endl;
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEddu  " << idim << " " << jdim << " --> " << (dEdu_[jdim]-dEdu[jdim])/epsilon << " " << ddEddu.mat[jdim*3+idim] << std::endl;
    }
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEdudw " << idim << " " << jdim << " --> " << (dEdw_[jdim]-dEdw[jdim])/epsilon << " " << ddEdudw.mat[jdim*3+idim] << std::endl;
    }
  }
  for(unsigned int idim=0;idim<3;idim++){
    CVector3 dEdu_, dEdw_;
    CMatrix3 ddEddu_, ddEddw_, ddEdudw_;
    double energy_ = 0.0;
    CVector3 w(0,0,0);
    w[idim] = epsilon;
    CMatrix3 dR; dR.SetRotMatrix_Cartesian(w.x,w.y,w.z);
    CMatrix3 R_ = dR*R;
    EdEddE_Exforce(energy_,  dEdu_,dEdw_,   ddEddu_,ddEddw_,ddEdudw_,
                   cg,pex,fex,u,R_);
    std::cout << "dEdw " << idim << " --> " << (energy_-energy)/epsilon << " " << dEdw[idim] << std::endl;
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEdudw " << idim << " " << jdim << " --> " << (dEdu_[jdim]-dEdu[jdim])/epsilon << " " << ddEdudw.mat[idim*3+jdim] << std::endl;
    }
    for(unsigned int jdim=0;jdim<3;jdim++){
      std::cout << " ddEddw  " << idim << " " << jdim << " --> " << (dEdw_[jdim]-dEdw[jdim])/epsilon << " " << ddEddw.mat[jdim*3+idim] << std::endl;
    }
  }
}




void CheckDiff_Joint()
{
  double energy = 0.0;
  CVector3 dEdu0(0,0,0);
  CVector3 dEdw0(0,0,0);
  CVector3 dEdu1(0,0,0);
  CVector3 dEdw1(0,0,0);
  CMatrix3 ddEdu0du0;
  CMatrix3 ddEdu0dw0;
  CMatrix3 ddEdu0du1;
  CMatrix3 ddEdu0dw1;
  CMatrix3 ddEdw0dw0;
  CMatrix3 ddEdw0du1;
  CMatrix3 ddEdw0dw1;
  CMatrix3 ddEdu1du1;
  CMatrix3 ddEdu1dw1;
  CMatrix3 ddEdw1dw1;
  double epsilon = 1.0e-5;
  const double trans_stiff = 1.0e+3;
  const double rot_stiff = 1.0e+3;
  CVector3 pj  = rand_vec(1.0);
  ////
  CVector3 cg0 = rand_vec(1.0);
  CVector3 cp0 = rand_vec(1.0);
  CVector3 u0  = rand_vec(1.0);
  CMatrix3  R0  = rand_rot();
  ////
  CVector3 cg1 = rand_vec(1.0);
  CVector3 cp1 = rand_vec(1.0);
  CVector3 u1  = rand_vec(1.0);
  CMatrix3  R1  = rand_rot();
  ////
  EdEddE_Joint(energy,
               dEdu0, dEdw0, dEdu1,dEdw1,
               ddEdu0du0, ddEdu0dw0, ddEdu0du1, ddEdu0dw1, ddEdw0dw0, ddEdw0du1, ddEdw0dw1, ddEdu1du1, ddEdu1dw1, ddEdw1dw1,
               trans_stiff, rot_stiff, pj,
               cg0, u0, R0,
               cg1, u1, R1);
  for(unsigned int kdim=0;kdim<3;kdim++){
    CVector3 dEdu0_, dEdw0_, dEdu1_, dEdw1_;
    CMatrix3 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_;
    double energy_ = 0.0;
    CVector3 u0_ = u0;
    u0_[kdim] += epsilon;
    EdEddE_Joint(energy_, dEdu0_,dEdw0_, dEdu1_,dEdw1_,
                 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_,
                 trans_stiff, rot_stiff, pj,
                 cg0, u0_, R0,
                 cg1, u1,  R1);
    std::cout << "dEdu0: " << kdim << " -->  " << (energy_-energy)/epsilon << " " << dEdu0[kdim] << std::endl;
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdu0du0 " << kdim << " " << idim << " --> " << (dEdu0_[idim]-dEdu0[idim])/epsilon << " " << ddEdu0du0.mat[idim*3+kdim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdw0du0 " << kdim << " " << idim << " --> " << (dEdw0_[idim]-dEdw0[idim])/epsilon << " " << ddEdu0dw0.mat[idim*3+kdim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << "*ddEdu1du0 " << kdim << " " << idim << " --> " << (dEdu1_[idim]-dEdu1[idim])/epsilon << " " << ddEdu0du1.mat[idim*3+kdim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << "*ddEdw1du0 " << kdim << " " << idim << " --> " << (dEdw1_[idim]-dEdw1[idim])/epsilon << " " << ddEdu0dw1.mat[idim*3+kdim] << std::endl;
    }
  }
  for(unsigned int jdim=0;jdim<3;jdim++){
    CVector3 dEdu0_, dEdw0_, dEdu1_, dEdw1_;
    CMatrix3 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_;
    double energy_ = 0.0;
    CVector3 w(0,0,0);
    w[jdim] = epsilon;
    CMatrix3 dR; dR.SetRotMatrix_Cartesian(w.x,w.y,w.z);
    CMatrix3 R0_ = dR*R0;
    EdEddE_Joint(energy_, dEdu0_,dEdw0_, dEdu1_,dEdw1_,
                 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_,
                 trans_stiff, rot_stiff, pj,
                 cg0, u0, R0_,
                 cg1, u1, R1);
    std::cout << "dEdw0: " << jdim << " -->  " << (energy_-energy)/epsilon << " " << dEdw0[jdim] << std::endl;
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdu0dw0 " << jdim << " " << idim << " --> " << (dEdu0_[idim]-dEdu0[idim])/epsilon << " " << ddEdu0dw0.mat[jdim*3+idim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdw0dw0 " << jdim << " " << idim << " --> " << (dEdw0_[idim]-dEdw0[idim])/epsilon << " " << ddEdw0dw0.mat[idim*3+jdim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdu1dw0 " << jdim << " " << idim << " --> " << (dEdu1_[idim]-dEdu1[idim])/epsilon << " " << ddEdw0du1.mat[idim*3+jdim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdw1dw0 " << jdim << " " << idim << " --> " << (dEdw1_[idim]-dEdw1[idim])/epsilon << " " << ddEdw0dw1.mat[idim*3+jdim] << std::endl;
    }
  }
  
  for(unsigned int jdim=0;jdim<3;jdim++){
    CVector3 dEdu0_, dEdw0_, dEdu1_, dEdw1_;
    CMatrix3 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_;
    double energy_ = 0.0;
    CVector3 u1_ = u1;
    u1_[jdim] += epsilon;
    EdEddE_Joint(energy_, dEdu0_,dEdw0_, dEdu1_,dEdw1_,
                 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_,
                 trans_stiff, rot_stiff, pj,
                 cg0, u0, R0,
                 cg1, u1_,R1);
    std::cout << "dEdu1: " << jdim << " -->  " << (energy_-energy)/epsilon << " " << dEdu1[jdim] << std::endl;
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdu0du1 " << jdim << " " << idim << " --> " << (dEdu0_[idim]-dEdu0[idim])/epsilon << " " << ddEdu0du1.mat[jdim*3+idim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdw0du1 " << jdim << " " << idim << " --> " << (dEdw0_[idim]-dEdw0[idim])/epsilon << " " << ddEdw0du1.mat[jdim*3+idim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdu1du1 " << jdim << " " << idim << " --> " << (dEdu1_[idim]-dEdu1[idim])/epsilon << " " << ddEdu1du1.mat[idim*3+jdim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << "*ddEdw1du1 " << jdim << " " << idim << " --> " << (dEdw1_[idim]-dEdw1[idim])/epsilon << " " << ddEdu1dw1.mat[idim*3+jdim] << std::endl;
    }
  }
  
  for(unsigned int jdim=0;jdim<3;jdim++){
    CVector3 dEdu0_, dEdw0_, dEdu1_, dEdw1_;
    CMatrix3 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_;
    double energy_ = 0.0;
    CVector3 w(0,0,0);
    w[jdim] = epsilon;
    CMatrix3 dR; dR.SetRotMatrix_Cartesian(w.x,w.y,w.z);
    CMatrix3 R1_ = dR*R1;
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
      std::cout << " ddEdu0dw1 " << jdim << " " << idim << " --> " << (dEdu0_[idim]-dEdu0[idim])/epsilon << " " << ddEdu0dw1.mat[jdim*3+idim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdw0dw1 " << jdim << " " << idim << " --> " << (dEdw0_[idim]-dEdw0[idim])/epsilon << " " << ddEdw0dw1.mat[jdim*3+idim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdu1dw1 " << jdim << " " << idim << " --> " << (dEdu1_[idim]-dEdu1[idim])/epsilon << " " << ddEdu1dw1.mat[jdim*3+idim] << std::endl;
    }
    for(unsigned int idim=0;idim<3;idim++){
      std::cout << " ddEdw1dw1 " << jdim << " " << idim << " --> " << (dEdw1_[idim]-dEdw1[idim])/epsilon << " " << ddEdw1dw1.mat[idim*3+jdim] << std::endl;
    }
  }
  
}




///////////////////////////////////////////////////////////////////////////


CRigidBodyAssembly_Static::CRigidBodyAssembly_Static()
{
  nitr = 30;
  damping_ratio = 0.01;
  ////
  n = CVector3(0,1,0); // normal direction of floor (should be an unit vector)
  gravity = CVector3(0,-10,0); // gravity
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
  n = CVector3(0,1,0); // normal direction of floor (should be an unit vector)
  gravity = CVector3(0,-10,0); // gravity
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
    CVector3 pc(contact_points[icp*3+0],contact_points[icp*3+1],contact_points[icp*3+2]);
    rb.aCP.push_back(pc);
  }
  aRigidBody.push_back(rb);
}

void CRigidBodyAssembly_Static::AddJoint(const double position[3],
              const int body_index1,
              const int body_index2){
  CJoint j(body_index1, body_index2,
           std::vector<double>(position,position+3) );
//  j.p = CVector3(position[0], position[1], position[2]);
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
               const CMatrix3& m,
               bool isnt_inverse)
{
  if( isnt_inverse ){
    for(unsigned int idim=0;idim<3;idim++){
      for(unsigned int jdim=0;jdim<3;jdim++){
        M(i0+idim,j0+jdim) += m.mat[idim*3+jdim];
      }
    }
  }
  else{
    for(unsigned int idim=0;idim<3;idim++){
      for(unsigned int jdim=0;jdim<3;jdim++){
        M(i0+idim,j0+jdim) += m.mat[jdim*3+idim];
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
 const CVector3 n,
 const CVector3 gravity,
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
    CVector3 cg = rb.cg;
    {
      double e = 0;
      CVector3 du, dw;
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
      CVector3 cp = rb.aCP[icp];
      double e = 0;
      CVector3 du, dw;
      CMatrix3 ddu,ddw,dudw;
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
          ddE(irb*6+0+idim,irb*6+0+jdim) += ddu.mat[idim*3+jdim];
        }
      }
      AddMatrix(ddE, irb*6+0, irb*6+0, ddu,  true );
      AddMatrix(ddE, irb*6+0, irb*6+3, dudw, false);
      ////
      AddMatrix(ddE, irb*6+3, irb*6+0, dudw, true );
      AddMatrix(ddE, irb*6+3, irb*6+3, ddw,  true );
    }
    for(std::size_t iexf=0;iexf<rb.aExForce.size();iexf++){
      CVector3 pex = rb.aExForce[iexf].first;
      CVector3 fex = rb.aExForce[iexf].second;
      double e = 0;
      CVector3 du, dw;
      CMatrix3 ddu,ddw,dudw;
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
          ddE(irb*6+0+idim,irb*6+0+jdim) += ddu.mat[idim*3+jdim];
        }
      }
      AddMatrix(ddE, irb*6+0, irb*6+0, ddu,  true );
      AddMatrix(ddE, irb*6+0, irb*6+3, dudw, false);
      ////
      AddMatrix(ddE, irb*6+3, irb*6+0, dudw, true );
      AddMatrix(ddE, irb*6+3, irb*6+3, ddw,  true );
    }
  }
  for(const auto & joint : aJoint){
    CVector3 pj = joint.p;
    int irb0 = joint.irb0;
    int irb1 = joint.irb1;
    const CRigidBody& rb0 = aRigidBody[irb0];
    const CRigidBody& rb1 = aRigidBody[irb1];
    {
      double e = 0;
      CVector3 du0, dw0, du1, dw1;
      CMatrix3 du0du0, du0dw0, du0du1, du0dw1, dw0dw0, dw0du1, dw0dw1, du1du1, du1dw1, dw1dw1;
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
 const CVector3 n,
 const CVector3 gravity,
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
      CVector3 w;
      w.x = x(irb*6+3+0);
      w.y = x(irb*6+3+1);
      w.z = x(irb*6+3+2);
      CMatrix3 dR; dR.SetRotMatrix_Cartesian(w.x,w.y,w.z);
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
 const CVector3 n,
 const CVector3 gravity,
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
 const CVector3 n,
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
      const CVector3& cp = irb.aCP[icp];
      CVector3 Rv = rb.R * (cp-rb.cg);
      CVector3 cq = Rv + rb.cg + rb.u;
      rb.aCForce[icp] = ((cq*n)*cont_stiff)*n;
    }
  }  
  for(auto & joint : aJoint){
    CVector3 pj = joint.p;
    int irb0 = joint.irb0;
    int irb1 = joint.irb1;
    const CRigidBody& rb0 = aRigidBody[irb0];
    const CRigidBody& rb1 = aRigidBody[irb1];
    
    CVector3 cg0 = rb0.cg;
    CVector3 u0  = rb0.u;
    CMatrix3  R0  = rb0.R;
    
    CVector3 cg1 = rb1.cg;
    CVector3 u1  = rb1.u;
    CMatrix3  R1  = rb1.R;
    
    CVector3 trans_f; // translation_force
    {
      CVector3 Rv0 = R0 * (pj-cg0);
      CVector3 qj0 = Rv0 + cg0 + u0; // after deformation joint pos relative to rigid body 0
      CVector3 Rv1 = R1 * (pj-cg1);
      CVector3 qj1 = Rv1 + cg1 + u1; // after deformation joint pos relative to rigid body 1
      trans_f = trans_stiff*(qj0 - qj1);
    }
    
    CVector3 torque_f; // rotation_force
    {
      CVector3 av(0,0,0);
      for(int i=0;i<3;i++){
        CVector3 r0( R0.mat[0*3+i], R0.mat[1*3+i], R0.mat[2*3+i] );
        CVector3 r1( R1.mat[0*3+i], R1.mat[1*3+i], R1.mat[2*3+i] );
        av += (r0^r1);
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
    CVector3 pj = joint.p;
    int irb0 = joint.irb0;
    int irb1 = joint.irb1;
    const CRigidBody& rb0 = aRigidBody[irb0];
    const CRigidBody& rb1 = aRigidBody[irb1];
    
    CVector3 cg0 = rb0.cg;
    CVector3 u0  = rb0.u;
    CMatrix3  R0  = rb0.R;
    
    CVector3 cg1 = rb1.cg;
    CVector3 u1  = rb1.u;
    CMatrix3  R1  = rb1.R;
    
    CVector3 trans_f; // translation_force
    {
      CVector3 Rv0 = R0 * (pj-cg0);
      CVector3 qj0 = Rv0 + cg0 + u0; // after deformation joint pos relative to rigid body 0
      CVector3 Rv1 = R1 * (pj-cg1);
      CVector3 qj1 = Rv1 + cg1 + u1; // after deformation joint pos relative to rigid body 1
      trans_f = trans_stiff*(qj0 - qj1);
    }
    
    CVector3 torque_f; // rotation_force
    {
      CVector3 av(0,0,0);
      for(unsigned int i=0;i<3;i++){
        CVector3 r0( R0.mat[0*3+i], R0.mat[1*3+i], R0.mat[2*3+i] );
        CVector3 r1( R1.mat[0*3+i], R1.mat[1*3+i], R1.mat[2*3+i] );
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
    CRigidBody rb(1.0, CVector3(0,1,0).stlvec() );
    aRigidBody.push_back(rb);
  }
  { // making plane 1
    CRigidBody rb(0.1, CVector3(-1, 0.5, -1).stlvec());
    rb.addCP(CVector3(-1.1,0,-1.1).stlvec());
    aRigidBody.push_back(rb);
  }
  { // making plane 2
    CRigidBody rb(0.1, CVector3(-1,0.5,+1).stlvec() );
    rb.addCP( CVector3(-1.1,0,+1.1).stlvec() );
    aRigidBody.push_back(rb);
  }
  { // making plane 3
    CRigidBody rb(0.1, CVector3(+1,0.5,-1).stlvec() );
    rb.addCP( CVector3(+1.1,0,-1.1).stlvec() );
    aRigidBody.push_back(rb);
  }
  { // making plane 4
    CRigidBody rb(0.1, CVector3(+1,0.5,+1).stlvec() );
    rb.addCP( CVector3(+1.1,0,+1.1).stlvec() );
    aRigidBody.push_back(rb);
  }
  
  
  { // joint 0
    CJoint jt(0,1, CVector3(-1,+1,-1).stlvec() );
    aJoint.push_back(jt);
  }
  { // joint 1
    CJoint jt(0,2, CVector3(-1,+1,+1).stlvec() );
    aJoint.push_back(jt);
  }
  { // joint 2
    CJoint jt(0,3, CVector3(+1,+1,-1).stlvec() );
    aJoint.push_back(jt);
  }
  { // joint 3
    CJoint jt(0,4, CVector3(+1,+1,+1).stlvec() );
    aJoint.push_back(jt);
  }
}


static void myGlVertex3d(const CVector3& v){
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



