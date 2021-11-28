/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/geo_quad.h"

#include <cmath>

#include "delfem2/geo_edge.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

template<typename T>
DFM2_INLINE delfem2::CVec3<T> delfem2::Nearst_Origin3_Quad3(
  double &s0,
  double &s1,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  const CVec3<T> &q2,
  const CVec3<T> &q3) {
  double dist_min = -1;
  CVec3<T> q_min;
  for (int ip = 0; ip < 5; ++ip) {
    double t0 = 0, t1 = 0;
    if (ip == 0) {
      t0 = 0.0;
      t1 = 0.0;
    } else if (ip == 1) {
      t0 = 1.0;
      t1 = 0.0;
    } else if (ip == 2) {
      t0 = 1.0;
      t1 = 1.0;
    } else if (ip == 3) {
      t0 = 0.0;
      t1 = 1.0;
    } else if (ip == 4) {
      t0 = 0.5;
      t1 = 0.5;
    }
    CVec3<T> q;
    for (int itr = 0; itr < 4; ++itr) {
      CVec3<T> pq = (1 - t0) * (1 - t1) * q0 + t0 * (1 - t1) * q1 + t0 * t1 * q2 + (1 - t0) * t1 * q3;
      CVec3<T> dqt0 = -(1 - t1) * q0 + (1 - t1) * q1 + t1 * q2 - t1 * q3;
      CVec3<T> dqt1 = -(1 - t0) * q0 - t0 * q1 + t0 * q2 + (1 - t0) * q3;
      CVec3<T> ddqt0t1 = q0 - q1 + q2 - q3;
      double f0 = -dqt0.dot(pq);
      double f1 = -dqt1.dot(pq);
      double A00 = dqt0.dot(dqt0);
      double A11 = dqt1.dot(dqt1);
      double A01 = dqt1.dot(dqt0) + ddqt0t1.dot(pq);
      double det = A00 * A11 - A01 * A01;
      double detinv = 1.0 / det;
      double B00 = +A11 * detinv;
      double B11 = +A00 * detinv;
      double B01 = -A01 * detinv;
      double d0 = B00 * f0 + B01 * f1;
      double d1 = B01 * f0 + B11 * f1;
      t0 += d0;
      t1 += d1;
    }
    double tol = 1.0e-4;
    if (t0 > -tol && t0 < 1.0 + tol && t1 > -tol && t1 < 1.0 + tol) {
      double d0 = q.norm();
      if (dist_min < 0 || d0 < dist_min) {
        dist_min = d0;
        s0 = t0;
        s1 = t1;
        q_min = q;
      }
    }
  }
  if (dist_min > 0) { return q_min; }
  //
  const CVec3<T> q01 = Nearest_Origin3_Edge3(q0, q1);
  const double d01 = q01.norm();
  if (dist_min < 0 || d01 < dist_min) {
    dist_min = d01;
    s0 = Distance(q01, q0) / Distance(q0, q1);
    s1 = 0.0;
    q_min = q01;
  }
  //
  CVec3<T> q12 = Nearest_Origin3_Edge3(q1, q2);
  const double d12 = q12.norm();
  if (dist_min < 0 || d12 < dist_min) {
    dist_min = d12;
    s0 = 1.0;
    s1 = Distance(q12, q1) / Distance(q1, q2);
    q_min = q12;
  }
  //
  CVec3<T> q23 = Nearest_Origin3_Edge3(q2, q3);
  const double d23 = q23.norm();
  if (dist_min < 0 || d23 < dist_min) {
    dist_min = d23;
    s0 = Distance(q23, q3) / Distance(q2, q3);
    s1 = 1.0;
    q_min = q23;
  }
  //
  CVec3<T> q30 = Nearest_Origin3_Edge3(q3, q0);
  const double d30 = q30.norm();
  if (dist_min < 0 || d30 < dist_min) {
    dist_min = d30;
    s0 = 0.0;
    s1 = Distance(q30, q0) / Distance(q3, q0);
    q_min = q30;
  }
  return q_min;
}
#ifdef DFM2_STATIC_LIBRARY
template DFM2_INLINE delfem2::CVec3d delfem2::Nearst_Origin3_Quad3(
  double &s0,
  double &s1,
  const CVec3d &q0,
  const CVec3d &q1,
  const CVec3d &q2,
  const CVec3d &q3);
#endif


// ==================================================


template<typename T>
void delfem2::iteration_intersection_Line_Quad
  (double &t0, double &t1,
   const CVec3<T> &src,
   const CVec3<T> &u,
   const CVec3<T> &v,
   const CVec3<T> &q0,
   const CVec3<T> &q1,
   const CVec3<T> &q2,
   const CVec3<T> &q3) {
  CVec3<T> q = (1 - t0) * (1 - t1) * q0 + t0 * (1 - t1) * q1 + t0 * t1 * q2 + (1 - t0) * t1 * q3;
  CVec3<T> pq = q - src;
  CVec3<T> dqt0 = -(1 - t1) * q0 + (1 - t1) * q1 + t1 * q2 - t1 * q3;
  CVec3<T> dqt1 = -(1 - t0) * q0 - t0 * q1 + t0 * q2 + (1 - t0) * q3;
  CVec3<T> ddqt0t1 = q0 - q1 + q2 - q3;
  double f0 = -u.dot(pq);
  double f1 = -v.dot(pq);
  double A00 = u.dot(dqt0);
  double A01 = u.dot(dqt1);
  double A10 = v.dot(dqt0);
  double A11 = v.dot(dqt1);
  double det = A00 * A11 - A01 * A10;
  double detinv = 1.0 / det;
  double B00 = +A11 * detinv;
  double B01 = -A01 * detinv;
  double B10 = -A10 * detinv;
  double B11 = +A00 * detinv;
  double d0 = B00 * f0 + B01 * f1;
  double d1 = B10 * f0 + B11 * f1;
  t0 += d0;
  t1 += d1;
}

template<typename T>
bool delfem2::intersection_Point_Quad(
  CVec3<T> &psec, double &s0, double &s1,
  const CVec3<T> &src, const CVec3<T> &dir,
  const CVec3<T> &q0, const CVec3<T> &q1, const CVec3<T> &q2, const CVec3<T> &q3) {
  CVec3<T> u, v;
  GetVertical2Vector(dir, u, v);
  //
  double dist_min = -1;
  CVec3<T> q_min;
  for (int ip = 0; ip < 5; ++ip) {
    double t0 = 0, t1 = 0;
    if (ip == 0) {
      t0 = 0.0;
      t1 = 0.0;
    } else if (ip == 1) {
      t0 = 1.0;
      t1 = 0.0;
    } else if (ip == 2) {
      t0 = 1.0;
      t1 = 1.0;
    } else if (ip == 3) {
      t0 = 0.0;
      t1 = 1.0;
    } else if (ip == 4) {
      t0 = 0.5;
      t1 = 0.5;
    }
    for (int itr = 0; itr < 4; ++itr) {
      iteration_intersection_Line_Quad(t0, t1,
                                       src, u, v,
                                       q0, q1, q2, q3);
    }
    CVec3<T> q = (1 - t0) * (1 - t1) * q0 + t0 * (1 - t1) * q1 + t0 * t1 * q2 + (1 - t0) * t1 * q3;
    double tol = 1.0e-4;
    //    std::cout << t0 << " " << t1 << std::endl;
    if (t0 > -tol && t0 < 1.0 + tol && t1 > -tol && t1 < 1.0 + tol) {
      double d0 = (q - src).norm();
      if (dist_min < 0 || d0 < dist_min) {
        dist_min = d0;
        s0 = t0;
        s1 = t1;
        q_min = q;
      }
    }
  }
  //  std::cout << dist_min << std::endl;
  if (dist_min > 0) {
    psec = q_min;
    return true;
  }
  return false;
}
#ifdef DFM2_STATIC_LIBRARY
template bool delfem2::intersection_Point_Quad(
  CVec3d &psec, 
  double &s0, double &s1,
  const CVec3d &src, 
  const CVec3d &dir,
  const CVec3d &q0, 
  const CVec3d &q1, 
  const CVec3d &q2, 
  const CVec3d &q3);
#endif


/*
CVector3 nearest_Origin_Tet
(double& r0, double& r1, double& r2,
 const CVector3& q0,
 const CVector3& q1,
 const CVector3& q2,
 const CVector3& q3)
{
  CVector3 p_min = q0;
  {
    bool res = barycentricCoord_Origin_Tet(r0, r1, r2, q0, q1, q2, q3);
    p_min = r0*q0 + r1*q1 + r2*q2 + (1-r0-r1-r2)*q3;
    if( r0>0 && r1>0 && r2>0 && (1-r0-r1-r2)>0 ){ return p_min; }
  }
  ////////////////////////
  double r3;
  { // face123
    r0 = 0;
    p_min = nearest_Orgin_PlaneTri(r1,r2,r3, q1,q2,q3);
    if( r1>0 && r2>0 && r3>0 ){ return p_min; }
  }
  { // face230
    r1 = 0;
    p_min = nearest_Orgin_PlaneTri(r2,r3,r0, q2,q3,q0);
    if( r2>0 && r3>0 && r0>0 ){ return p_min; }
  }
  { // face301
    r2 = 0;
    p_min = nearest_Orgin_PlaneTri(r3,r0,r1, q3,q0,q1);
    if( r3>0 && r0>0 && r1>0 ){ return p_min; }
  }
  { // face012
    r3 = 0;
    p_min = nearest_Orgin_PlaneTri(r0,r1,r2, q0,q1,q2);
    if( r0>0 && r1>0 && r2>0 ){ return p_min; }
  }
  ////////////////////////
  double d_min = q0.Length();
  double s0,s1,s2,s3;
  { // edge01
    CVector3 p01 = nearest_Origin_LineSeg(s0,s1,q0,q1);
    double d01 = p01.Length();
    if( d01<d_min ){ d_min=d01; p_min=p01; r0=s0; r1=s1; r2=0; r3=0; }
  }
  { // edge02
    CVector3 p02 = nearest_Origin_LineSeg(s0,s2,q0,q2);
    double d02 = p02.Length();
    if( d02<d_min ){ d_min=d02; p_min=p02; r0=s0; r1=0; r2=s2; r3=0; }
  }
  { // edge03
    CVector3 p03 = nearest_Origin_LineSeg(s0,s3,q0,q3);
    double d03 = p03.Length();
    if( d03<d_min ){ d_min=d03; p_min=p03; r0=s0; r1=0; r2=0; r3=s3; }
  }
  { // edge12
    CVector3 p12 = nearest_Origin_LineSeg(s1,s2,q1,q2);
    double d12 = p12.Length();
    if( d12<d_min ){ d_min=d12; p_min=p12; r0=0; r1=s1; r2=s2; r3=0; }
  }
  { // edge13
    CVector3 p13 = nearest_Origin_LineSeg(s1,s3,q1,q3);
    double d13 = p13.Length();
    if( d13<d_min ){ d_min=d13; p_min=p13; r0=0; r1=s1; r2=0; r3=s3; }
  }
  { // edge23
    CVector3 p23 = nearest_Origin_LineSeg(s2,s3,q2,q3);
    double d23 = p23.Length();
    if( d23<d_min ){ d_min=d23; p_min=p23; r0=0; r1=0; r2=s2; r3=s3; }
  }
  return p_min;
}
 */