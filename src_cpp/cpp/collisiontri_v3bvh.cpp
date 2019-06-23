#include "delfem2/collisiontri_v3bvh.h"

// distance VF
double DistanceFaceVertex
(const CVector3& p0, const CVector3& p1, const CVector3& p2,
 const CVector3& p3,
 double& w0, double& w1)
{
  CVector3 v20 =p0-p2;
  CVector3 v21 =p1-p2;
  double t0 = Dot(v20,v20);
  double t1 = Dot(v21,v21);
  double t2 = Dot(v20,v21);
  double t3 = Dot(v20,p3-p2);
  double t4 = Dot(v21,p3-p2);
  double det = t0*t1-t2*t2;
  double invdet = 1.0/det;
  w0 = (+t1*t3-t2*t4)*invdet;
  w1 = (-t2*t3+t0*t4)*invdet;
  const double w2 = 1-w0-w1;
  CVector3 pw = w0*p0 + w1*p1 + w2*p2;
  return (pw-p3).Length();
}

//　distance EE
double DistanceEdgeEdge
(const CVector3& p0, const CVector3& p1,
 const CVector3& q0, const CVector3& q1,
 double& ratio_p, double& ratio_q)
{
  const CVector3& vp =p1-p0;
  const CVector3& vq =q1-q0;
  if( Cross(vp,vq).Length() < 1.0e-10 ){ // handling parallel edge
    CVector3 pq0 = p0-q0;
    CVector3 nvp = vp; nvp.SetNormalizedVector();
    CVector3 vert = pq0 - Dot(pq0,nvp)*nvp;
    double dist = vert.Length();
    double lp0 = Dot(p0,nvp);
    double lp1 = Dot(p1,nvp);
    double lq0 = Dot(q0,nvp);
    double lq1 = Dot(q1,nvp);
    double p_min  = ( lp0 < lp1 ) ? lp0 : lp1;
    double p_max  = ( lp0 > lp1 ) ? lp0 : lp1;
    double q_min  = ( lq0 < lq1 ) ? lq0 : lq1;
    double q_max  = ( lq0 > lq1 ) ? lq0 : lq1;
    double lm;
    if(      p_max < q_min ){ lm = (p_max+q_min)*0.5; }
    else if( q_max < p_min ){ lm = (q_max+p_min)*0.5; }
    else if( p_max < q_max ){ lm = (p_max+q_min)*0.5; }
    else{                     lm = (q_max+p_min)*0.5; }
    ratio_p = (lm-lp0)/(lp1-lp0);
    ratio_q = (lm-lq0)/(lq1-lq0);
    return dist;
  }
  double t0 = Dot(vp,vp);
  double t1 = Dot(vq,vq);
  double t2 = Dot(vp,vq);
  double t3 = Dot(vp,q0-p0);
  double t4 = Dot(vq,q0-p0);
  double det = t0*t1-t2*t2;
  double invdet = 1.0/det;
  ratio_p = (+t1*t3-t2*t4)*invdet;
  ratio_q = (+t2*t3-t0*t4)*invdet;
  CVector3 pc = p0 + ratio_p*vp;
  CVector3 qc = q0 + ratio_q*vq;
  return (pc-qc).Length();
}

// EEの距離が所定の距離以下にあるかどうか
bool IsContact_EE_Proximity
(int ino0,        int ino1,        int jno0,        int jno1,
 const CVector3& p0, const CVector3& p1, const CVector3& q0, const CVector3& q1,
 const double delta)
{
  if( ino0 == jno0 || ino0 == jno1 || ino1 == jno0 || ino1 == jno1 ) return false;
  if( q0.x+delta < p0.x && q0.x+delta < p1.x && q1.x+delta < p0.x && q1.x+delta < p1.x ) return false;
  if( q0.x-delta > p0.x && q0.x-delta > p1.x && q1.x-delta > p0.x && q1.x-delta > p1.x ) return false;
  if( q0.y+delta < p0.y && q0.y+delta < p1.y && q1.y+delta < p0.y && q1.y+delta < p1.y ) return false;
  if( q0.y-delta > p0.y && q0.y-delta > p1.y && q1.y-delta > p0.y && q1.y-delta > p1.y ) return false;
  if( q0.z+delta < p0.z && q0.z+delta < p1.z && q1.z+delta < p0.z && q1.z+delta < p1.z ) return false;
  if( q0.z-delta > p0.z && q0.z-delta > p1.z && q1.z-delta > p0.z && q1.z-delta > p1.z ) return false;
  double ratio_p, ratio_q;
  double dist = DistanceEdgeEdge(p0, p1, q0, q1, ratio_p, ratio_q);
  if( dist > delta ) return false;
  if( ratio_p < 0 ) return false;
  if( ratio_p > 1 ) return false;
  if( ratio_q < 0 ) return false;
  if( ratio_q > 1 ) return false;
  const CVector3& pm = (1-ratio_p)*p0 + ratio_p*p1;
  const CVector3& qm = (1-ratio_q)*q0 + ratio_q*q1;
  if( (pm-qm).Length() > delta ) return false;
  return true;
}

// 三次関数を評価する関数
inline double EvaluateCubic
(double r2, // 入力の値
 double k0, double k1, double k2, double k3) // 三次関数の係数
{
  return k0 + k1*r2 + k2*r2*r2 + k3*r2*r2*r2;
}

// 二分法における三次関数の根を探す範囲を狭める関数
static void BisectRangeCubicRoot
(int& icnt,                  // (in out)何回幅を狭めたかというカウンタ
 double& r0, double& r1, // (in,out)入力の範囲から階を探して、狭められた後の範囲を返す
 double v0, double v1, // (in)入力の範囲の両端における値
 double k0, double k1, double k2, double k3) // 三次関数の係数
{
  icnt--;
  if( icnt <= 0 ) return;
  double r2 = 0.5*(r0+r1); // r2はr0とr1の中点
  double v2 = EvaluateCubic(r2, k0,k1,k2,k3); // v2はr2における値
  if( v0*v2 < 0 ){ r1 = r2; } // r0とr2の間で符号が変化する
  else{            r0 = r2; } // r1とr2の間で符号が変化する
  BisectRangeCubicRoot(icnt,r0,r1,v0,v2,k0,k1,k2,k3); // r0とr1の間でさらに範囲を狭める
}

// 三次関数の根を探す関数
static double FindRootCubic
(double r0, double r1,
 double v0, double v1,
 double k0, double k1, double k2, double k3)
{
  int icnt=15; // １５回範囲を狭める
  BisectRangeCubicRoot(icnt, r0,r1, v0,v1, k0,k1,k2,k3);
  return 0.5*(r0+r1);
}

// ４つの点が同一平面上にならぶような補間係数を探す
double FindCoplanerInterp
(const CVector3& s0, const CVector3& s1, const CVector3& s2, const CVector3& s3,
 const CVector3& e0, const CVector3& e1, const CVector3& e2, const CVector3& e3)
{
  const CVector3 x1 = s1-s0;
  const CVector3 x2 = s2-s0;
  const CVector3 x3 = s3-s0;
  const CVector3 v1 = e1-e0-x1;
  const CVector3 v2 = e2-e0-x2;
  const CVector3 v3 = e3-e0-x3;
  // 三次関数の係数の計算
  const double k0 = ScalarTripleProduct(x3,x1,x2);
  const double k1 = ScalarTripleProduct(v3,x1,x2)+ScalarTripleProduct(x3,v1,x2)+ScalarTripleProduct(x3,x1,v2);
  const double k2 = ScalarTripleProduct(v3,v1,x2)+ScalarTripleProduct(v3,x1,v2)+ScalarTripleProduct(x3,v1,v2);
  const double k3 = ScalarTripleProduct(v3,v1,v2);
  double r0=-0.0;
  double r1=+1.0;
  const double f0 = EvaluateCubic(r0,k0,k1,k2,k3);
  const double f1 = EvaluateCubic(r1,k0,k1,k2,k3);
  double det = k2*k2-3*k1*k3;
  if( fabs(k3) < 1.0e-10 && fabs(k2) > 1.0e-10 ){ // quadric function、二次関数
    double r2 = -k1/(2*k2); // 極値をとるr
    const double f2 = EvaluateCubic(r2, k0,k1,k2,k3);
    if( r2 > 0 && r2 < 1 ){
      if(      f0*f2 < 0 ){
        return FindRootCubic(r0,r2, f0,f2, k0,k1,k2,k3);
        
      }
      else if( f2*f1 < 0 ){
        return FindRootCubic(r2,r1, f2,f1, k0,k1,k2,k3);
      }
    }
  }
  if( det > 0 && fabs(k3) > 1.0e-10 ) // cubic function with two extream value、三次関数で極値がある場合
  {
    double r3 = (-k2-sqrt(det))/(3*k3); // 極値をとる小さい方のr
    const double f3 = EvaluateCubic(r3, k0,k1,k2,k3);
    if( r3 > 0 && r3 < 1 ){
      if(      f0*f3 < 0 ){
        return FindRootCubic(r0,r3, f0,f3, k0,k1,k2,k3);
      }
      else if( f3*f1 < 0 ){
        return FindRootCubic(r3,r1, f3,f1, k0,k1,k2,k3);
      }
    }
    double r4 = (-k2+sqrt(det))/(3*k3); // 極値をとる大きい方のr
    const double f4 = EvaluateCubic(r4, k0,k1,k2,k3);
    if( r3 > 0 && r3 < 1 && r4 > 0 && r4 < 1 ){
      if( f3*f4 < 0 ){
        return FindRootCubic(r3,r4, f3,f4, k0,k1,k2,k3);
      }
    }
    if( r4 > 0 && r4 < 1 ){
      if(      f0*f4 < 0 ){
        return FindRootCubic(r0,r4, f0,f4, k0,k1,k2,k3);
      }
      else if( f4*f1 < 0 ){
        return FindRootCubic(r4,r1, f4,f1, k0,k1,k2,k3);
      }
    }
  }
  // monotonus function、0と１の間で短調増加関数
  if( f0*f1 > 0 ){ return -1; } // 根がない場合
  return FindRootCubic(r0,r1, f0,f1, k0,k1,k2,k3);
}

// CCDのFVで接触する要素を検出
bool IsContact_FV_CCD2
(int ino0,        int ino1,        int ino2,        int ino3,
 const CVector3& p0, const CVector3& p1, const CVector3& p2, const CVector3& p3,
 const CVector3& q0, const CVector3& q1, const CVector3& q2, const CVector3& q3)
{
  { // CSAT
    CVector3 n = Cross(p1-p0,p2-p0);
    double t0 = Dot(p0-p3,n);
    double t1 = Dot(q0-q3,n);
    double t2 = Dot(q1-q3,n);
    double t3 = Dot(q2-q3,n);
    if( t0*t1 > 0 && t0*t2 > 0 && t0*t3 > 0 ){ return false; }
  }
  double r0,r1;
  double dist = DistanceFaceVertex(p0, p1, p2, p3, r0,r1);
  {
    double vn0 = (p0-q0).Length();
    double vn1 = (p1-q1).Length();
    double vn2 = (p2-q2).Length();
    double vn3 = (p3-q3).Length();
    double vnt = ( vn0 > vn1 ) ? vn0 : vn1;
    vnt = ( vn2 > vnt ) ? vn2 : vnt;
    double max_app = (vnt+vn3);
    ////
    const double r2 = 1-r0-r1;
    if( dist > max_app ) return false;
    if( r0 < 0 || r0 > 1 || r1 < 0 || r1 > 1 || r2 < 0 || r2 > 1 ){
      double dist01 = (nearest_LineSeg_Point(p3, p0, p1)-p3).Length();
      double dist12 = (nearest_LineSeg_Point(p3, p1, p2)-p3).Length();
      double dist20 = (nearest_LineSeg_Point(p3, p2, p0)-p3).Length();
      if( dist01 > max_app && dist12 > max_app && dist20 > max_app ){ return false; }
    }
  }
  double t = FindCoplanerInterp(p0,p1,p2,p3, q0,q1,q2,q3);
  if( t < 0 || t > 1 ) return false;
  CVector3 p0m = (1-t)*p0 + t*q0;
  CVector3 p1m = (1-t)*p1 + t*q1;
  CVector3 p2m = (1-t)*p2 + t*q2;
  CVector3 p3m = (1-t)*p3 + t*q3;
  double w0, w1;
  DistanceFaceVertex(p0m, p1m, p2m, p3m, w0,w1);
  double w2 = 1-w0-w1;
  if( w0 < 0 || w0 > 1 ) return false;
  if( w1 < 0 || w1 > 1 ) return false;
  if( w2 < 0 || w2 > 1 ) return false;
  return true;
}


bool isIntersectTriPair
(CVector3& P0, CVector3& P1,
 int itri, int jtri,
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aXYZ)
{
  const int i0 = aTri[itri*3+0];
  const int i1 = aTri[itri*3+1];
  const int i2 = aTri[itri*3+2];
  const int j0 = aTri[jtri*3+0];
  const int j1 = aTri[jtri*3+1];
  const int j2 = aTri[jtri*3+2];
  if( i0 == j0 || i0 == j1 || i0 == j2 ) return false;
  if( i1 == j0 || i1 == j1 || i1 == j2 ) return false;
  if( i2 == j0 || i2 == j1 || i2 == j2 ) return false;
  const CVector3 p0(aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]);
  const CVector3 p1(aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]);
  const CVector3 p2(aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2]);
  const CVector3 q0(aXYZ[j0*3+0], aXYZ[j0*3+1], aXYZ[j0*3+2]);
  const CVector3 q1(aXYZ[j1*3+0], aXYZ[j1*3+1], aXYZ[j1*3+2]);
  const CVector3 q2(aXYZ[j2*3+0], aXYZ[j2*3+1], aXYZ[j2*3+2]);
  const CVector3 np = Normal(p0,p1,p2);
  const CVector3 nq = Normal(q0,q1,q2);
  double dp0 = (p0-q0)*nq;
  double dp1 = (p1-q0)*nq;
  double dp2 = (p2-q0)*nq;
  double dq0 = (q0-p0)*np;
  double dq1 = (q1-p0)*np;
  double dq2 = (q2-p0)*np;
  if( ((dp0>0) == (dp1>0)) && ((dp1>0) == (dp2>0)) ) return false;
  if( ((dq0>0) == (dq1>0)) && ((dq1>0) == (dq2>0)) ) return false;
  const CVector3 p01 = (1.0/(dp0-dp1))*(dp0*p1-dp1*p0);
  const CVector3 p12 = (1.0/(dp1-dp2))*(dp1*p2-dp2*p1);
  const CVector3 p20 = (1.0/(dp2-dp0))*(dp2*p0-dp0*p2);
  const CVector3 q01 = (1.0/(dq0-dq1))*(dq0*q1-dq1*q0);
  const CVector3 q12 = (1.0/(dq1-dq2))*(dq1*q2-dq2*q1);
  const CVector3 q20 = (1.0/(dq2-dq0))*(dq2*q0-dq0*q2);
  const CVector3 vz = Cross(np,nq);
  CVector3 ps,pe;
  if(      dp0*dp1>0 ){ ps=p20; pe=p12; }
  else if( dp1*dp2>0 ){ ps=p01; pe=p20; }
  else{                 ps=p12; pe=p01; }
  if( ps*vz>pe*vz ){ CVector3 pt=ps; ps=pe; pe=pt; }
  double zps = ps*vz;
  double zpe = pe*vz;
  assert( zps<=zpe );
  ////
  CVector3 qs,qe;
  if(      dq0*dq1>0 ){ qs=q20; qe=q12; }
  else if( dq1*dq2>0 ){ qs=q01; qe=q20; }
  else{                 qs=q12; qe=q01; }
  if( qs*vz>qe*vz ){ CVector3 qt=qs; qs=qe; qe=qt; }
  double zqs = qs*vz;
  double zqe = qe*vz;
  assert( zqs<=zqe );
  ////
  if( zps>zqe || zqs>zpe ) return false;
  CVector3 P[4];
  int icnt = 0;
  if( zps>zqs && zps<zqe ){ P[icnt]=ps; icnt++; }
  if( zpe>zqs && zpe<zqe ){ P[icnt]=pe; icnt++; }
  if( zqs>zps && zqs<zpe ){ P[icnt]=qs; icnt++; }
  if( zqe>zps && zqe<zpe ){ P[icnt]=qe; icnt++; }
  if( icnt != 2 ) return false;
  P0 = P[0];
  P1 = P[1];
  return true;
}
