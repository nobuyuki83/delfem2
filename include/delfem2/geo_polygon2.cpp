/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>

#include "delfem2/geo_polygon2.h"

// ================================================

DFM2_INLINE void delfem2::makeSplineLoop(
    const std::vector<double>& aCV,
    std::vector<double>& aVecCurve)
{
  aVecCurve.resize(0);
  const int nCV = (int)aCV.size()/2;
  unsigned int ndiv = 5;
  for(int icv=0;icv<nCV;icv++){
    int icv0=icv;   if( icv0 >= nCV ){ icv0-=nCV; }
    int icv1=icv+1; if( icv1 >= nCV ){ icv1-=nCV; }
    int icv2=icv+2; if( icv2 >= nCV ){ icv2-=nCV; }
    const double p0[2] = { aCV[icv0*2+0], aCV[icv0*2+1] };
    const double p1[2] = { aCV[icv1*2+0], aCV[icv1*2+1] };
    const double p2[2] = { aCV[icv2*2+0], aCV[icv2*2+1] };
    for(unsigned int idiv=0;idiv<ndiv;idiv++){
      const double t = 1.0-(double)idiv/ndiv;
      const double w[3] = {0.5*t*t, -t*t + t + 0.5, 0.5*(1-t)*(1-t) };
      const double px = w[0]*p0[0] + w[1]*p1[0] + w[2]*p2[0];
      const double py = w[0]*p0[1] + w[1]*p1[1] + w[2]*p2[1];
      aVecCurve.push_back(px);
      aVecCurve.push_back(py);
    }
  }
}


template <typename T>
double delfem2::Length_Polygon(
    const std::vector<CVec2<T>>& aP)
{
  if( aP.size() < 2 ){ return 0; }
  const size_t np = aP.size();
  double len = 0;
  for(unsigned int ip0=0;ip0<np;ip0++){
    const unsigned int ip1 = (ip0+1)%np;
    len += (aP[ip0]-aP[ip1]).norm();
  }
  return len;
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::Length_Polygon(const std::vector<CVec2d>& aP);
#endif

// ---------------------------------------

template <typename T>
double delfem2::Area_Polygon(
    const std::vector<CVec2<T>>& aP)
{
  CVec2<T> vtmp(0,0);
  const int ne = aP.size();
  double area_loop = 0.0;
  for(int ie=0;ie<ne;ie++){
    area_loop += Area_Tri(vtmp, aP[(ie+0)%ne], aP[(ie+1)%ne]);
  }
  return area_loop;
}


template <typename T>
void makeRandomLoop(
    unsigned int nCV,
    std::vector<double>& aCV)
{
  aCV.clear();
  for(unsigned int icv=0;icv<nCV;icv++){
    /*
     {
     aCV.push_back((double)rand()/(RAND_MAX+1.0));
     aCV.push_back((double)rand()/(RAND_MAX+1.0));
     }
     */
    { // polar coordinate random position
      double tht = icv*3.1415*2.0/nCV;
      double r = (double)rand()/(RAND_MAX+1.0);
      double px = r*sin(tht);
      double py = r*cos(tht);
      aCV.push_back(px);
      aCV.push_back(py);
    }
  }
}

// ------------------------------------------------------------

template <typename T>
std::vector<delfem2::CVec2<T>> delfem2::Polygon_Resample_Polygon(
    const std::vector<CVec2<T>>& stroke0,
    double l)
{
  std::vector<CVec2<T>> stroke;
  if( stroke0.empty() ) return stroke;
  stroke.push_back( stroke0[0] );
  int jcur = 0;
  double rcur = 0;
  double lcur = l;
  for(;;){
    if( jcur > (int)stroke0.size()-1 ) break;
    int ip0 = jcur;
    int ip1 = (jcur+1)%stroke0.size();
    double lenj  = (stroke0[ip1]-stroke0[ip0]).Length();
    double lenjr = lenj*(1.0-rcur);
    if( lenjr > lcur ){ // put point in this segment
      rcur += lcur/lenj;
      stroke.push_back( (1-rcur)*stroke0[ip0] + rcur*stroke0[ip1] );
      lcur = l;
    }
    else{ // next segment
      lcur -= lenjr;
      rcur =0;
      jcur++;
    }
  }
  //  stroke.push_back( stroke0.back() );
  return stroke;
}

template <typename T>
void delfem2::CgArea_Polygon(
    CVec2<T>& cg,
    T& area,
    const std::vector<CVec2<T>>& aVec2D)
{
  area = 0;
  const size_t nseg = aVec2D.size();
  cg = CVec2<T>(0.0, 0.0);
  for(unsigned int iseg=0;iseg<nseg;iseg++){
    unsigned int ip0 = (iseg+0)%nseg;
    unsigned int ip1 = (iseg+1)%nseg;
    const T x0 = aVec2D[ip0].p[0];
    const T y0 = aVec2D[ip0].p[1];
    const T x1 = aVec2D[ip1].p[0];
    const T y1 = aVec2D[ip1].p[1];
    const T ai = x0*y1 - x1*y0;
    area += ai;
    cg.p[0] += ai*(x0+x1)/3;
    cg.p[1] += ai*(y0+y1)/3;
  }
  cg.p[0] /= area;
  cg.p[1] /= area;
  area *= 0.5;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CgArea_Polygon(
    CVec2d& cg,
    double& area,
    const std::vector<CVec2d>& aVec2D);
template void delfem2::CgArea_Polygon(
    CVec2f& cg,
    float& area,
    const std::vector<CVec2f>& aVec2D);
#endif

// https://en.wikipedia.org/wiki/List_of_moments_of_inertia
template <typename T>
T delfem2::RotationalMomentPolar_Polygon2(
    const std::vector<CVec2<T>>& aVec2,
    const CVec2<T>& pivot)
{
  const size_t ne = aVec2.size();
  T sum_I = 0.0;
  for(unsigned int ie=0;ie<ne;++ie){
    const unsigned int ip0 = ie;
    const unsigned int ip1 = (ie+1)%ne;
    const CVec2<T> p0 = aVec2[ip0] - pivot;
    const CVec2<T> p1 = aVec2[ip1] - pivot;
    T a0 = (p0^p1)*static_cast<T>(0.5);
    sum_I += a0*(p0.dot(p0)+p0.dot(p1)+p1.dot(p1));
  }
  return sum_I*static_cast<T>(1.0/6.0);
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::RotationalMomentPolar_Polygon2(
    const std::vector<CVec2d>& aVec2,
    const CVec2d& pivot);
template float delfem2::RotationalMomentPolar_Polygon2(
    const std::vector<CVec2f>& aVec2,
    const CVec2f& pivot);
#endif


template <typename T>
void delfem2::SecondMomentOfArea_Polygon(
    CVec2<T>& cg,  double& area,
    CVec2<T>& pa1, double& I1,
    CVec2<T>& pa2, double& I2,
    const std::vector<CVec2<T>>& aVec2D)
{
  area = 0;
  const size_t nseg = aVec2D.size();
  cg = CVec2<T>(0.0, 0.0);
  for(unsigned int iseg=0;iseg<nseg;iseg++){
    unsigned int ip0 = (iseg+0)%nseg;
    unsigned int ip1 = (iseg+1)%nseg;
    double x0 = aVec2D[ip0].p[0];
    double y0 = aVec2D[ip0].p[1];
    double x1 = aVec2D[ip1].p[0];
    double y1 = aVec2D[ip1].p[1];
    double ai = x0*y1 - x1*y0;
    area += ai;
    cg.p[0] += ai*(x0+x1)/3.0;
    cg.p[1] += ai*(y0+y1)/3.0;
  }
  cg.p[0] /= area;
  cg.p[1] /= area;
  area *= 0.5;
  // ----------------
  double Ix=0, Iy=0, Ixy=0;
  for(unsigned int iseg=0;iseg<nseg;iseg++){
    unsigned int ip0 = (iseg+0)%nseg;
    unsigned int ip1 = (iseg+1)%nseg;
    double x0 = aVec2D[ip0].p[0]-cg.p[0];
    double y0 = aVec2D[ip0].p[1]-cg.p[1];
    double x1 = aVec2D[ip1].p[0]-cg.p[0];
    double y1 = aVec2D[ip1].p[1]-cg.p[1];
    double ai = x0*y1 - x1*y0;
    Ix  += ai*(y0*y0 + y0*y1 + y1*y1)/12.0;
    Iy  += ai*(x0*x0 + x0*x1 + x1*x1)/12.0;
    Ixy += ai*(x0*y0 + 0.5*x0*y1 + 0.5*x1*y0 + x1*y1)/12.0;
  }
  if( fabs(Ix-Iy)+fabs(Ixy) < 1.0e-20 ){
    pa1 = CVec2<T>(1,0);
    pa2 = CVec2<T>(0,1);
    return;
  }
  double phi = 0.5*atan2(-2*Ixy,Ix-Iy);
  pa1 = CVec2<T>(+cos(phi), +sin(phi));
  pa2 = CVec2<T>(-sin(phi), +cos(phi));
  I1 = 0.5*(Ix+Iy)+0.5*sqrt( (Ix-Iy)*(Ix-Iy) + 4*Ixy*Ixy );
  I2 = 0.5*(Ix+Iy)-0.5*sqrt( (Ix-Iy)*(Ix-Iy) + 4*Ixy*Ixy );
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::SecondMomentOfArea_Polygon(
    CVec2d& cg,  double& area,
    CVec2d& pa1, double& I1,
    CVec2d& pa2, double& I2,
    const std::vector<CVec2d>& aVec2D);
#endif

// ----------------------------------------

template <typename T>
void delfem2::JArray_FromVecVec_XY(
    std::vector<int>& aIndXYs,
    std::vector<int>& loopIP0,
    std::vector<CVec2<T>>& aXY,
    const std::vector< std::vector<double> >& aVecAry0)
{
  const int nloop = (int)aVecAry0.size();
  aIndXYs.resize(nloop+1);
  aIndXYs[0] = 0;
  int npo_sum = 0;
  for(int iloop=0;iloop<(int)nloop;iloop++){
    const int npo = (int)aVecAry0[iloop].size()/2;
    aIndXYs[iloop+1] = aIndXYs[iloop]+npo;
    npo_sum += npo;
  }
  aXY.resize(npo_sum);
  npo_sum = 0;
  for(int iloop=0;iloop<(int)nloop;iloop++){
    const int nxys = (int)aVecAry0[iloop].size()/2;
    for(int ixys=0;ixys<nxys;ixys++){
      aXY[npo_sum].p[0] = aVecAry0[iloop][ixys*2+0];
      aXY[npo_sum].p[1] = aVecAry0[iloop][ixys*2+1];
      npo_sum++;
    }
  }
  loopIP0.resize(aXY.size());
  for(unsigned int ip=0;ip<aXY.size();++ip){
    loopIP0[ip] = ip;
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::JArray_FromVecVec_XY(
    std::vector<int>& aIndXYs, std::vector<int>& loopIP0, std::vector<CVec2d>& aXY,
    const std::vector< std::vector<double> >& aVecAry0);
#endif

// --------------------------

template <typename T>
void delfem2::ResamplingLoop(
    std::vector<int>& loopIP1_ind,
    std::vector<int>& loopIP1,
    std::vector<CVec2<T>>& aVec2,
    double max_edge_length)
{
  assert( aVec2.size() == loopIP1.size() );
  const std::vector<int> loopIP0_ind = loopIP1_ind;
  const std::vector<int> loopIP0 = loopIP1;
  assert( loopIP0.size() >= 2 );
  const std::size_t nloop = loopIP0_ind.size()-1;
  std::vector< std::vector<int> > aPoInEd(loopIP0.size());
  {
    for(unsigned int iloop=0;iloop<nloop;++iloop){
	  assert(loopIP0_ind[iloop + 1] > loopIP0_ind[iloop]);
      const unsigned int np = loopIP0_ind[iloop+1]-loopIP0_ind[iloop];
      for(unsigned int ip=0;ip<np;ip++){
        const int iipo0 = loopIP0_ind[iloop]+(ip+0)%np; assert( iipo0>=0 && iipo0<(int)loopIP0.size() );
        const int iipo1 = loopIP0_ind[iloop]+(ip+1)%np; assert( iipo1>=0 && iipo1<(int)loopIP0.size() );
        const int ipo0 = loopIP0[iipo0]; assert(ipo0>=0&&ipo0<(int)aVec2.size());
        const int ipo1 = loopIP0[iipo1]; assert(ipo1>=0&&ipo1<(int)aVec2.size());
        const CVec2<T> po0 = aVec2[ipo0]; // never use reference here because aVec2 will resize afterward
        const CVec2<T> po1 = aVec2[ipo1]; // never use reference here because aVec2 will resize afterward
        const int nadd = (int)( Distance( po0, po1 ) / max_edge_length);
        if( nadd == 0 ) continue;
        for(int iadd=0;iadd<nadd;++iadd){
          double r2 = (double)(iadd+1)/(nadd+1);
          CVec2<T> v2 = (1-r2)*po0 + r2*po1;
          const unsigned int ipo2 = static_cast<unsigned int>(aVec2.size());
          aVec2.push_back(v2);
          assert( iipo0>=0 && iipo0<(int)aPoInEd.size() );
          aPoInEd[ iipo0 ].push_back(ipo2);
        }
      }
    }
  }
  ////
  loopIP1_ind.resize(nloop+1);
  loopIP1_ind[0] = 0;
  for(unsigned int iloop=0;iloop<nloop;++iloop){
    const int nbar0 = loopIP0_ind[iloop+1]-loopIP0_ind[iloop];
    int nbar1 = nbar0;
    for(int ibar=0;ibar<nbar0;ibar++){
      const int iip_loop = loopIP0_ind[iloop]+ibar;
      nbar1 += static_cast<int>(aPoInEd[iip_loop].size());
    }
    loopIP1_ind[iloop+1] = loopIP1_ind[iloop] + nbar1;
  }
  // adding new vertices on the outline
  loopIP1.resize(loopIP1_ind[nloop]);
  unsigned int ivtx0 = 0;
  for(unsigned int iloop=0;iloop<nloop;iloop++){
    for(int iip_loop=loopIP0_ind[iloop];iip_loop<loopIP0_ind[iloop+1];iip_loop++){
      const int ip_loop = loopIP0[iip_loop];
      loopIP1[ivtx0] = ip_loop;
      ivtx0++;
      for(std::size_t iadd=0;iadd<aPoInEd[ip_loop].size();iadd++){
        loopIP1[ivtx0] = aPoInEd[iip_loop][iadd];
        ivtx0++;
      }
    }
  }
  assert( loopIP1.size() == aVec2.size() );
  assert( loopIP1.size() == ivtx0 );
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::ResamplingLoop(
    std::vector<int>& loopIP1_ind, std::vector<int>& loopIP1,
    std::vector<CVec2d>& aVec2,
    double max_edge_length);
#endif

template <typename T>
bool delfem2::IsInclude_Loop(
    const double co[],
    const int ixy_stt,
    const int ixy_end,
    const std::vector<CVec2<T>>& aXY)
{
  int inum_cross = 0;
  for(int itr=0;itr<10;itr++){
    const double dir[2] = { cos((itr+1)*23.0), sin((itr+1)*23.0) }; // random direction
    const double codir[2] = { co[0]+dir[0], co[1]+dir[1] };
    bool is_fail = false;
    inum_cross = 0;
    for(int ixys=ixy_stt;ixys<ixy_end;ixys++){
      const int ipo0 = ixys;
      int ipo1 = ixys+1;
      if( ipo1 == ixy_end ){ ipo1 = ixy_stt; }
      const double p0[2] = {aXY[ipo0].p[0],aXY[ipo0].p[1]};
      const double p1[2] = {aXY[ipo1].p[0],aXY[ipo1].p[1]};
      const double area0 = delfem2::Area_Tri2(co,codir,p0);
      const double area1 = delfem2::Area_Tri2(co,p1,codir);
      double r1 =  area0 / (area0+area1);
      double r0 =  area1 / (area0+area1);
      if( fabs(area0+area1) < 1.0e-20 ){
        is_fail = true;
        break;
      }
      if( fabs(r0) < 1.0e-3 || fabs(r1) < 1.0e-3 ){
        is_fail = true;
        break;
      }
      if( r0*r1 < 0 ){ continue; }
      double po2[2] = {
        r0*aXY[ipo0].p[0]+r1*aXY[ipo1].p[0],
        r0*aXY[ipo0].p[1]+r1*aXY[ipo1].p[1] };
      //      const double area2 = TriArea2D(co,codir,po2);
      double d = (po2[0]-co[0])*dir[0] + (po2[1]-co[1])*dir[1];
      if( d > 0 ) inum_cross++;
    }
    if( is_fail ){ continue; }
    if( inum_cross%2 == 0 ){ return false; }
    else if( inum_cross%2 == 1 ){ return true; }
  }
  return false;
}

template <typename T>
bool delfem2::CheckInputBoundaryForTriangulation(
    const std::vector<int>& loopIP_ind,
    const std::vector<CVec2<T>>& aXY)
{
  // ------------------------------------
  // enter Input check section

  assert( loopIP_ind.size() >= 2 );
  const size_t nloop = loopIP_ind.size()-1;
  
  { // make sure every loop has at least 3 points
    for(unsigned int iloop=0;iloop<nloop;iloop++){
      if( loopIP_ind[iloop+1]-loopIP_ind[iloop] < 3 ) return false;
    }
  }
  {
    // ------------------------------
    // check inclusion of loops
    for(unsigned int iloop=1;iloop<nloop;iloop++){
      for(int ipo=loopIP_ind[iloop];ipo<loopIP_ind[iloop+1];ipo++){
        const double pi[2] = {aXY[ipo].p[0],aXY[ipo].p[1]};
        if( !IsInclude_Loop(pi,
                            loopIP_ind[0],loopIP_ind[1],
                            aXY) ) return false;
      }
    }
    // check inclusion
    for(unsigned int iloop=1;iloop<nloop;iloop++){
      for(unsigned int jloop=0;jloop<nloop;jloop++){
        if( iloop == jloop ) continue;
        for(int jpo=loopIP_ind[jloop];jpo<loopIP_ind[jloop+1];jpo++){
          const double pj[2] = {aXY[jpo].p[0],aXY[jpo].p[1]};
          if( IsInclude_Loop(pj,
                             loopIP_ind[iloop],loopIP_ind[iloop+1],
                             aXY) ) return false;
        }
      }
    }
  }
  { // check intersection
    bool is_intersect = false;
    for(unsigned int iloop=0;iloop<nloop;iloop++){
      const int nei = loopIP_ind[iloop+1]-loopIP_ind[iloop];
      for(int ie=0;ie<nei;ie++){
        const int i0 = loopIP_ind[iloop] + (ie+0)%nei;
        const int i1 = loopIP_ind[iloop] + (ie+1)%nei;
        const double pi0[2] = {aXY[i0].p[0],aXY[i0].p[1]};
        const double pi1[2] = {aXY[i1].p[0],aXY[i1].p[1]};
        const double xmax_i = ( pi1[0] > pi0[0] ) ? pi1[0] : pi0[0];
        const double xmin_i = ( pi1[0] < pi0[0] ) ? pi1[0] : pi0[0];
        const double ymax_i = ( pi1[1] > pi0[1] ) ? pi1[1] : pi0[1];
        const double ymin_i = ( pi1[1] < pi0[1] ) ? pi1[1] : pi0[1];
        for(int je=ie+1;je<nei;je++){
          const int j0 = loopIP_ind[iloop] + (je+0)%nei;
          const int j1 = loopIP_ind[iloop] + (je+1)%nei;
          if( i0 == j0 || i0 == j1 || i1 == j0 || i1 == j1 ){ continue; }
          const double pj0[2] = {aXY[j0].p[0],aXY[j0].p[1]};
          const double pj1[2] = {aXY[j1].p[0],aXY[j1].p[1]};
          const double xmax_j = ( pj1[0] > pj0[0] ) ? pj1[0] : pj0[0];
          const double xmin_j = ( pj1[0] < pj0[0] ) ? pj1[0] : pj0[0];
          const double ymax_j = ( pj1[1] > pj0[1] ) ? pj1[1] : pj0[1];
          const double ymin_j = ( pj1[1] < pj0[1] ) ? pj1[1] : pj0[1];
          if( xmin_j > xmax_i || xmax_j < xmin_i ){ continue; }
          if( ymin_j > ymax_i || ymax_j < ymin_i ){ continue; }
          if( IsCrossLines(pi0,pi1,  pj0,pj1) ){
            is_intersect = true;
            break;
          }
        }
        if( is_intersect ) break;
        for(unsigned int jloop=iloop+1;jloop<nloop;jloop++){
          const int nbar_j = loopIP_ind[jloop+1]-loopIP_ind[jloop];
          for(int jbar=0;jbar<nbar_j;jbar++){
            const int jpo0 = loopIP_ind[jloop] + jbar;
            int jpo1 = loopIP_ind[jloop] + jbar+1;
            if( jbar == nbar_j-1 ){ jpo1 = loopIP_ind[jloop]; }
            const double pj0[2] = {aXY[jpo0].p[0],aXY[jpo0].p[1]};
            const double pj1[2] = {aXY[jpo1].p[1],aXY[jpo1].p[1]};
            const double xmax_j = ( pj1[0] > pj0[0] ) ? pj1[0] : pj0[0];
            const double xmin_j = ( pj1[0] < pj0[0] ) ? pj1[0] : pj0[0];
            const double ymax_j = ( pj1[1] > pj0[1] ) ? pj1[1] : pj0[1];
            const double ymin_j = ( pj1[1] < pj0[1] ) ? pj1[1] : pj0[1];
            if( xmin_j > xmax_i || xmax_j < xmin_i ) continue;  // åçˆÇ™Ç†ÇËÇ¶Ç»Ç¢ÉpÉ^Å[ÉìÇèúäO
            if( ymin_j > ymax_i || ymax_j < ymin_i ) continue;  // è„Ç…ìØÇ∂
            if( IsCrossLines(pi0,pi1,  pj0,pj1) ){
              is_intersect = true;
              break;
            }
          }
          if( is_intersect ) break;
        }
        if( is_intersect ) break;
      }
      if( is_intersect ) break;
    }
    if( is_intersect ) return false;
  }
  // end of input check section
  return true;
}
#ifdef DFM2_STATIC_LIBRARY
template bool delfem2::CheckInputBoundaryForTriangulation(
    const std::vector<int>& loopIP_ind, const std::vector<CVec2d>& aXY);
#endif

// -------------------------------------------

template <typename T>
std::vector<delfem2::CVec2<T>> delfem2::Polygon_Invert(
    const std::vector<CVec2<T>>& aP)
{
  std::vector<CVec2<T>> res;
  for(int ip=(int)aP.size()-1;ip>=0;--ip){
    res.push_back(aP[ip]);
  }
  return res;
}

template <typename T>
std::vector<double> delfem2::XY_Polygon(
    const std::vector<CVec2<T>>& aP)
{
  std::vector<double> res;
  res.reserve(aP.size()*2);
  for(const auto & ip : aP){
    res.push_back(ip.p[0]);
    res.push_back(ip.p[1]);
  }
  return res;
}

template <typename T>
void delfem2::FixLoopOrientation(
    std::vector<int>& loopIP,
    const std::vector<int>& loopIP_ind,
    const std::vector<CVec2<T>>& aXY)
{
  const std::vector<int> loop_old = loopIP;
  assert( loopIP_ind.size()>1 );
  const size_t nloop = loopIP_ind.size()-1;
  int ivtx0 = 0;
  for(unsigned int iloop=0;iloop<nloop;iloop++){
    double area_loop = 0;
    { // area of this loop
      CVec2<T> vtmp(0,0);
      const int nbar = loopIP_ind[iloop+1]-loopIP_ind[iloop];
      for(int ibar=0;ibar<nbar;ibar++){
        const int iipo0 = loopIP_ind[iloop]+(ibar+0)%nbar;
        const int iipo1 = loopIP_ind[iloop]+(ibar+1)%nbar;
        const int ipo0 = loop_old[iipo0];
        const int ipo1 = loop_old[iipo1];
        area_loop += Area_Tri(vtmp, aXY[ipo0], aXY[ipo1]);
      }
    }
    const int nbar0 = loopIP_ind[iloop+1]-loopIP_ind[iloop];
    if( (area_loop > 0) == (iloop == 0) ){ // outer loop
      for(int ibar=0;ibar<nbar0;ibar++){
        const int iipo = loopIP_ind[iloop] + ibar;
        const int ipo = loop_old[iipo];
        loopIP[ivtx0] = ipo;
        ivtx0++;
      }
    }
    else{
      for(int ibar=0;ibar<nbar0;ibar++){ // inner loop
        const int iipo = loopIP_ind[iloop+1] - 1 - ibar;
        const int ipo = loop_old[iipo];
        loopIP[ivtx0] = ipo;
        ivtx0++;
      }
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::FixLoopOrientation(
    std::vector<int>& loopIP,
    const std::vector<int>& loopIP_ind, const std::vector<CVec2d>& aXY);
#endif



