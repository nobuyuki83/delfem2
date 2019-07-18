#include <stdio.h>
#include <deque>
#include <set>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#include "delfem2/mat3.h"
#include "delfem2/mshtopo.h"
#include "delfem2/msh.h"
#include "delfem2/cad3d.h"

#include "delfem2/v23m3q.h"
#include "delfem2/dyntri_v2.h"

#include "delfem2/gl_funcs.h"
#include "delfem2/gl_color.h"
#include "delfem2/gl_v23q.h"


void CCad3D_Vertex::Draw(bool is_selected, int ielem, double view_height) const
{
  ::glDisable(GL_LIGHTING);
  if( is_selected ){
    double s = view_height*0.3;
    ::glDisable(GL_LIGHTING);
    if( !isConst[0] ){
      ::glColor3d(1, 0, 0);
      DrawArrow(pos, CVector3(+s, 0, 0));
      DrawArrow(pos, CVector3(-s, 0, 0));
    }
    if( !isConst[1] ){
      ::glColor3d(0, 1, 0);
      DrawArrow(pos, CVector3(0, +s, 0));
      DrawArrow(pos, CVector3(0, -s, 0));
    }
    if( !isConst[2] ){
      ::glColor3d(0, 0, 1);
      DrawArrow(pos, CVector3(0, 0, +s));
      DrawArrow(pos, CVector3(0, 0, -s));
    }
  }
  else{
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0,0,1);
    ::glMatrixMode(GL_MODELVIEW);
    ::glPushMatrix();
    ::myGlTranslate(pos);
    ::glDisable(GL_CULL_FACE);
    ::glScaled(view_height*0.03, view_height*0.03, view_height*0.03);
    DrawSphere(32, 32);
    ::glPopMatrix();
  }
  ////
  ::myGlColorDiffuse(CColor::Red());
//  DrawArrow(pos, +norm*view_height*0.1);
//  ::DrawArrow(pos, -norm*0.15);
}

void CCad3D_Edge::DrawLine(bool is_picked, double view_height) const
{
  ::glDisable(GL_LIGHTING);
  if( is_picked ){
    ::glColor3d(1,1,0);
  }
  else if( is_sim ){
    ::glColor3d(0,0,0);
  }
  else{
    ::glColor3d(0.8,0.8,0.8);
  }
  ::glLineWidth(2);
  ::glBegin(GL_LINE_STRIP);
  for(int ip=0;ip<aP.size();++ip){ ::myGlVertex(aP[ip]); }
  ::glEnd();
}

void CCad3D_Edge::DrawHandler(int ielem_picked, double view_height) const
{
  ::glColor3d(0,1,0);
  DrawCylinder(p0, q0, view_height*0.01);
  DrawCylinder(p1, q1, view_height*0.01);
  {
    if( ielem_picked == 2 ){ ::glColor3d(1,0.0,0); }
    else{                    ::glColor3d(0,0.9,0); }
    ::glMatrixMode(GL_MODELVIEW);
    ::glPushMatrix();
    ::myGlTranslate(q0);
    ::glScaled(view_height*0.02, view_height*0.02, view_height*0.02);
    DrawSphere(32, 32);
    ::glPopMatrix();
  }
  {
    if( ielem_picked == 3 ){ ::glColor3d(1,0.0,0); }
    else{                    ::glColor3d(0,0.9,0); }
    ::glMatrixMode(GL_MODELVIEW);
    ::glPushMatrix();
    ::myGlTranslate(q1);
    ::glScaled(view_height*0.02, view_height*0.02, view_height*0.02);
    DrawSphere(32, 32);
    ::glPopMatrix();
  }
  ::glEnd();
}


bool CCad3D_Edge::isPick(double& ratio, const CVector2& sp0, const float mMV[16], const float mPj[16]) const
{
  const int np = (int)aP.size();
  for(int ie=0;ie<np-1;ie++){
    int ip0 = ie;
    int ip1 = ie+1;
    CVector3 p0 = aP[ip0];
    CVector3 p1 = aP[ip1];
    CVector2 s0 = screenXYProjection(p0, mMV, mPj);
    CVector2 s1 = screenXYProjection(p1, mMV, mPj);
    double dist = GetDist_LineSeg_Point(sp0,s0,s1);
    if( dist < 0.03 ){
      ratio = (ip0+0.5)/(np-1.0);
      return true;
    }
  }
  return false;
}

void FaceCenterNormal(
  CVector3& cg, 
  CVector3& nf,
  const std::vector< std::pair<int,bool> >& aIE,
  const std::vector<CCad3D_Edge>& aEdge)
{
  const  int nIE = aIE.size();
  cg.SetZero();
  double len_tot = 0.0;
  for (int iie = 0; iie<nIE; ++iie){
    int ie0 = aIE[(iie+0)%nIE].first;
    bool dir0 = aIE[(iie+0)%nIE].second;    
    CVector3 pA = dir0 ? aEdge[ie0].p0 : aEdge[ie0].p1;
    CVector3 pB = dir0 ? aEdge[ie0].p1 : aEdge[ie0].p0;
    double lenAB = Distance(pA, pB);
    cg += (pA+pB)*(0.5*lenAB);
    len_tot += lenAB;
  }
  cg /= len_tot;
  ///////
  nf.SetZero();
  for (int iie = 0; iie<nIE; ++iie){
    int ie0 = aIE[(iie+0)%nIE].first;
//    int ie1 = aIE[(iie+1)%nIE].first;
    bool dir0 = aIE[(iie+0)%nIE].second;
    CVector3 pA = dir0 ? aEdge[ie0].p0 : aEdge[ie0].p1;
    CVector3 pB = dir0 ? aEdge[ie0].p1 : aEdge[ie0].p0;
    nf += ((pA-cg)^(pB-cg));
  }
  nf.SetNormalizedVector();
}

///////////////////////////////////////

void CCad3D_Face::Initialize
(const std::vector<CCad3D_Vertex>& aVertex,
 const std::vector<CCad3D_Edge>& aEdge,
 double elen)
{
  aPInfo.resize(0);
  std::vector<double> aXYZ_B0;
  std::vector<double> aXYZ_B1;
  const int ne = (int)aIE.size();
  for(int iie=0;iie<aIE.size();++iie){
    int ie0 = aIE[iie].first;
    assert( ie0>=0 && ie0<aEdge.size() );
    const CCad3D_Edge& e0 = aEdge[ie0];
    const bool dir0 = aIE[iie].second;
    int iv0 = (dir0) ? e0.iv0 : e0.iv1;
    {
      CVector3 p0 = (dir0) ? e0.p0 : e0.p1;
      aXYZ_B1.push_back(p0.x);
      aXYZ_B1.push_back(p0.y);
      aXYZ_B1.push_back(p0.z);
    }
    const int nep = (int)e0.aP.size();
    for(int iep=0;iep<nep-1;++iep){
      int iep0 = (dir0) ? iep : nep-1-iep;
      double ratio = (double)iep0/(nep-1.0);
      CVector3 pep = (1-ratio)*e0.p0 + ratio*e0.p1;
      aXYZ_B0.push_back(pep.x);
      aXYZ_B0.push_back(pep.y);
      aXYZ_B0.push_back(pep.z);
      CFacePointInfo pinfo;
      if( iep==0 ){
        pinfo.itype = 0;
        pinfo.iv = iv0;
      }
      else{
        pinfo.itype = 1;
        pinfo.iv = -1;
      }
      pinfo.ie = ie0;
      pinfo.iep = iep0;
      aPInfo.push_back(pinfo);
    }
    { // for debug
      int iie1 = (iie+ne-1)%ne; // back
      int ie1 = aIE[iie1].first;
      assert( ie1>=0 && ie1<aEdge.size() );
      const CCad3D_Edge& e1 = aEdge[ie1];
      bool dir1 = aIE[iie1].second;
      int iv1 = (dir1) ? e1.iv1 : e1.iv0;
      assert( iv0 == iv1 );
    }
  }
  ////
  CVector3 cg, norm; FaceCenterNormal(cg,norm, aIE, aEdge);
  CVector3 axis_x, axis_y; GetVertical2Vector(norm, axis_x, axis_y);
  std::vector<double> aXY_B0;
  for(int ixyz=0;ixyz<aXYZ_B0.size()/3;++ixyz){
    CVector3 p(aXYZ_B0[ixyz*3+0],aXYZ_B0[ixyz*3+1],aXYZ_B0[ixyz*3+2]);
    aXY_B0.push_back((p-cg)*axis_x);
    aXY_B0.push_back((p-cg)*axis_y);
  }
  std::vector<double> aXY_B1;
  for(int ixyz=0;ixyz<aXYZ_B1.size()/3;++ixyz){
    CVector3 p(aXYZ_B1[ixyz*3+0],aXYZ_B1[ixyz*3+1],aXYZ_B1[ixyz*3+2]);
    aXY_B1.push_back((p-cg)*axis_x);
    aXY_B1.push_back((p-cg)*axis_y);
  }
  std::vector<double> aXY_out;
  {
    std::vector< std::vector<double> > aaXY;
    aaXY.push_back(aXY_B0);
    /////
    std::vector<int> loopIP_ind,loopIP;
    std::vector<CVector2> aVec2;
    double elen = 0.05;
    {
      JArray_FromVecVec_XY(loopIP_ind,loopIP, aVec2,
                           aaXY);
      if( !CheckInputBoundaryForTriangulation(loopIP_ind,aVec2) ){
        return;
      }
      FixLoopOrientation(loopIP,
                         loopIP_ind,aVec2);
    }
    {
      std::vector<CEPo2> aPo2D;
      std::vector<ETri> aETri;
      Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                     loopIP_ind,loopIP);
      if( elen > 1.0e-10 ){
        CInputTriangulation_Uniform param(1.0);
        std::vector<int> aFlgPnt(aVec2.size()), aFlgTri(aETri.size());
        MeshingInside(aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
                      aVec2.size(), 0, elen, param);
      }
      MeshTri2D_Export(aXY_out,aTri, aVec2,aETri);
    }
  }
  const int nxy_bound = (int)aXY_B0.size()/2;
  for(int ip=nxy_bound;ip<aXY_out.size()/2;++ip){
    double x0 = aXY_out[ip*2+0];
    double y0 = aXY_out[ip*2+1];
    CFacePointInfo pinfo;
    pinfo.itype = 2;
    pinfo.aW0.resize(aXY_B0.size()/2);
    pinfo.aW1.resize(aXY_B1.size()/2);
    MeanValueCoordinate2D(pinfo.aW0.data(),x0,y0,aXY_B0.data(),aXY_B0.size()/2);
    MeanValueCoordinate2D(pinfo.aW1.data(),x0,y0,aXY_B1.data(),aXY_B1.size()/2);
    aPInfo.push_back(pinfo);
  }
  MovePoints(aVertex,aEdge);
}

void CCad3D_Face::MovePoints
(const std::vector<CCad3D_Vertex>& aVertex,
 const std::vector<CCad3D_Edge>& aEdge)
{
  aXYZ.resize(aPInfo.size()*3);
  for(int ip=0;ip<aPInfo.size();++ip){
    if( aPInfo[ip].itype == 0 ){
      int iv0 = aPInfo[ip].iv;
      aPInfo[ip].n = aVertex[iv0].norm;
      aXYZ[ip*3+0] = aVertex[iv0].pos.x;
      aXYZ[ip*3+1] = aVertex[iv0].pos.y;
      aXYZ[ip*3+2] = aVertex[iv0].pos.z;
    }
    else if( aPInfo[ip].itype == 1 ){
      int ie0 = aPInfo[ip].ie;
      int iep0 = aPInfo[ip].iep;
      CVector3 ne = aEdge[ie0].getNorm();
      CVector3 te = aEdge[ie0].GetTangentInEdge((double)iep0/(aEdge[ie0].aP.size()-1));
      CVector3 nep = ne^te;
      nep.SetNormalizedVector();
      aPInfo[ip].n = nep;
      aXYZ[ip*3+0] = aEdge[ie0].aP[iep0].x;
      aXYZ[ip*3+1] = aEdge[ie0].aP[iep0].y;
      aXYZ[ip*3+2] = aEdge[ie0].aP[iep0].z;
    }
  }
  //////
  if( aIE.size() == 3 ){
    CVector3 aP[9] = {
      aEdge[aIE[0].first].getVtxPos(aIE[0].second,0),
      aEdge[aIE[0].first].getVtxPos(aIE[0].second,1),
      aEdge[aIE[0].first].getVtxPos(aIE[0].second,2),
      ////
      aEdge[aIE[1].first].getVtxPos(aIE[1].second,0),
      aEdge[aIE[1].first].getVtxPos(aIE[1].second,1),
      aEdge[aIE[1].first].getVtxPos(aIE[1].second,2),
      ////
      aEdge[aIE[2].first].getVtxPos(aIE[2].second,0),
      aEdge[aIE[2].first].getVtxPos(aIE[2].second,1),
      aEdge[aIE[2].first].getVtxPos(aIE[2].second,2),
    };
    for(int ip=0;ip<aPInfo.size();++ip){
      if( aPInfo[ip].itype != 2 ){ continue; }
      const std::vector<double>& aW1 = aPInfo[ip].aW1;
      assert( aW1.size() == 3 );
      CVector3 p =  getPointCoonsTri_CubicBezierEdge(aW1[0],aW1[1],aW1[2],aP);
      aXYZ[ip*3+0] = p.x;
      aXYZ[ip*3+1] = p.y;
      aXYZ[ip*3+2] = p.z;
    }
  }
  else if( aIE.size() == 4 ){
    CVector3 aP[12] = {
      aEdge[aIE[0].first].getVtxPos(aIE[0].second,0),
      aEdge[aIE[0].first].getVtxPos(aIE[0].second,1),
      aEdge[aIE[0].first].getVtxPos(aIE[0].second,2),
      ////
      aEdge[aIE[1].first].getVtxPos(aIE[1].second,0),
      aEdge[aIE[1].first].getVtxPos(aIE[1].second,1),
      aEdge[aIE[1].first].getVtxPos(aIE[1].second,2),
      ////
      aEdge[aIE[2].first].getVtxPos(aIE[2].second,0),
      aEdge[aIE[2].first].getVtxPos(aIE[2].second,1),
      aEdge[aIE[2].first].getVtxPos(aIE[2].second,2),
      ////
      aEdge[aIE[3].first].getVtxPos(aIE[3].second,0),
      aEdge[aIE[3].first].getVtxPos(aIE[3].second,1),
      aEdge[aIE[3].first].getVtxPos(aIE[3].second,2),
    };
    for(int ip=0;ip<aPInfo.size();++ip){
      if( aPInfo[ip].itype != 2 ){ continue; }
      const std::vector<double>& aW1 = aPInfo[ip].aW1;
      assert( aW1.size() == 4 );
      const double u = aW1[1] + aW1[2];
      const double v = aW1[2] + aW1[3];
//      CVector3 p =  getPointCoons_CubicBezier(u,v,aP);
      CVector3 p =  getPointHermetianQuad(u,v,aP);
      aXYZ[ip*3+0] = p.x;
      aXYZ[ip*3+1] = p.y;
      aXYZ[ip*3+2] = p.z;
    }
  }
  else{
    for(int ip=0;ip<aPInfo.size();++ip){
      if( aPInfo[ip].itype != 2 ){ continue; }
      aXYZ[ip*3+0] = 0;
      aXYZ[ip*3+1] = 0;
      aXYZ[ip*3+2] = 0;
      const std::vector<double>& aW = aPInfo[ip].aW0;
      for(int jp=0;jp<aW.size();++jp){
        aXYZ[ip*3+0] += aW[jp]*aXYZ[jp*3+0];
        aXYZ[ip*3+1] += aW[jp]*aXYZ[jp*3+1];
        aXYZ[ip*3+2] += aW[jp]*aXYZ[jp*3+2];
      }
      /*
      CVector3 pi(aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]);
      for(int jp=0;jp<aW.size();++jp){
        CVector3 pj(aXYZ[jp*3+0],aXYZ[jp*3+1],aXYZ[jp*3+2]);
        const CVector3 nj = aPInfo[jp].n;
        CVector3 dp = ((pj-pi)*nj)*nj*0.8; // control per edge ?
        aXYZ[ip*3+0] += aW[jp]*dp.x;
        aXYZ[ip*3+1] += aW[jp]*dp.y;
        aXYZ[ip*3+2] += aW[jp]*dp.z;
      }
       */
    }
  }
  aNorm.resize(aXYZ.size());
  Normal_MeshTri3D(aNorm.data(),
                   aXYZ.data(), aXYZ.size()/3,
                   aTri.data(), aTri.size()/3);
}

void CCad3D_Face::DrawFace() const
{
  DrawMeshTri3D_FaceNorm(aXYZ, aTri, aNorm);
}

void CCad3D_Face::DrawBackFace() const
{
  ::glEnable(GL_LIGHTING);
  ::myGlColorDiffuse(CColor::Gray(0.9));
  DrawMeshTri3D_FaceNorm_XYsym(aXYZ, aTri);
}

void CCad3D_Face::DrawEdge() const
{
  ::glLineWidth(1);
  DrawMeshTri3D_Edge(aXYZ, aTri);
}


////////////////////////////////////////////////////

int AddPointEdge
(int ie_div, double ratio_edge,
 std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen)
{
  if( ie_div < 0 || ie_div >= aEdge.size() ) return -1;
  if (ratio_edge < 0.01||ratio_edge > 0.99) return -1;
  const int iv_new = (int)aVertex.size();
  {
    CVector3 nv;
    for (int ifc = 0; ifc<aFace.size(); ++ifc){       
      const CCad3D_Face& fc = aFace[ifc];       
      for (int ie = 0; ie<fc.aIE.size(); ++ie){
        if (fc.aIE[ie].first!=ie_div) continue;         
        CVector3 cg, nf;
        FaceCenterNormal(cg, nf, fc.aIE, aEdge);         
        nv += nf;        
      }
    }
    nv.SetNormalizedVector();
    //////
    CVector3 p = aEdge[ie_div].GetPosInEdge(ratio_edge);
    CCad3D_Vertex v(p);
    {
      int ien = aEdge[ie_div].inorm;
      v.isConst[ien] = true;
    }
    v.norm = nv;
    aVertex.push_back(v);
  }
  const int iv0 = aEdge[ie_div].iv0;
  const int iv1 = aEdge[ie_div].iv1;
  {
    aEdge[ie_div].iv0 = iv0;
    aEdge[ie_div].iv1 = iv_new;
    aEdge[ie_div].Initialize(aVertex,elen);
  }
  const int ie_new = (int)aEdge.size();
  aEdge.push_back( CCad3D_Edge(iv_new,iv1,aEdge[ie_div].is_sim,aEdge[ie_div].inorm) );
  aEdge[ie_new].Initialize(aVertex,elen);
  
  for(int ifc=0;ifc<aFace.size();++ifc){
    CCad3D_Face& fc = aFace[ifc];
    for(int ie=0;ie<fc.aIE.size();++ie){
      if(fc.aIE[ie].first!=ie_div) continue;
      if(fc.aIE[ie].second){
        fc.aIE.insert(fc.aIE.begin()+ie+1,std::make_pair(ie_new,true));
      }
      else{
        fc.aIE.insert(fc.aIE.begin()+ie,std::make_pair(ie_new,false));
      }
      fc.Initialize(aVertex, aEdge, elen);
      break;
    }
  }
  return iv_new;
}

void ConectEdge
(int iv0, int iv1, int iface_div, int inorm_new,
 std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen)
{
  if( iface_div < 0 || iface_div >= (int)aFace.size() ) return;
  int iie0 = aFace[iface_div].findIIE_CP(iv0,aEdge);
  int iie1 = aFace[iface_div].findIIE_CP(iv1,aEdge);
  if( iie0 == -1 || iie1 == -1 ) return;
  { // move iv0 and iv1
    CVector3 p0 = aVertex[iv0].pos;
    CVector3 p1 = aVertex[iv1].pos;
    CVector3 mid = (p0+p1)*0.5;
    CVector3 n(0,0,0); n[inorm_new] =1;
    aVertex[iv0].pos = p0-((p0-mid)*n)*n;
    aVertex[iv1].pos = p1-((p1-mid)*n)*n;
  }
  if( inorm_new >= 0 && inorm_new < 3 ){
    aVertex[iv0].isConst[inorm_new] = true;
    aVertex[iv1].isConst[inorm_new] = true;
  }
  for(int iie=0;iie<aFace[iface_div].aIE.size();++iie){
    int ie0 = aFace[iface_div].aIE[iie].first;
    int jv0 = aEdge[ie0].iv0;
    int jv1 = aEdge[ie0].iv1;
    if( (jv0==iv0&&jv1==iv1) || (jv0==iv1&&jv1==iv0) ) return;
  }
  ////////////
  const int ie_new = (int)aEdge.size();
  aEdge.push_back( CCad3D_Edge(iv0,iv1,false,inorm_new) );
  aEdge[ie_new].Initialize(aVertex,elen);
  
  const std::vector< std::pair<int,bool> > aIE = aFace[iface_div].aIE;
  const int nie = (int)aIE.size();
  { // modify exisiting
    std::vector< std::pair<int,bool> > aIE0;
    aIE0.push_back( std::make_pair(ie_new,true) );
    for(int iie=iie1;iie%nie!=iie0;++iie){
      aIE0.push_back( aIE[iie%nie] );
    }
    aFace[iface_div].aIE = aIE0;
    aFace[iface_div].Initialize(aVertex, aEdge, elen);
  }
  { // make new
    std::vector< std::pair<int,bool> > aIE0;
    aIE0.push_back( std::make_pair(ie_new,false) );
    for(int iie=iie0;iie%nie!=iie1;++iie){
      aIE0.push_back( aIE[iie%nie] );
    }
    CCad3D_Face face(aIE0);
    face.Initialize(aVertex, aEdge, elen);
    aFace.push_back(face);
  }
  for(int ie=0;ie<aEdge.size();++ie){
    aEdge[ie].MovePoints(aVertex);
  }
  for(int ifc=0;ifc<aFace.size();++ifc){
    aFace[ifc].MovePoints(aVertex,aEdge);
  }
}

void MakeItSmooth
(std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace)
{
  for(int iv=0;iv<aVertex.size();++iv){
    aVertex[iv].norm.SetZero();
  }
  for(int ifc=0;ifc<aFace.size();++ifc){
    const std::vector< std::pair<int,bool> >& aIE = aFace[ifc].aIE;
    int nIE = (int)aIE.size();
    CVector3 nf,cg; FaceCenterNormal(cg,nf,aIE,aEdge);
    nf.SetNormalizedVector();
    for(int iie=0;iie<nIE;++iie){
      int ie0 = aIE[iie].first;
      bool dir0 = aIE[iie].second;
      int ipA = dir0 ? aEdge[ie0].iv0 : aEdge[ie0].iv1;
      aVertex[ipA].norm += nf;
    }
  }
  for(int iv=0;iv<aVertex.size();++iv){
//    if( aVertex[iv].isConst[0] ){ aVertex[iv].norm.x = 0; }
//    if( aVertex[iv].isConst[1] ){ aVertex[iv].norm.y = 0; }
//    if( aVertex[iv].isConst[2] ){ aVertex[iv].norm.z = 0; }
    if( aVertex[iv].norm.Length() < 0.1 ){
      aVertex[iv].norm.SetZero();
      continue;
    }
    aVertex[iv].norm.SetNormalizedVector();
  }
  for(int ie=0;ie<aEdge.size();ie++){
    aEdge[ie].MovePoints(aVertex); // ie0+0
  }
  for(int ifc=0;ifc<aFace.size();ifc++){
    aFace[ifc].MovePoints(aVertex,aEdge); // ie0+0
  }
}

void findEdgeGroup
(std::vector< std::pair<int,bool> >& aIE,
 int iedge0,
 std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge)
{
  aIE.clear();
  if( iedge0 < 0 || iedge0 >= aEdge.size() ) return;
  std::deque< std::pair<int,bool> > deqIE;
  deqIE.push_back( std::make_pair(iedge0,true) );
  bool is_loop = false;
  for(;;){
    const int ie0 = deqIE.back().first;
    int iv0 = deqIE.back().second ? aEdge[ie0].iv1 : aEdge[ie0].iv0;
    int ine0 = aEdge[ie0].inorm;
    if( iv0 == aEdge[iedge0].iv0 ){ is_loop = true; break; }
    int ndeqIE = (int)deqIE.size(); // prev
    for(int ie=0;ie<aEdge.size();++ie){
      if( ie == ie0 ) continue;
      if(      aEdge[ie].iv0 == iv0 && aEdge[ie].inorm == ine0){
        deqIE.push_back( std::make_pair(ie,true ) );
        break;
      }
      else if( aEdge[ie].iv1 == iv0 && aEdge[ie].inorm == ine0 ){
        deqIE.push_back( std::make_pair(ie,false) );
        break;
      }
    }
    if( deqIE.size() == ndeqIE ) break; // couldn't find new one
  }
  if( is_loop ){ aIE.assign(deqIE.begin(), deqIE.end()); return; }
  ///
  for(;;){
    const int ie0 = deqIE.front().first;
    int iv0 = deqIE.front().second ? aEdge[ie0].iv0 : aEdge[ie0].iv1;
    int ine0 = aEdge[ie0].inorm;
    assert( iv0 != aEdge[iedge0].iv1 ); // this should not be loop
    int ndeqIE = (int)deqIE.size(); // prev
    for(int ie=0;ie<aEdge.size();++ie){
      if( ie == ie0 ) continue;
      if(      aEdge[ie].iv0 == iv0 && aEdge[ie].inorm == ine0){
        deqIE.push_front( std::make_pair(ie,false ) );
        break;
      }
      else if( aEdge[ie].iv1 == iv0 && aEdge[ie].inorm == ine0 ){
        deqIE.push_front( std::make_pair(ie,true) );
        break;
      }
    }
    if( deqIE.size() == ndeqIE ) break; // couldn't find new one
  }
  aIE.assign(deqIE.begin(), deqIE.end());
}

void AddSphere_ZSym
(std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen)
{
  int icp0 = (int)aVertex.size();
  aVertex.push_back( CCad3D_Vertex(CVector3(-1, 0, 0)) ); // icp0+0
  aVertex.push_back( CCad3D_Vertex(CVector3( 0,+1, 0)) ); // icp0+1
  aVertex.push_back( CCad3D_Vertex(CVector3(+1, 0, 0)) ); // icp0+2
  aVertex.push_back( CCad3D_Vertex(CVector3( 0,-1, 0)) ); // icp0+3
  aVertex.push_back( CCad3D_Vertex(CVector3( 0, 0,+1)) ); // icp0+4
  ////
  ////
  int ie0 = (int)aEdge.size();
  CCad3D_Edge e0(icp0+0,icp0+1,true,2);
  CCad3D_Edge e1(icp0+1,icp0+2,true,2);
  CCad3D_Edge e2(icp0+2,icp0+3,true,2);
  CCad3D_Edge e3(icp0+3,icp0+0,true,2);
  CCad3D_Edge e4(icp0+4,icp0+0,false,1);
  CCad3D_Edge e5(icp0+4,icp0+1,false,0);
  CCad3D_Edge e6(icp0+2,icp0+4,false,1);
  CCad3D_Edge e7(icp0+3,icp0+4,false,0);
  aEdge.push_back(e0);
  aEdge.push_back(e1);
  aEdge.push_back(e2);
  aEdge.push_back(e3);
  aEdge.push_back(e4);
  aEdge.push_back(e5);
  aEdge.push_back(e6);
  aEdge.push_back(e7);
  /////
  int ifc0 = (int)aFace.size();
  { // face041
    std::vector< std::pair<int,bool> > aIE;
    aIE.push_back( std::make_pair(ie0+0,false) );
    aIE.push_back( std::make_pair(ie0+4,false ) );
    aIE.push_back( std::make_pair(ie0+5,true) );
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face041
    std::vector< std::pair<int,bool> > aIE;
    aIE.push_back( std::make_pair(ie0+1,false) );
    aIE.push_back( std::make_pair(ie0+5,false) );
    aIE.push_back( std::make_pair(ie0+6,false) );
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face041
    std::vector< std::pair<int,bool> > aIE;
    aIE.push_back( std::make_pair(ie0+2,false) );
    aIE.push_back( std::make_pair(ie0+6,true ) );
    aIE.push_back( std::make_pair(ie0+7,false) );
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face041
    std::vector< std::pair<int,bool> > aIE;
    aIE.push_back( std::make_pair(ie0+3,false) );
    aIE.push_back( std::make_pair(ie0+7,true ) );
    aIE.push_back( std::make_pair(ie0+4,true) );
    aFace.push_back(CCad3D_Face(aIE));
  }
  {
    for(int iv=0;iv<aVertex.size();++iv){
      aVertex[iv].isConst[0] = false;
      aVertex[iv].isConst[1] = false;
      aVertex[iv].isConst[2] = false;
    }
    for(int ie=0;ie<aEdge.size();++ie){
      int iv0 = aEdge[ie].iv0;
      int iv1 = aEdge[ie].iv1;
      int inorm = aEdge[ie].inorm;
      if( inorm < 0 || inorm >= 3 ){ continue; }
      aVertex[iv0].isConst[inorm] = true;
      aVertex[iv1].isConst[inorm] = true;
    }
  }
  elen = 0.1;
  for(int ie=ie0;ie<ie0+8;ie++){
    aEdge[ie].Initialize(aVertex,elen); // ie0+0
  }
  for(int ifc=ifc0;ifc<ifc0+4;ifc++){
    aFace[ifc].Initialize(aVertex,aEdge, elen); // ie0+0
  }
  MakeItSmooth(aVertex,aEdge,aFace);
}

void AddTorus_XSym
(std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen)
{
  int icp0 = (int)aVertex.size();
  aVertex.push_back( CCad3D_Vertex(CVector3(0,-1.0, 0.0)) ); // icp0+0
  aVertex.push_back( CCad3D_Vertex(CVector3(0,-0.3, 0.0)) ); // icp0+0
  aVertex.push_back( CCad3D_Vertex(CVector3(0, 0.0,+1.0)) ); // icp0+1
  aVertex.push_back( CCad3D_Vertex(CVector3(0, 0.0,+0.3)) ); // icp0+1
  aVertex.push_back( CCad3D_Vertex(CVector3(0,+1.0, 0.0)) ); // icp0+2
  aVertex.push_back( CCad3D_Vertex(CVector3(0,+0.3, 0.0)) ); // icp0+2
  aVertex.push_back( CCad3D_Vertex(CVector3(0, 0.0,-1.0)) ); // icp0+3
  aVertex.push_back( CCad3D_Vertex(CVector3(0, 0.0,-0.3)) ); // icp0+3
  aVertex[icp0+0].norm = CVector3(0,-1,0);
  aVertex[icp0+1].norm = CVector3(0,+1,0);
  aVertex[icp0+2].norm = CVector3(0,0,+1);
  aVertex[icp0+3].norm = CVector3(0,0,-1);
  aVertex[icp0+4].norm = CVector3(0,+1,0);
  aVertex[icp0+5].norm = CVector3(0,-1,0);
  aVertex[icp0+6].norm = CVector3(0,0,-1);
  aVertex[icp0+7].norm = CVector3(0,0,+1);
  /////
  int ie0 = (int)aEdge.size();
  aEdge.push_back( CCad3D_Edge(icp0+0,icp0+2,true,0) ); // 0
  aEdge.push_back( CCad3D_Edge(icp0+2,icp0+4,true,0) ); // 1
  aEdge.push_back( CCad3D_Edge(icp0+4,icp0+6,true,0) ); // 2
  aEdge.push_back( CCad3D_Edge(icp0+6,icp0+0,true,0) ); // 3
  aEdge.push_back( CCad3D_Edge(icp0+3,icp0+1,true,0) ); // 4
  aEdge.push_back( CCad3D_Edge(icp0+5,icp0+3,true,0) ); // 5
  aEdge.push_back( CCad3D_Edge(icp0+7,icp0+5,true,0) ); // 6
  aEdge.push_back( CCad3D_Edge(icp0+1,icp0+7,true,0) ); // 7
  aEdge.push_back( CCad3D_Edge(icp0+1,icp0+0,false,2) ); // 8
  aEdge.push_back( CCad3D_Edge(icp0+3,icp0+2,false,1) ); // 9
  aEdge.push_back( CCad3D_Edge(icp0+4,icp0+5,false,2) ); // 10
  aEdge.push_back( CCad3D_Edge(icp0+6,icp0+7,false,1) ); // 11
  for(int ie=ie0+8;ie<ie0+12;++ie){
    aEdge[ie].r0 = 1;
    aEdge[ie].r1 = 1;
  }
  for(int ie=ie0;ie<ie0+12;ie++){
    aEdge[ie].Initialize(aVertex,elen); // ie0+0
  }
  /////
  int ifc0 = (int)aFace.size();
  { // face041
    std::vector< std::pair<int,bool> > aIE;
    aIE.push_back( std::make_pair(ie0+0,false) );
    aIE.push_back( std::make_pair(ie0+8,false) );
    aIE.push_back( std::make_pair(ie0+4,false) );
    aIE.push_back( std::make_pair(ie0+9,true ) );
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face041
    std::vector< std::pair<int,bool> > aIE;
    aIE.push_back( std::make_pair(ie0+1,false) );
    aIE.push_back( std::make_pair(ie0+9,false) );
    aIE.push_back( std::make_pair(ie0+5,false) );
    aIE.push_back( std::make_pair(ie0+10,false ) );
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face041
    std::vector< std::pair<int,bool> > aIE;
    aIE.push_back( std::make_pair(ie0+2,false) );
    aIE.push_back( std::make_pair(ie0+10,true) );
    aIE.push_back( std::make_pair(ie0+6,false) );
    aIE.push_back( std::make_pair(ie0+11,false ) );
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face041
    std::vector< std::pair<int,bool> > aIE;
    aIE.push_back( std::make_pair(ie0+3,false) );
    aIE.push_back( std::make_pair(ie0+11,true) );
    aIE.push_back( std::make_pair(ie0+7,false) );
    aIE.push_back( std::make_pair(ie0+8,true ) );
    aFace.push_back(CCad3D_Face(aIE));
  }
  {
    for(int iv=0;iv<aVertex.size();++iv){
      aVertex[iv].isConst[0] = false;
      aVertex[iv].isConst[1] = false;
      aVertex[iv].isConst[2] = false;
    }
    for(int ie=0;ie<aEdge.size();++ie){
      int iv0 = aEdge[ie].iv0;
      int iv1 = aEdge[ie].iv1;
      int inorm = aEdge[ie].inorm;
      if( inorm < 0 || inorm >= 3 ){ continue; }
      aVertex[iv0].isConst[inorm] = true;
      aVertex[iv1].isConst[inorm] = true;
    }
  }
  for(int ifc=ifc0;ifc<ifc0+4;ifc++){
    aFace[ifc].Initialize(aVertex,aEdge, elen); // ie0+0
  }
//  MakeItSmooth(aVertex,aEdge,aFace);
}


void AddSphere_XSym
(std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen)
{
  int icp0 = (int)aVertex.size();
  aVertex.push_back( CCad3D_Vertex(CVector3(0,-1, 0)) ); // icp0+0
  aVertex.push_back( CCad3D_Vertex(CVector3(0, 0,+1)) ); // icp0+1
  aVertex.push_back( CCad3D_Vertex(CVector3(0,+1, 0)) ); // icp0+2
  aVertex.push_back( CCad3D_Vertex(CVector3(0, 0,-1)) ); // icp0+3
  aVertex.push_back( CCad3D_Vertex(CVector3(1, 0, 0)) ); // icp0+4
  ////
  ////
  int ie0 = (int)aEdge.size();
  CCad3D_Edge e0(icp0+0,icp0+1,true,0);
  CCad3D_Edge e1(icp0+1,icp0+2,true,0);
  CCad3D_Edge e2(icp0+2,icp0+3,true,0);
  CCad3D_Edge e3(icp0+3,icp0+0,true,0);
  CCad3D_Edge e4(icp0+4,icp0+0,false,2);
  CCad3D_Edge e5(icp0+4,icp0+1,false,1);
  CCad3D_Edge e6(icp0+2,icp0+4,false,2);
  CCad3D_Edge e7(icp0+3,icp0+4,false,1);
  aEdge.push_back(e0);
  aEdge.push_back(e1);
  aEdge.push_back(e2);
  aEdge.push_back(e3);
  aEdge.push_back(e4);
  aEdge.push_back(e5);
  aEdge.push_back(e6);
  aEdge.push_back(e7);
  /////
  int ifc0 = (int)aFace.size();
  { // face041
    std::vector< std::pair<int,bool> > aIE;
    aIE.push_back( std::make_pair(ie0+0,false) );
    aIE.push_back( std::make_pair(ie0+4,false ) );
    aIE.push_back( std::make_pair(ie0+5,true) );
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face041
    std::vector< std::pair<int,bool> > aIE;
    aIE.push_back( std::make_pair(ie0+1,false) );
    aIE.push_back( std::make_pair(ie0+5,false) );
    aIE.push_back( std::make_pair(ie0+6,false) );
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face041
    std::vector< std::pair<int,bool> > aIE;
    aIE.push_back( std::make_pair(ie0+2,false) );
    aIE.push_back( std::make_pair(ie0+6,true ) );
    aIE.push_back( std::make_pair(ie0+7,false) );
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face041
    std::vector< std::pair<int,bool> > aIE;
    aIE.push_back( std::make_pair(ie0+3,false) );
    aIE.push_back( std::make_pair(ie0+7,true ) );
    aIE.push_back( std::make_pair(ie0+4,true) );
    aFace.push_back(CCad3D_Face(aIE));
  }
  {
    for(int iv=0;iv<aVertex.size();++iv){
      aVertex[iv].isConst[0] = false;
      aVertex[iv].isConst[1] = false;
      aVertex[iv].isConst[2] = false;
    }
    for(int ie=0;ie<aEdge.size();++ie){
      int iv0 = aEdge[ie].iv0;
      int iv1 = aEdge[ie].iv1;
      int inorm = aEdge[ie].inorm;
      if( inorm < 0 || inorm >= 3 ){ continue; }
      aVertex[iv0].isConst[inorm] = true;
      aVertex[iv1].isConst[inorm] = true;
    }
  }
  elen = 0.1;
  for(int ie=ie0;ie<ie0+8;ie++){
    aEdge[ie].Initialize(aVertex,elen); // ie0+0
  }
  for(int ifc=ifc0;ifc<ifc0+4;ifc++){
    aFace[ifc].Initialize(aVertex,aEdge, elen); // ie0+0
  }
  MakeItSmooth(aVertex,aEdge,aFace);
}



void AddCube
(std::vector<CCad3D_Vertex>& aVertex,
std::vector<CCad3D_Edge>& aEdge,
std::vector<CCad3D_Face>& aFace,
double elen)
{
  int iv0 = (int)aVertex.size();
  aVertex.push_back(CCad3D_Vertex(CVector3(-1, -1, -1))); // icp0+0
  aVertex.push_back(CCad3D_Vertex(CVector3(-1, -1, +1))); // icp0+1
  aVertex.push_back(CCad3D_Vertex(CVector3(-1, +1, -1))); // icp0+2
  aVertex.push_back(CCad3D_Vertex(CVector3(-1, +1, +1))); // icp0+3
  aVertex.push_back(CCad3D_Vertex(CVector3(+1, -1, -1))); // icp0+4
  aVertex.push_back(CCad3D_Vertex(CVector3(+1, -1, +1))); // icp0+5
  aVertex.push_back(CCad3D_Vertex(CVector3(+1, +1, -1))); // icp0+6
  aVertex.push_back(CCad3D_Vertex(CVector3(+1, +1, +1))); // icp0+7
  ////
  int ie0 = (int)aEdge.size();
  aEdge.push_back(CCad3D_Edge(iv0+0, iv0+1, false, 0)); // 0
  aEdge.push_back(CCad3D_Edge(iv0+1, iv0+3, false, 0)); // 1
  aEdge.push_back(CCad3D_Edge(iv0+3, iv0+2, false, 0)); // 2
  aEdge.push_back(CCad3D_Edge(iv0+2, iv0+0, false, 0)); // 3
  /////
  aEdge.push_back(CCad3D_Edge(iv0+6, iv0+4, false, 0)); // 4
  aEdge.push_back(CCad3D_Edge(iv0+7, iv0+6, false, 0)); // 5
  aEdge.push_back(CCad3D_Edge(iv0+5, iv0+7, false, 0)); // 6
  aEdge.push_back(CCad3D_Edge(iv0+4, iv0+5, false, 0)); // 7
  /////
  aEdge.push_back(CCad3D_Edge(iv0+0, iv0+4, false, 1)); // 8
  aEdge.push_back(CCad3D_Edge(iv0+5, iv0+1, false, 1)); // 9
  aEdge.push_back(CCad3D_Edge(iv0+2, iv0+6, false, 1)); // 10
  aEdge.push_back(CCad3D_Edge(iv0+7, iv0+3, false, 1)); // 11
  /////  
  int ifc0 = (int)aFace.size();
  { // face0132
    std::vector< std::pair<int, bool> > aIE;
    aIE.push_back(std::make_pair(ie0+0, true));
    aIE.push_back(std::make_pair(ie0+1, true));
    aIE.push_back(std::make_pair(ie0+2, true));
    aIE.push_back(std::make_pair(ie0+3, true));
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face4567
    std::vector< std::pair<int, bool> > aIE;
    aIE.push_back(std::make_pair(ie0+4, false));
    aIE.push_back(std::make_pair(ie0+5, false));
    aIE.push_back(std::make_pair(ie0+6, false));
    aIE.push_back(std::make_pair(ie0+7, false));
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face0451
    std::vector< std::pair<int, bool> > aIE;
    aIE.push_back(std::make_pair(ie0+0, false));
    aIE.push_back(std::make_pair(ie0+8, true));
    aIE.push_back(std::make_pair(ie0+7, true));
    aIE.push_back(std::make_pair(ie0+9, true));
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face041
    std::vector< std::pair<int, bool> > aIE;
    aIE.push_back(std::make_pair(ie0+2,  false));
    aIE.push_back(std::make_pair(ie0+11, false));
    aIE.push_back(std::make_pair(ie0+5,  true));
    aIE.push_back(std::make_pair(ie0+10, false));
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face041     
    std::vector< std::pair<int, bool> > aIE;
    aIE.push_back(std::make_pair(ie0+3,  false));
    aIE.push_back(std::make_pair(ie0+10, true));
    aIE.push_back(std::make_pair(ie0+4,  true));
    aIE.push_back(std::make_pair(ie0+8,  false));
    aFace.push_back(CCad3D_Face(aIE));
  }
  { // face041     
    std::vector< std::pair<int, bool> > aIE;
    aIE.push_back(std::make_pair(ie0+1,  false));
    aIE.push_back(std::make_pair(ie0+9,  false));
    aIE.push_back(std::make_pair(ie0+6,  true));
    aIE.push_back(std::make_pair(ie0+11, true));
    aFace.push_back(CCad3D_Face(aIE));
  }
  {
    for (int iv = iv0; iv<iv0+8; ++iv){
      aVertex[iv].isConst[0] = false;
      aVertex[iv].isConst[1] = false;
      aVertex[iv].isConst[2] = false;
    }
    for (int ie = ie0; ie<ie0+12; ++ie){
      int iv0 = aEdge[ie].iv0;
      int iv1 = aEdge[ie].iv1;
      int inorm = aEdge[ie].inorm;
      if (inorm < 0||inorm>=3){ continue; }
      aVertex[iv0].isConst[inorm] = true;
      aVertex[iv1].isConst[inorm] = true;
    }
  }
  
  elen = 0.1;
  for (int ie = ie0; ie<ie0+12; ie++){
    aEdge[ie].Initialize(aVertex, elen); // ie0+0
  }  
  for (int ifc = ifc0; ifc<ifc0+6; ifc++){
    aFace[ifc].Initialize(aVertex, aEdge, elen); // ie0+0
  }
  
  MakeItSmooth(aVertex, aEdge, aFace);  
}

bool FindFittingPoint
(CVector2& p2d_near,
 CVector2& p2d_norm,
 const CVector2& p2d_org,
 const std::vector<CVector2>& aP2D,
 bool isConstX, bool isConstY,
 double half_view_height)
{
  bool isHit = false;
  if( isConstX &&  isConstY ){ return false; }
  else if( isConstX && !isConstY ){
    for(int iq=0;iq<aP2D.size()-1;++iq){
      CVector2 q0 = aP2D[iq+0];
      CVector2 q1 = aP2D[iq+1];
      if( (q0.x-p2d_org.x)*(q1.x-p2d_org.x) < 0 ){
        p2d_near.x = p2d_org.x;
        p2d_near.y = ((q0+q1)*0.5).y;
        p2d_norm.x = (q1-q0).y;
        p2d_norm.y = (q0-q1).x;
        isHit = true;
        break;
      }
    }
  }
  else if( !isConstX && isConstY ){
    for(int iq=0;iq<aP2D.size()-1;++iq){
      CVector2 q0 = aP2D[iq+0];
      CVector2 q1 = aP2D[iq+1];
      if( (q0.y-p2d_org.y)*(q1.y-p2d_org.y) < 0 ){
        p2d_near.x = ((q0+q1)*0.5).x;
        p2d_near.y = p2d_org.y;
        p2d_norm.x = (q1-q0).y;
        p2d_norm.y = (q0-q1).x;
        isHit = true;
        break;
      }
    }
  }
  else{
    double min_dist = -1;
    for(int iq=0;iq<aP2D.size();++iq){
      double len = (aP2D[iq]-p2d_org).Length();
      if( min_dist < 0 || len < min_dist ){
        min_dist = len;
        p2d_near = aP2D[iq];
        p2d_norm = aP2D[iq]-p2d_org;
        isHit = true;
      }
    }
  }
  if( !isHit ) return false;
  double dist = (p2d_near-p2d_org).Length();
  if( dist > 0.4 ) return false;
  p2d_norm.SetNormalizedVector();
  return true;
}

std::vector<int> getPointsInEdges
(const std::vector< std::pair<int,bool > >& aIE_picked,
 const std::vector<CCad3D_Edge>& aEdge)
{
  std::vector<int> aIP;
  for(int iie=0;iie<aIE_picked.size()+1;++iie){
    int iv0;
    if( iie != aIE_picked.size() ){
      int ie0 = aIE_picked[iie].first;
      iv0 = aIE_picked[iie].second ? aEdge[ie0].iv0 : aEdge[ie0].iv1;
    }
    else{
      int ie0 = aIE_picked[iie-1].first;
      iv0 = aIE_picked[iie-1].second ? aEdge[ie0].iv1 : aEdge[ie0].iv0;
    }
    aIP.push_back(iv0);
  }
  if( aIP.front() == aIP.back() ){ aIP.pop_back(); }
  return aIP;
}

bool MovePointsAlongSketch
(std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 const std::vector<CVector2>& aStroke,
 const std::vector< std::pair<int,bool > >& aIE_picked,
 const CVector3& plane_org, int inorm,
 float mMV[16], float mPj[16], double view_height)
{
  // resampling
  std::vector<CVector2> aStroke1 = Polyline_Resample_Polyline(aStroke,0.025);
  //
  CVector3 plane_nrm(0,0,0); plane_nrm[inorm] = 1;
  CVector3 plane_ex(0,0,0); plane_ex[(inorm+1)%3] = 1;
  CVector3 plane_ey(0,0,0); plane_ey[(inorm+2)%3] = 1;
  std::vector<CVector2> aP2D;
  for(int ist=0;ist<aStroke1.size();++ist){
    CVector2 sp0 = aStroke1[ist];
    CVector3 src = screenUnProjection(CVector3(sp0.x,sp0.y,0), mMV, mPj);
    CVector3 dir = screenUnProjection(CVector3(0,0,1), mMV, mPj);
    CVector3 p = intersection_Plane_Line(plane_org, plane_nrm, src,dir);
    aP2D.push_back(CVector2((p-plane_org)*plane_ex,(p-plane_org)*plane_ey));
  }
  bool is_moved = false;
  std::vector<int> aIP = getPointsInEdges(aIE_picked,aEdge);
  for(int iip=0;iip<aIP.size();++iip){
    int iv0 = aIP[iip];
    CCad3D_Vertex& v = aVertex[iv0];
    CVector2 p2d_org((v.pos-plane_org)*plane_ex, (v.pos-plane_org)*plane_ey);
    const bool isConstX = v.isConst[(inorm+1)%3];
    const bool isConstY = v.isConst[(inorm+2)%3];
    CVector2 p2d_near, p2d_norm;
    bool res = FindFittingPoint(p2d_near,p2d_norm,
                                p2d_org, aP2D, isConstX,isConstY,view_height*0.2);
    if( res ){
      CVector3 p3d_near = plane_org + p2d_near.x*plane_ex + p2d_near.y*plane_ey;
      CVector3 n3d_near = p2d_norm.x*plane_ex + p2d_norm.y*plane_ey;
      v.pos = p3d_near;
      v.norm = n3d_near;
      is_moved = true;
    }
  }
  ////
  if( is_moved ){
    for(int ie=0;ie<aEdge.size();++ie){
      aEdge[ie].MovePoints(aVertex);
    }
    for(int ifc=0;ifc<aFace.size();++ifc){
      aFace[ifc].MovePoints(aVertex,aEdge);
    }
  }
  return is_moved;
}


void DivideFace
(int ifc,
 const CVector3& org, int inorm,
 std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen)
{
  if( inorm == -1 ) return;
//  const CVector3& plane_ex = CVector3::Axis((inorm+1)%3);
//  const CVector3& plane_ey = CVector3::Axis((inorm+2)%3);
  const CVector3 n01 = CVector3::Axis(inorm);
  const std::vector< std::pair<int,bool> > aIE = aFace[ifc].aIE;
  const int nie = (int)aIE.size();
  std::set<int> setIV_new;
  for(int iie=0;iie<nie;++iie){
    const int ie0 = aIE[iie].first;
    double ratio;
    if (!aEdge[ie0].GetParameterIntersection(ratio, org, n01)) continue;
    int iv0 = -1;
    if(      fabs(ratio-1.0)<1.0e-5 ){ iv0 = aEdge[ie0].iv1; }
    else if( fabs(ratio-0.0)<1.0e-5 ){ iv0 = aEdge[ie0].iv0; }
    else if( ratio < 0.01 || ratio > 0.99 ) continue;
    else{ iv0 = AddPointEdge(ie0, ratio, aVertex, aEdge, aFace, elen); }
    setIV_new.insert(iv0);
  }
  if( setIV_new.size() != 2 ) return;
  int iv0 = *(setIV_new.begin());
  int iv1 = *(++setIV_new.begin());
  {
    CVector3 p0 = aVertex[iv0].pos;
    CVector3 p1 = aVertex[iv1].pos;
    CVector3 n0 = aVertex[iv0].norm;
    CVector3 n1 = aVertex[iv1].norm;
    CVector3 v = Cross(n0+n1,n01);
    if( v*(p1-p0) > 0 ){ ConectEdge(iv0,iv1, ifc, inorm, aVertex, aEdge, aFace, elen); }
    else{                ConectEdge(iv1,iv0, ifc, inorm, aVertex, aEdge, aFace, elen); }
  }
}

void BuildTriMesh
(std::vector<double>& aXYZ,
 std::vector<unsigned int>& aTri,
 std::vector<int>& aTriSurRel,
 std::vector<double>& aNorm,
 std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 int isym)
{
  std::vector<int> aIsSymVtx(aVertex.size(),0);
  for(int ie=0;ie<aEdge.size();++ie){
    const CCad3D_Edge& e = aEdge[ie];
    if( e.is_sim ){
      assert( e.inorm == isym );
      int iv0 = e.iv0;
      int iv1 = e.iv1;
      aIsSymVtx[iv0] = 1;
      aIsSymVtx[iv1] = 1;
    }
  }
  aXYZ.resize(0);
  for(int iv=0;iv<aVertex.size();++iv){
    CCad3D_Vertex& v = aVertex[iv];
    int iq0 = (int)aXYZ.size()/3;
    aXYZ.push_back(+v.pos.x);
    aXYZ.push_back(+v.pos.y);
    aXYZ.push_back(+v.pos.z);
    v.iq_right = iq0;
    if( aIsSymVtx[iv] == 1 ){
      aVertex[iv].iq_left = iq0;
    }
    else{
      if( isym == 2 ){
        aXYZ.push_back(+v.pos.x);
        aXYZ.push_back(+v.pos.y);
        aXYZ.push_back(-v.pos.z);
      }
      else if( isym == 0 ){
        aXYZ.push_back(-v.pos.x);
        aXYZ.push_back(+v.pos.y);
        aXYZ.push_back(+v.pos.z);
      }
      aVertex[iv].iq_left = iq0+1;
    }
  }
  for(int ie=0;ie<aEdge.size();++ie){
    CCad3D_Edge& e = aEdge[ie];
    const int iv0 = e.iv0;
    const int iv1 = e.iv1;
    assert( iv0>=0 && iv0<aVertex.size() );
    assert( iv1>=0 && iv1<aVertex.size() );
    const int np = (int)e.aP.size();
    assert(np>=2);
    e.aIQ_RightLeft.resize(np*2);
    e.aIQ_RightLeft[0*2+0] = aVertex[iv0].iq_right;
    e.aIQ_RightLeft[0*2+1] = aVertex[iv0].iq_left;
    e.aIQ_RightLeft[(np-1)*2+0] = aVertex[iv1].iq_right;
    e.aIQ_RightLeft[(np-1)*2+1] = aVertex[iv1].iq_left;
    for(int ip=1;ip<np-1;++ip){
      int iq0 = (int)aXYZ.size()/3;
      aXYZ.push_back(+e.aP[ip].x);
      aXYZ.push_back(+e.aP[ip].y);
      aXYZ.push_back(+e.aP[ip].z);
      e.aIQ_RightLeft[ip*2+0] = iq0;
      if( e.is_sim ){
        e.aIQ_RightLeft[ip*2+1] = iq0;
      }
      else{
        if( isym == 2 ){
          aXYZ.push_back(+e.aP[ip].x);
          aXYZ.push_back(+e.aP[ip].y);
          aXYZ.push_back(-e.aP[ip].z);
        }
        else if( isym == 0 ){
          aXYZ.push_back(-e.aP[ip].x);
          aXYZ.push_back(+e.aP[ip].y);
          aXYZ.push_back(+e.aP[ip].z);
        }
        e.aIQ_RightLeft[ip*2+1] = iq0+1;
      }
    }
  }
  for(int ifc=0;ifc<aFace.size();++ifc){
    CCad3D_Face& fc = aFace[ifc];
    int np = (int)fc.aPInfo.size();
    for(int ip=0;ip<np;++ip){
      CCad3D_Face::CFacePointInfo& pinfo = fc.aPInfo[ip];
      if( pinfo.itype == 0 ){
        int iv0 = pinfo.iv;
        pinfo.iq_right = aVertex[iv0].iq_right;
        pinfo.iq_left  = aVertex[iv0].iq_left;
      }
      else if( pinfo.itype == 1 ){
        int ie0 = pinfo.ie;
        int iep0 = pinfo.iep;
        pinfo.iq_right = aEdge[ie0].aIQ_RightLeft[iep0*2+0];
        pinfo.iq_left  = aEdge[ie0].aIQ_RightLeft[iep0*2+1];
      }
      else if( pinfo.itype == 2 ){
        int iq0 = (int)aXYZ.size()/3;
        aXYZ.push_back(+fc.aXYZ[ip*3+0]);
        aXYZ.push_back(+fc.aXYZ[ip*3+1]);
        aXYZ.push_back(+fc.aXYZ[ip*3+2]);
        ////
        if( isym == 2 ){
          aXYZ.push_back(+fc.aXYZ[ip*3+0]);
          aXYZ.push_back(+fc.aXYZ[ip*3+1]);
          aXYZ.push_back(-fc.aXYZ[ip*3+2]);
        }
        else if( isym == 0 ){
          aXYZ.push_back(-fc.aXYZ[ip*3+0]);
          aXYZ.push_back(+fc.aXYZ[ip*3+1]);
          aXYZ.push_back(+fc.aXYZ[ip*3+2]);
        }
        pinfo.iq_right = iq0;
        pinfo.iq_left  = iq0+1;
      }
    }
  }
  aTri.resize(0);
  for(int ifc=0;ifc<aFace.size();++ifc){
    CCad3D_Face& fc = aFace[ifc];
    for(int it=0;it<fc.aTri.size()/3;++it){
      int ip0 = fc.aTri[it*3+0];
      int ip1 = fc.aTri[it*3+1];
      int ip2 = fc.aTri[it*3+2];
      aTri.push_back(fc.aPInfo[ip0].iq_right);
      aTri.push_back(fc.aPInfo[ip1].iq_right);
      aTri.push_back(fc.aPInfo[ip2].iq_right);
      ///
      aTri.push_back(fc.aPInfo[ip0].iq_left);
      aTri.push_back(fc.aPInfo[ip2].iq_left);
      aTri.push_back(fc.aPInfo[ip1].iq_left);
    }
  }
  for(int ie=0;ie<aEdge.size();++ie){
    const CCad3D_Edge& e = aEdge[ie];
    if( !e.is_sim ){ continue; }
    for(int ip=0;ip<e.aP.size();++ip){
      int iq0 = e.aIQ_RightLeft[ip*2+0];
      aXYZ[iq0*3+2] += (double)rand()/(RAND_MAX+1.0)*1.0e-5;
    }
  }
  makeSurroundingRelationship(aTriSurRel,
                              aTri.data(),aTri.size()/3,
                              MESHELEM_TRI,
                              (int)aXYZ.size()/3);
  aNorm.resize(aXYZ.size());
  Normal_MeshTri3D(aNorm.data(),
                   aXYZ.data(), aXYZ.size()/3,
                   aTri.data(), aTri.size()/3);

}

void UpdateTriMesh
(std::vector<double>& aXYZ, std::vector<double>& aNorm,
 /////
 const std::vector<unsigned int>& aTri,
 const std::vector<CCad3D_Vertex>& aVertex,
 const std::vector<CCad3D_Edge>& aEdge,
 const std::vector<CCad3D_Face>& aFace,
 int isym)
{
  for(int ifc=0;ifc<aFace.size();++ifc){
    const CCad3D_Face& fc = aFace[ifc];
    for(int ip=0;ip<fc.aPInfo.size();++ip){
      int iq0 = fc.aPInfo[ip].iq_right;
      int iq1 = fc.aPInfo[ip].iq_left;
      assert( iq0 < (int)aXYZ.size()/3 );
      assert( iq1 < (int)aXYZ.size()/3 );
      aXYZ[iq0*3+0] = +fc.aXYZ[ip*3+0];
      aXYZ[iq0*3+1] = +fc.aXYZ[ip*3+1];
      aXYZ[iq0*3+2] = +fc.aXYZ[ip*3+2];
      if( isym == 0 ){
        aXYZ[iq1*3+0] = -fc.aXYZ[ip*3+0];
        aXYZ[iq1*3+1] = +fc.aXYZ[ip*3+1];
        aXYZ[iq1*3+2] = +fc.aXYZ[ip*3+2];
      }
      else{
        aXYZ[iq1*3+0] = +fc.aXYZ[ip*3+0];
        aXYZ[iq1*3+1] = +fc.aXYZ[ip*3+1];
        aXYZ[iq1*3+2] = -fc.aXYZ[ip*3+2];
      }
    }
  }
  for(int ie=0;ie<aEdge.size();++ie){
    const CCad3D_Edge& e = aEdge[ie];
    if( !e.is_sim ){ continue; }
    for(int ip=0;ip<e.aP.size();++ip){
      int iq0 = e.aIQ_RightLeft[ip*2+0];
      aXYZ[iq0*3+2] += (double)rand()/(RAND_MAX+1.0)*1.0e-5;
    }
  }
  aNorm.resize(aXYZ.size());
  Normal_MeshTri3D(aNorm.data(),
                   aXYZ.data(), aXYZ.size()/3,
                   aTri.data(), aTri.size()/3);

}


void CCad3D::Pick
(const CVector3& src_pick, const CVector3& dir_pick,
 CVector2 sp0, float mMV[16], float mPj[16],
 double view_height)
{
  if (ivtx_picked>=0&&ivtx_picked<(int)aVertex.size()){
    ielem_vtx_picked = 0;
    for (int iaxis = 0; iaxis<3; iaxis++){
      if( aVertex[ivtx_picked].isConst[iaxis] ) continue;
      CVector3 axis(0, 0, 0);
      axis[iaxis] = 1;
      if (isPick_AxisHandler(sp0, aVertex[ivtx_picked].pos,
                             axis, 0.5,
                             mMV, mPj, 0.03)){
        ielem_vtx_picked = iaxis+1;
        break;
      }
    }
  }
  if( ielem_vtx_picked == 0 || (ivtx_picked<0||ivtx_picked>=(int)aVertex.size()) ){
    ivtx_picked = -1;
    for(int icp=0;icp<aVertex.size();++icp){
      const CVector3& pos = aVertex[icp].pos;
      CVector3 pn = nearest_Line_Point(pos, src_pick, dir_pick);
      if( (pn-pos).Length() < 0.05 ){
        ivtx_picked = icp;
        ielem_vtx_picked = 0;
        break;
      }
    }
  }
  if( ivtx_picked>=0&&ivtx_picked<(int)aVertex.size() ){
    iedge_picked = -1;
    plane_inorm = -1;
    iface_picked = -1;
    aIE_picked.clear();
    return;
  }
  
  // edge pick
  if( iedge_picked != -1 ){ // edge was picked
    ielem_edge_picked = 0;
    {
      CVector2 sp = screenXYProjection(aEdge[iedge_picked].q0, mMV, mPj);
      if( (sp0-sp).Length() < 0.05 ){
        ielem_edge_picked = 2;
        iface_picked = -1;
        return;
      }
    }
    {
      CVector2 sp = screenXYProjection(aEdge[iedge_picked].q1, mMV, mPj);
      if( (sp0-sp).Length() < 0.05 ){
        ielem_edge_picked = 3;
        iface_picked = -1;
        return;
      }
    }
    if( plane_inorm>=0 && plane_inorm<3 ){ // plane pick
      CVector3 plane_ex = CVector3::Axis((plane_inorm+1)%3);
      CVector3 plane_ey = CVector3::Axis((plane_inorm+2)%3);
      CVector3 aP[4] = {
        plane_org-plane_sizeX*plane_ex-plane_sizeY*plane_ey,
        plane_org+plane_sizeX*plane_ex-plane_sizeY*plane_ey,
        plane_org+plane_sizeX*plane_ex+plane_sizeY*plane_ey,
        plane_org-plane_sizeX*plane_ex+plane_sizeY*plane_ey };
      CVector2 sp[4]  = {
        screenXYProjection(aP[0], mMV, mPj),
        screenXYProjection(aP[1], mMV, mPj),
        screenXYProjection(aP[2], mMV, mPj),
        screenXYProjection(aP[3], mMV, mPj) };
      double d01 = GetDist_LineSeg_Point(sp0, sp[0],sp[1]);
      double d12 = GetDist_LineSeg_Point(sp0, sp[1],sp[2]);
      double d23 = GetDist_LineSeg_Point(sp0, sp[2],sp[3]);
      double d30 = GetDist_LineSeg_Point(sp0, sp[3],sp[0]);
      if( d01 < 0.05 || d12 < 0.05 || d23 < 0.05 || d30 < 0.05 ) {
        ielem_edge_picked = 1;
        iface_picked = -1;
        return;
      }
    }
  }
  /////
  plane_inorm = -1;
  iedge_picked = -1;
  aIE_picked.clear();
  {
    std::map<double, std::pair<int, double> > mapDepthEdge;
    for(int ie=0;ie<aEdge.size();++ie){
      double ratio_edge;
      bool res = aEdge[ie].isPick(ratio_edge, sp0, mMV, mPj);
      if( res ){
        CVector3 p = aEdge[ie].GetPosInEdge(ratio_edge);
        double depth = -p*dir_pick;
        mapDepthEdge.insert( std::make_pair(depth, std::make_pair(ie,ratio_edge) ) );
      }
    }
    if( !mapDepthEdge.empty() ){
      iedge_picked = mapDepthEdge.begin()->second.first;
      ratio_edge_picked = mapDepthEdge.begin()->second.second;
      ielem_vtx_picked = 0;
      iface_picked = -1;
    }
  }
  if( iedge_picked != -1 ){ // make plane
    findEdgeGroup(aIE_picked, iedge_picked, aVertex, aEdge);
    plane_inorm = aEdge[iedge_picked].inorm;
    const CVector3 n = CVector3::Axis(plane_inorm);
    const CVector3 plane_ex = CVector3::Axis((plane_inorm+1)%3);
    const CVector3 plane_ey = CVector3::Axis((plane_inorm+2)%3);
    int iv0 = aEdge[iedge_picked].iv0;
    int iv1 = aEdge[iedge_picked].iv1;
    plane_org = (aVertex[iv0].pos+aVertex[iv1].pos)*0.5;
    double minX=0, maxX=0, minY=0, maxY=0;
    for(int iie=0;iie<aIE_picked.size();++iie){
      int ie = aIE_picked[iie].first;
      std::vector<CVector3>& aP = aEdge[ie].aP;
      for(int ip=0;ip<aP.size();++ip){
        const CVector3& p = aP[ip];
        double x0 = (p-plane_org)*plane_ex;
        double y0 = (p-plane_org)*plane_ey;
        if( iie==0 && ip==0 ){
          minX = x0;
          minY = y0;
          maxX = x0;
          maxY = y0;
          continue;
        }
        minX = (x0<minX) ? x0 : minX;
        maxX = (x0>maxX) ? x0 : maxX;
        minY = (y0<minY) ? y0 : minY;
        maxY = (y0>maxY) ? y0 : maxY;
      }
    }
    plane_org += ((minX+maxX)*plane_ex + (minY+maxY)*plane_ey)*0.5;
    plane_sizeX = 0.5*(maxX-minX) + view_height*0.3;
    plane_sizeY = 0.5*(maxY-minY) + view_height*0.3;
    return;
  }
  
  iface_picked = -1;
  for(int ifc=0;ifc<aFace.size();++ifc){
    if( aFace[ifc].isPick(src_pick, dir_pick) ){
      iface_picked = ifc;
      return;
    }
  }
}

void CCad3D::DrawFace_LeftRight() const
{
  ::glEnable(GL_LIGHTING);
  ::myGlColorDiffuse(color_face);
  DrawMeshTri3D_FaceNorm(aXYZ, aTri, aNorm);
}

void CCad3D::DrawFace_RightSelected(bool is_edge) const
{
  {
    float specular[4] ={ 0.797357f, 0.723991f, 0.208006f, 1.0f};
    float shine =83.2f ;
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT, GL_SHININESS, shine);
  }
  
  ::glEnable(GL_LIGHTING);
  for(int iface=0;iface<aFace.size();++iface){
    if( iface == iface_picked ){
      ::myGlColorDiffuse(color_face_selected);
    }
    else{
      ::myGlColorDiffuse(color_face);
    }
    aFace[iface].DrawFace();
    if( is_edge ){ aFace[iface].DrawEdge(); }
  }
}

void CCad3D::DrawVtxEdgeHandler(double view_height) const
{

  for(int icp=0;icp<aVertex.size();++icp){
    aVertex[icp].Draw(icp==ivtx_picked,ielem_vtx_picked, view_height);
  }
  
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
//  glEnable(GL_MULTISAMPLE);
  
  myGlColorDiffuse(CColor::Blue());
  for(int ie=0;ie<aEdge.size();++ie){
    bool is_loop0 = std::find(aIE_picked.begin(), aIE_picked.end(), std::make_pair(ie,true ) ) != aIE_picked.end();
    bool is_loop1 = std::find(aIE_picked.begin(), aIE_picked.end(), std::make_pair(ie,false) ) != aIE_picked.end();
    bool is_loop = is_loop0 || is_loop1;
    aEdge[ie].DrawLine(is_loop, view_height);
  }
  if( iedge_picked != -1 ){
    aEdge[iedge_picked].DrawHandler(ielem_edge_picked, view_height);
  }

  ::glDisable(GL_LIGHTING);
  ::glDisable(GL_CULL_FACE);
  if( plane_inorm >= 0 && plane_inorm < 3 ){
    myGlColorDiffuse(CColor(0.0, 0.0, 1.0, 0.5));
    CVector3 plane_ex = CVector3::Axis((plane_inorm+1)%3);
    CVector3 plane_ey = CVector3::Axis((plane_inorm+2)%3);
    CVector3 p0 = plane_org-plane_sizeX*plane_ex-plane_sizeY*plane_ey;
    CVector3 p1 = plane_org+plane_sizeX*plane_ex-plane_sizeY*plane_ey;
    CVector3 p2 = plane_org+plane_sizeX*plane_ex+plane_sizeY*plane_ey;
    CVector3 p3 = plane_org-plane_sizeX*plane_ex+plane_sizeY*plane_ey;
    ::glBegin(GL_QUADS);
    ::myGlVertex(p0);
    ::myGlVertex(p1);
    ::myGlVertex(p2);
    ::myGlVertex(p3);
    ::glEnd();
    ////
    ::glLineWidth(5);
    if( ielem_edge_picked == 1 ){ myGlColorDiffuse(CColor(1.0, 0.0, 0.0, 0.9)); }
    else{                         myGlColorDiffuse(CColor(0.0, 0.0, 1.0, 0.9)); }
    DrawCylinder(p0, p1, view_height*0.01);
    DrawCylinder(p1, p2, view_height*0.01);
    DrawCylinder(p2, p3, view_height*0.01);
    DrawCylinder(p3, p0, view_height*0.01);
  }
  
  {
    ::glDisable(GL_DEPTH_TEST);
    ::glMatrixMode(GL_PROJECTION);
    ::glPushMatrix();
    ::glLoadIdentity();
    ::glMatrixMode(GL_MODELVIEW);
    ::glPushMatrix();
    ::glLoadIdentity();
    //////
    ::glDisable(GL_LIGHTING);
    ::glColor3d(1,0,0);
    ::glLineWidth(3);
    ::glBegin(GL_LINE_STRIP);
    for(int ist=0;ist<aStroke.size();++ist){
      ::myGlVertex(aStroke[ist]);
    }
    ::glEnd();
    //////
    ::glMatrixMode(GL_PROJECTION);
    ::glPopMatrix();
    ::glMatrixMode(GL_MODELVIEW);
    ::glPopMatrix();
    ::glEnable(GL_DEPTH_TEST);
  }
  
  glDisable(GL_LINE_SMOOTH);
  
}

void CCad3D::MouseUp(float mMV[16], float mPj[16], double view_height)
{
  
  if( imode_edit == EDIT_SKETCH ){
    if( aStroke.size() > 3 && iedge_picked != -1 ){
      bool res = MovePointsAlongSketch(aVertex,aEdge,aFace,
                                  aStroke,aIE_picked,
                                  plane_org,aEdge[iedge_picked].inorm,
                                  mMV,mPj,view_height);
      if( res ){
        UpdateTriMesh(aXYZ,aNorm, aTri,aVertex,aEdge,aFace,isym);
      }
      else{
        plane_inorm = -1;
        iedge_picked = -1;
        aIE_picked.clear();
      }
    }
    aStroke.clear();
  }
}

bool CCad3D::ReflectChangeForCurveAndSurface
(std::vector<int>& aIsMoved_Edge,
 const std::vector<int>& aIsMoved_Vtx)
{
  bool is_edit = false;
  for(int ie=0;ie<aEdge.size();++ie){
    int iv0 = aEdge[ie].iv0;
    int iv1 = aEdge[ie].iv1;
    if( aIsMoved_Vtx[iv0]==0 && aIsMoved_Vtx[iv1]==0 && aIsMoved_Edge[ie]==0 ) continue;
    aIsMoved_Edge[ie] = 1;
    aEdge[ie].MovePoints(aVertex);
    is_edit = true;
  }
  for(int ifc=0;ifc<aFace.size();++ifc){
    bool is_edit_face = false;
    for(int iie=0;iie<aFace[ifc].aIE.size();++iie){
      int ie0 = aFace[ifc].aIE[iie].first;
      if( aIsMoved_Edge[ie0] == 0 ) continue;
      is_edit_face = true;
      break;
    }
    if( !is_edit_face ) continue;
    aFace[ifc].MovePoints(aVertex,aEdge);
    is_edit = true;
  }
  if( is_edit ){
    UpdateTriMesh(aXYZ,aNorm, aTri,aVertex,aEdge,aFace,isym);
  }
  return is_edit;
}

bool CCad3D::MouseMotion
(const CVector3& src_pick, const CVector3& dir_pick,
 const CVector2& sp0, const CVector2& sp1,
 float mMV[16], float mPj[16])
{
  if( imode_edit == EDIT_MOVE ){
    std::vector<int> aIsMoved_Vtx(aVertex.size(),0);
    std::vector<int> aIsMoved_Edge(aEdge.size(),0);
    if( ivtx_picked>=0 && ivtx_picked < aVertex.size() ){ // move vtx
      if( ielem_vtx_picked <= 0 || ielem_vtx_picked > 3 ){ return false; }
      int iaxis = ielem_vtx_picked-1;
      CVector3 axis(0, 0, 0); axis[iaxis] = 1;
      CVector3 d0 = drag_AxisHandler(sp0, sp1, aVertex[ivtx_picked].pos, axis, 0.5, mMV, mPj);
      aVertex[ivtx_picked].pos += d0;
      aIsMoved_Vtx[ivtx_picked] = 1;
    }
    if( iedge_picked>=0 && iedge_picked<aEdge.size() && (ielem_edge_picked==2 || ielem_edge_picked==3) ){ // moved edge ctrl point
      const int ie0 = iedge_picked;
      CVector3 axis = CVector3::Axis(aEdge[ie0].inorm);
      CVector3 qe = intersection_Plane_Line(plane_org, axis, src_pick, dir_pick);
      const int iv = (ielem_edge_picked==2) ? aEdge[ie0].iv0 : aEdge[ie0].iv1;
      CVector3  qs = (ielem_edge_picked==2) ? aEdge[ie0].q0  : aEdge[ie0].q1;
      CVector3  pv = (ielem_edge_picked==2) ? aEdge[ie0].p0  : aEdge[ie0].p1;
      const double len0 = (aEdge[ie0].p0-aEdge[ie0].p1).Length();
      if( (qs-pv).Length() < 0.01*len0 ) return false;
      if( (qe-pv).Length() < 0.01*len0 ) return false;
      CMatrix3 R = Mat3_MinimumRotation(qs-pv, qe-pv);
      aVertex[iv].norm = R*(aVertex[iv].norm);
      double len1 = (qe-pv).Length();
      if( ielem_edge_picked==2 ){ aEdge[ie0].r0 = len1/len0; }
      else{                       aEdge[ie0].r1 = len1/len0; }
      aIsMoved_Vtx[iv] = 1;
      aIsMoved_Edge[iedge_picked] = 1;
    }
    if( iedge_picked>=0 && iedge_picked<aEdge.size() && ielem_edge_picked == 1 ){ // move vtx on the plane
      if( !aEdge[iedge_picked].is_sim  ){
        int iaxis = aEdge[iedge_picked].inorm;
        CVector3 axis(0, 0, 0); axis[iaxis] = 1;
        CVector2 spa0 = screenXYProjection(plane_org+axis, mMV, mPj);
        CVector2 spa1 = screenXYProjection(plane_org-axis, mMV, mPj);
        double r = (spa0-spa1)*(sp1-sp0)/(spa0-spa1).SqLength();
        CVector3 d = r*axis;
        plane_org += d;
        std::vector<int> aIP = getPointsInEdges(aIE_picked, aEdge);
        for(int iip=0;iip<aIP.size();++iip){
          int ip0 = aIP[iip];
          aVertex[ip0].pos += d;
          aIsMoved_Vtx[ip0] = 1;
        }
      }
    }
    return ReflectChangeForCurveAndSurface(aIsMoved_Edge,aIsMoved_Vtx);
  }
  else if( imode_edit == EDIT_SKETCH ){
    aStroke.push_back(sp0);
    return false;
  }
  else if( imode_edit == EDIT_ADD_CROSS_SECTION ){
    if( plane_inorm >= 0 && plane_inorm < 3 ){
      int iaxis = aEdge[iedge_picked].inorm;
      if( iaxis>=0 && iaxis<3 ){
        CVector3 axis(0, 0, 0); axis[iaxis] = 1;
        CVector2 spa0 = screenXYProjection(plane_org+axis, mMV, mPj);
        CVector2 spa1 = screenXYProjection(plane_org-axis, mMV, mPj);
        double r = (spa0-spa1)*(sp1-sp0)/(spa0-spa1).SqLength();
        CVector3 d = r*axis;
        plane_org += d;
        return false;
      }
    }
  }
  return false;
}

void CCad3D::MouseDown
(const CVector3& src_pick, const CVector3& dir_pick,
 const CVector2& sp0, float mMV[16], float mPj[16],
 double view_height)
{
  if( imode_edit == EDIT_ADD_POINT_EDGE ){
    Pick(src_pick, dir_pick, sp0, mMV, mPj, view_height);
    AddPointEdge(iedge_picked, ratio_edge_picked, aVertex, aEdge, aFace, elen);
    BuildTriMesh(aXYZ,aTri,aTriSurRel,aNorm, aVertex,aEdge,aFace, isym);
    this->iedge_picked = -1;
    this->aIE_picked.clear();
    this->plane_inorm = -1;
    imode_edit = EDIT_NONE;
  }
  else if( imode_edit == EDIT_SKETCH ){
    if( iedge_picked == -1 ){
      Pick(src_pick, dir_pick, sp0, mMV, mPj, view_height);
    }
  }
  else if( imode_edit == EDIT_ADD_CROSS_SECTION ){
    if( plane_inorm != -1 ){
      for(int ifc=0;ifc<aFace.size();++ifc){
        if( aFace[ifc].isPick(src_pick, dir_pick) ){
          DivideFace(ifc,plane_org,plane_inorm,
                     aVertex,aEdge,aFace, elen);
          BuildTriMesh(aXYZ,aTri,aTriSurRel,aNorm, aVertex,aEdge,aFace, isym);
          return;
        }
      }
    }
    Pick(src_pick, dir_pick, sp0, mMV, mPj, view_height);
  }
  else if( imode_edit == EDIT_MOVE ){
    Pick(src_pick, dir_pick, sp0, mMV, mPj, view_height);
  }
}

bool CCad3D::isSym(int iv) const{
  for(int ie=0;ie<aEdge.size();++ie){
    if( !aEdge[ie].is_sim ) continue;
    if( aEdge[ie].iv0 == iv ) return true;
    if( aEdge[ie].iv1 == iv ) return true;
  }
  return false;
}

void CCad3D::WriteFile(std::ofstream& fout) const
{
  fout << aVertex.size() << std::endl;
  for(int iv=0;iv<aVertex.size();++iv){
    fout << " " << iv << std::endl;
    aVertex[iv].WriteFile(fout);
  }
  ////
  fout << aEdge.size() << std::endl;
  for(int ie=0;ie<aEdge.size();++ie){
    fout << " " << ie << std::endl;
    aEdge[ie].WriteFile(fout);
  }
  ///
  fout << aFace.size() << std::endl;
  for(int ifc=0;ifc<aFace.size();++ifc){
    fout << " " << ifc << std::endl;
    aFace[ifc].WriteFile(fout);
  }
}

void CCad3D::ReadFile(std::ifstream& fin)
{
  int nv = 0;
  fin >> nv;
  aVertex.resize(nv);
  for(int iv=0;iv<aVertex.size();++iv){
    int iv0;
    fin >> iv0;
    assert( iv0 == iv );
    aVertex[iv].ReadFile(fin);
  }
  ////
  int ne;
  fin >> ne;
  aEdge.resize(ne);
  for(int ie=0;ie<aEdge.size();++ie){
    int ie0;
    fin >> ie0;
    assert( ie0 == ie );
    aEdge[ie].ReadFile(fin);
  }
  ////
  int nfc;
  fin >> nfc;
  aFace.resize(nfc);
  for(int ifc=0;ifc<aFace.size();++ifc){
    int ifc0;
    fin >> ifc0;
    assert( ifc0 == ifc );
    aFace[ifc].ReadFile(fin);
  }
  {
    for(int iv=0;iv<aVertex.size();++iv){
      aVertex[iv].isConst[0] = false;
      aVertex[iv].isConst[1] = false;
      aVertex[iv].isConst[2] = false;
    }
    for(int ie=0;ie<aEdge.size();++ie){
      int iv0 = aEdge[ie].iv0;
      int iv1 = aEdge[ie].iv1;
      int inorm = aEdge[ie].inorm;
      if( inorm < 0 || inorm >= 3 ){ continue; }
      aVertex[iv0].isConst[inorm] = true;
      aVertex[iv1].isConst[inorm] = true;
    }
  }
  for(int ie=0;ie<aEdge.size();ie++){
    aEdge[ie].Initialize(aVertex,elen); // ie0+0
  }
  for(int ifc=0;ifc<aFace.size();ifc++){
    aFace[ifc].Initialize(aVertex,aEdge, elen); // ie0+0
  }
  BuildTriMesh(aXYZ,aTri,aTriSurRel,aNorm, aVertex,aEdge,aFace, isym);
}







