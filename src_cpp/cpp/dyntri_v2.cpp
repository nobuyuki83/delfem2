#include <set>
#include <algorithm>

#include "delfem2/dyntri_v2.h"


/*
static inline double TriArea2D(const double v1[], const double v2[], const double v3[]){
  return 0.5*( (v2[0]-v1[0])*(v3[1]-v1[1]) - (v3[0]-v1[0])*(v2[1]-v1[1]) );
}
 */

static bool IsCrossLines(const double po_s0[], const double po_e0[],
                         const double po_s1[], const double po_e1[] )
{
  const double area1 = TriArea2D(po_s0,po_e0,po_s1);
  const double area2 = TriArea2D(po_s0,po_e0,po_e1);
  if( area1 * area2 > 0.0 ) return false;
  const double area3 = TriArea2D(po_s1,po_e1,po_s0);
  const double area4 = TriArea2D(po_s1,po_e1,po_e0);
  if( area3 * area4 > 0.0 ) return false;
  return true;
}

/*
double DistanceXY(const CVector3& p0, const CVector3& p1 )
{
  const double dx = p0.x-p1.x;
  const double dy = p0.y-p1.y;
  return sqrt(dx*dx+dy*dy);
}

double TriAreaXY(const CVector3& v1, const CVector3& v2, const CVector3& v3){
  return 0.5*( (v2.x-v1.x)*(v3.y-v1.y) - (v3.x-v1.x)*(v2.y-v1.y) );
}

double SquareDistanceXY(const CVector3& ipo0, const CVector3& ipo1)
{
  return	( ipo1.x - ipo0.x )*( ipo1.x - ipo0.x ) + ( ipo1.y - ipo0.y )*( ipo1.y - ipo0.y );
}

static inline int DetDelaunayXY
(const CVector3& p0,
 const CVector3& p1,
 const CVector3& p2,
 const CVector3& p3)
{
  const double area = TriAreaXY(p0,p1,p2);
  if( fabs(area) < 1.0e-10 ){
    return 3;
  }
  const double tmp_val = 1.0/(area*area*16.0);
  
  const double dtmp0 = SquareDistanceXY(p1,p2);
  const double dtmp1 = SquareDistanceXY(p0,p2);
  const double dtmp2 = SquareDistanceXY(p0,p1);
  
  const double etmp0 = tmp_val*dtmp0*(dtmp1+dtmp2-dtmp0);
  const double etmp1 = tmp_val*dtmp1*(dtmp0+dtmp2-dtmp1);
  const double etmp2 = tmp_val*dtmp2*(dtmp0+dtmp1-dtmp2);
  
  const CVector3 out_center(etmp0*p0.x + etmp1*p1.x + etmp2*p2.x,
                            etmp0*p0.y + etmp1*p1.y + etmp2*p2.y,
                            0.0);
  
  const double qradius = SquareDistanceXY(out_center,p0);
  const double qdistance = SquareDistanceXY(out_center,p3);
  
  //	assert( fabs( qradius - SquareLength(out_center,p1) ) < 1.0e-10*qradius );
  //	assert( fabs( qradius - SquareLength(out_center,p2) ) < 1.0e-10*qradius );
  
  const double tol = 1.0e-20;
  if( qdistance > qradius*(1.0+tol) ){ return 2; }
  else{
    if( qdistance < qradius*(1.0-tol) ){ return 0; }
    else{ return 1;	}
  }
  return 0;
}
 */

//////////////////////////////////////////////////////////////////////////////////

bool FindEdgePoint_AcrossEdge2
(int& itri0, int& inotri0, int& inotri1, double& ratio,
 const int& ipo0, const int& ipo1,
 std::vector<CEPo2>& po,
 std::vector<ETri>& tri,
 std::vector<CVector2>& aVec2)
{
  const unsigned int itri_ini = po[ipo0].e;
  const unsigned int inotri_ini = po[ipo0].d;
  unsigned int inotri_cur = inotri_ini;
  unsigned int itri_cur = itri_ini;
  for(;;){
    assert( tri[itri_cur].v[inotri_cur] == ipo0 );
    {
      const unsigned int inotri2 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
      const unsigned int inotri3 = (inotri_cur+2)%3;//  indexRot3[2][inotri_cur];
      double area0 = TriArea(aVec2[ipo0],
                             aVec2[ tri[itri_cur].v[inotri2] ],
                             aVec2[ipo1] );
      if( area0 > -1.0e-20 ){
        double area1 =  TriArea(aVec2[ipo0],
                                aVec2[ipo1],
                                aVec2[ tri[itri_cur].v[inotri3] ] );
        if( area1 > -1.0e-20 ){
          assert( area0 + area1 > 1.0e-20 );
          ratio = area0 / ( area0 + area1 );
          itri0 = itri_cur;
          inotri0 = inotri2;
          inotri1 = inotri3;
          return true;
        }
      }
    }
    {
      const unsigned int inotri2 = (inotri_cur+1)%3;// indexRot3[1][inotri_cur];
      //      if( tri[itri_cur].g2[inotri2] != -2 && tri[itri_cur].g2[inotri2] != -3 ){
      //        break;
      //      }
      unsigned int itri_nex = tri[itri_cur].s2[inotri2];
      const unsigned int* rel = relTriTri[ tri[itri_cur].r2[inotri2] ];
      const unsigned int inotri3 = rel[inotri_cur];
      assert( tri[itri_nex].v[inotri3] == ipo0 );
      if( itri_nex == itri_ini ){
        itri0 = 0;
        inotri0 = 0; inotri1 = 0;
        ratio = 0.0;
        return false;
      }
      itri_cur = itri_nex;
      inotri_cur = inotri3;
    }
  }
  
  inotri_cur = inotri_ini;
  itri_cur = itri_ini;
  for(;;){
    assert( tri[itri_cur].v[inotri_cur] == ipo0 );
    {
      const unsigned int inotri2 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
      const unsigned int inotri3 = (inotri_cur+2)%3; // indexRot3[2][inotri_cur];
      double area0 = TriArea(aVec2[ipo0],
                             aVec2[ tri[itri_cur].v[inotri2] ],
                             aVec2[ipo1] );
      if( area0 > -1.0e-20 ){
        double area1 =  TriArea(aVec2[ipo0],
                                aVec2[ipo1],
                                aVec2[ tri[itri_cur].v[inotri3] ] );
        if( area1 > -1.0e-20 ){
          assert( area0 + area1 > 1.0e-20 );
          ratio = area0 / ( area0 + area1 );
          itri0 = itri_cur;
          inotri0 = inotri2;
          inotri1 = inotri3;
          return true;
        }
      }
    }
    {
      const unsigned int inotri2 = (inotri_cur+2)%3; // indexRot3[2][inotri_cur];
      //      if( tri[itri_cur].g2[inotri2] != -2 && tri[itri_cur].g2[inotri2] != -3 ){
      //        break;
      //      }
      unsigned int itri_nex = tri[itri_cur].s2[inotri2];
      const unsigned int* rel = relTriTri[ tri[itri_cur].r2[inotri2] ];
      const unsigned int inotri3 = rel[inotri_cur];
      assert( tri[itri_nex].v[inotri3] == ipo0 );
      if( itri_nex == itri_ini ){
        assert(0);	// àÍé¸ÇµÇ»Ç¢ÇÕÇ∏
      }
      itri_cur = itri_nex;
      inotri_cur = inotri3;
    }
  }
  
  itri0 = 0;
  inotri0 = 0; inotri1 = 0;
  ratio = 0.0;
  
  return false;
}

bool DelaunayAroundPoint2
(int ipo0,
 std::vector<CEPo2>& aPo,
 std::vector<ETri>& aTri,
 std::vector<CVector2>& aVec2)
{
  assert(ipo0 < (int)aPo.size());
  if (aPo[ipo0].e==-1) return true;
  
  assert(aPo[ipo0].e>=0&&(int)aPo[ipo0].e < (int)aTri.size());
  assert(aTri[aPo[ipo0].e].v[aPo[ipo0].d]==ipo0);
  
  const int itri0 = aPo[ipo0].e;
  int inotri0 = aPo[ipo0].d;
  
  int itri_cur = itri0;
  int inotri_cur = aPo[ipo0].d;
  bool flag_is_wall = false;
  for (;;){
    assert(aTri[itri_cur].v[inotri_cur]==ipo0);
    if (aTri[itri_cur].s2[inotri_cur]>=0&&aTri[itri_cur].s2[inotri_cur]<(int)aTri.size()){
      assert(aTri[itri_cur].v[inotri_cur]==ipo0);
      // check opposing element
      const int itri_dia = aTri[itri_cur].s2[inotri_cur];
      const unsigned int* rel_dia = relTriTri[aTri[itri_cur].r2[inotri_cur]];
      const int inotri_dia = rel_dia[inotri_cur];
      assert(aTri[itri_dia].s2[inotri_dia]==itri_cur);
      const int ipo_dia = aTri[itri_dia].v[inotri_dia];
      if (DetDelaunay(aVec2[aTri[itri_cur].v[0]],
                      aVec2[aTri[itri_cur].v[1]],
                      aVec2[aTri[itri_cur].v[2]],
                      aVec2[ipo_dia])==0)
      {
        bool res = FlipEdge(itri_cur, inotri_cur, aPo, aTri);
        if( res ){
          inotri_cur = 2;
          assert(aTri[itri_cur].v[inotri_cur]==ipo0);
          if (itri_cur==itri0) inotri0 = inotri_cur;
          continue;
        }
        else{
          break;
        }
      }
    }
    MoveCCW(itri_cur, inotri_cur, flag_is_wall, aPo,aTri);
    if( flag_is_wall ) break;
    if( itri_cur == itri0 ) break;
    /*
     {	// next element
     const int inotri1 = indexRot3[1][inotri_cur];
     if (aTri[itri_cur].s2[inotri1]==-1){
     flag_is_wall = true;
     break;
     }
     const int itri_nex = aTri[itri_cur].s2[inotri1];
     const unsigned int* rel_nex = relTriTri[aTri[itri_cur].r2[inotri1]];
     const int inotri_nex = rel_nex[inotri_cur];
     assert(aTri[itri_nex].v[inotri_nex]==ipo0);
     if (itri_nex==itri0) break;	// finish if we reach starting elemnt
     itri_cur = itri_nex;
     inotri_cur = inotri_nex;
     }
     */
  }
  if (!flag_is_wall) return true;
  
  ////////////////////////////////
  // rotate counter clock-wise
  
  itri_cur = itri0;
  inotri_cur = inotri0;
  for (;;){
    assert(aTri[itri_cur].v[inotri_cur]==ipo0);
    
    if (aTri[itri_cur].s2[inotri_cur]>=0&&aTri[itri_cur].s2[inotri_cur]<(int)aTri.size()){
      // check elements in opposing side
      const int itri_dia = aTri[itri_cur].s2[inotri_cur];
      const unsigned int* rel_dia = relTriTri[aTri[itri_cur].r2[inotri_cur]];
      const int inotri_dia = rel_dia[inotri_cur];
      assert(aTri[itri_dia].s2[inotri_dia]==itri_cur);
      const int ipo_dia = aTri[itri_dia].v[inotri_dia];
      if (DetDelaunay(aVec2[aTri[itri_cur].v[0]],
                      aVec2[aTri[itri_cur].v[1]],
                      aVec2[aTri[itri_cur].v[2]],
                      aVec2[ipo_dia])==0)	// Delaunay condition is not satisfiled
      {
        FlipEdge(itri_cur, inotri_cur, aPo, aTri);
        itri_cur = itri_dia;
        inotri_cur = 1;
        assert(aTri[itri_cur].v[inotri_cur]==ipo0);
        continue;
      }
    }
    
    {
      const int inotri2 = (inotri_cur+2)%3;// indexRot3[2][inotri_cur];
      if (aTri[itri_cur].s2[inotri2]==-1){
        return true;
      }
      const int itri_nex = aTri[itri_cur].s2[inotri2];
      const unsigned int* rel_nex = relTriTri[aTri[itri_cur].r2[inotri2]];
      const int inotri_nex = rel_nex[inotri_cur];
      assert(aTri[itri_nex].v[inotri_nex]==ipo0);
      assert(itri_nex!=itri0);	// finsih if reach starting elemnet
      itri_cur = itri_nex;
      inotri_cur = inotri_nex;
    }
  }
  return true;
}


void LaplacianSmoothing2
(std::vector<CVector2>& aVec2,
 const std::vector<CEPo2>& aPo,
 const std::vector<ETri>& aTri,
 const std::vector<int>& aflg_isnt_move)
{
  /*
  const unsigned int noelTriEdge[3][2] = {
    { 1, 2 }, //edge 0
    { 2, 0 }, //edge 1
    { 0, 1 }, //edge 2
  };
   */
  for(int ipoin=0;ipoin<(int)aPo.size();ipoin++){
    if( ipoin < (int)aflg_isnt_move.size() ){
      if( aflg_isnt_move[ipoin] == 1 ) continue;
    }
    const unsigned int itri_ini = aPo[ipoin].e;
    const unsigned int inoel_c_ini = aPo[ipoin].d;
    assert( itri_ini < aTri.size() );
    assert( inoel_c_ini < 3 );
    assert( aTri[itri_ini].v[inoel_c_ini] == ipoin );
    int itri0 = itri_ini;
    int inoel_c0 = inoel_c_ini;
    int inoel_b0 = (inoel_c0+1)%3;
    bool is_bound_flg = false;
    CVector2 vec_delta = aVec2[ipoin];
    unsigned int ntri_around = 1;
    for(;;){
      assert( itri0 < (int)aTri.size() );
      assert( inoel_c0 < 3 );
      assert( aTri[itri0].v[inoel_c0] == ipoin );
      {
        vec_delta.x += aVec2[ aTri[itri0].v[inoel_b0] ].x;
        vec_delta.y += aVec2[ aTri[itri0].v[inoel_b0] ].y;
        ntri_around++;
      }
      if( aTri[itri0].s2[inoel_b0] >= 0 ){
        const int itri1 = aTri[itri0].s2[inoel_b0];
        const int rel01 = (int)aTri[itri0].r2[inoel_b0];
        const int inoel_c1 = relTriTri[rel01][inoel_c0];
        unsigned int inoel_b1 = relTriTri[rel01][ (inoel_c0+2)%3 ];
        assert( itri1 < (int)aTri.size() );
        assert( aTri[itri1].s2[ relTriTri[rel01][inoel_b0] ] == itri0 );
        assert( aTri[itri1].v[inoel_c1] == ipoin );
        if( itri1 == (int)itri_ini ) break;
        itri0 = itri1;
        inoel_c0 = inoel_c1;
        inoel_b0 = inoel_b1;
      }
      else{
        is_bound_flg = true;
        break;
      }
    }
    if( is_bound_flg ) continue;
    aVec2[ipoin].x = vec_delta.x / ntri_around;
    aVec2[ipoin].y = vec_delta.y / ntri_around;
  }
}


void LaplaceDelaunaySmoothing2
(std::vector<CVector2>& aVec2,
 std::vector<CEPo2>& aPo,
 std::vector<ETri>& aTri,
 const std::vector<int>& aflg_isnt_move )
{
  for(unsigned int ipoin=0;ipoin<aPo.size();ipoin++){	// ì_é¸ÇËÇÃì_ÇíTçıÇµÇƒí≤Ç◊ÇÈÅB
    if( ipoin < aflg_isnt_move.size() ){
      if( aflg_isnt_move[ipoin] == 1 ) continue;
    }
    const unsigned int itri_ini = aPo[ipoin].e;
    const unsigned int inoel_c_ini = aPo[ipoin].d;
    assert( itri_ini < aTri.size() );
    assert( inoel_c_ini < 3 );
    assert( aTri[itri_ini].v[inoel_c_ini] == (int)ipoin );
    unsigned int itri0= itri_ini;
    unsigned int inoel_c0 = inoel_c_ini;
    unsigned int inoel_b0 = (inoel_c0+1)%3;
    bool is_bound_flg = false;
    CVector2 vec_delta = aVec2[ipoin];
    unsigned int ntri_around = 1;
    for(;;){
      assert( itri0 < aTri.size() );
      assert( inoel_c0 < 3 );
      assert( aTri[itri0].v[inoel_c0] == (int)ipoin );
      {
        vec_delta.x += aVec2[ aTri[itri0].v[inoel_b0] ].x;
        vec_delta.y += aVec2[ aTri[itri0].v[inoel_b0] ].y;
        ntri_around++;
      }
      if( aTri[itri0].s2[inoel_b0] >= 0 ){
        unsigned int itri1 = aTri[itri0].s2[inoel_b0];
        const int rel01 = (int)aTri[itri0].r2[inoel_b0];
        unsigned int inoel_c1 = relTriTri[rel01][inoel_c0];
        unsigned int inoel_b1 = relTriTri[rel01][ (inoel_c0+2)%3 ];
        assert( itri1 < aTri.size() );
        assert( aTri[itri1].s2[ relTriTri[rel01][inoel_b0] ] == (int)itri0 );
        assert( aTri[itri1].v[inoel_c1] == (int)ipoin );
        if( itri1 == itri_ini ) break;
        itri0 = itri1;
        inoel_c0 = inoel_c1;
        inoel_b0 = inoel_b1;
      }
      else{	// Ç±ÇÃì_ÇÕã´äEè„ÇÃì_ÇæÇ©ÇÁìÆÇ©ÇµÇƒÇÕÇ»ÇÁÇ»Ç¢ÅB
        is_bound_flg = true;
        break;
      }
    }
    if( is_bound_flg ) continue;
    aVec2[ipoin].x = vec_delta.x / ntri_around;
    aVec2[ipoin].y = vec_delta.y / ntri_around;
    DelaunayAroundPoint2(ipoin,aPo,aTri,aVec2);
  }
}

void MeshingInside2
(std::vector<CEPo2>& aPo2D,
 std::vector<ETri>& aTri,
 std::vector<CVector2>& aVec2,
 const std::vector<int>& aVtxInd,
 const double len,
 const CMeshDensity& mesh_density)
{
  std::vector<int> aflag_isnt_move;
  {
    aflag_isnt_move.resize( aPo2D.size(), 0 );
    for(unsigned int iver=0;iver<aVtxInd.size();iver++){
      const int ivec = aVtxInd[iver];
      aflag_isnt_move[ivec] = 1;
    }
  }
  
  {
    double ratio = 3.0;
    for(;;){
      int nadd = 0;
      for(int itri=0;itri<(int)aTri.size();itri++){
        const double area = TriArea(aVec2[aTri[itri].v[0]],
                                    aVec2[aTri[itri].v[1]],
                                    aVec2[aTri[itri].v[2]]);
        const double pcnt[2] = {
          (aVec2[aTri[itri].v[0]].x + aVec2[aTri[itri].v[1]].x + aVec2[aTri[itri].v[2]].x)/3.0,
          (aVec2[aTri[itri].v[0]].y + aVec2[aTri[itri].v[1]].y + aVec2[aTri[itri].v[2]].y)/3.0
        };
        double len2 = len*mesh_density.edgeLengthRatio(pcnt[0], pcnt[1]);
        if( area > len2 * len2 * ratio ){
          const int ipo0 = (int)aPo2D.size();
          aPo2D.resize( aPo2D.size()+1 );
          aVec2.resize( aVec2.size()+1 );
          aVec2[ipo0].x = (aVec2[aTri[itri].v[0]].x+aVec2[aTri[itri].v[1]].x+aVec2[aTri[itri].v[2]].x)/3.0;
          aVec2[ipo0].y = (aVec2[aTri[itri].v[0]].y+aVec2[aTri[itri].v[1]].y+aVec2[aTri[itri].v[2]].y)/3.0;
          InsertPoint_Elem(ipo0,itri,aPo2D,aTri);
          DelaunayAroundPoint2(ipo0,aPo2D,aTri,aVec2);
          nadd++;
        }
      }
      LaplacianSmoothing2(aVec2, aPo2D,aTri,aflag_isnt_move);
      //			LaplaceDelaunaySmoothing(aPo2D,aTri);
      if( nadd != 0 ){ ratio *= 0.8; }
      else{ ratio *= 0.5; }
      if( ratio < 0.65 ) break;
    }
  }
  
  LaplaceDelaunaySmoothing2(aVec2, aPo2D,aTri,aflag_isnt_move);
}


bool TriangulateOuterLoop2
(std::vector<CEPo2>& aPo2D,
 std::vector<ETri>& aTri_in,
 std::vector<CVector2>& aVec2,
 const std::vector<int>& aPtrVtxInd,
 const std::vector<int>& aVtxInd)
{
  assert( aPo2D.size() == aVec2.size() );
  
  std::vector<ETri> aTri;
  std::vector<int> aPoDel;
  { // super triangle
    double max_len;
    double center[2];
    {
      double bound_2d[4];
      bound_2d[0] = aVec2[0].x;
      bound_2d[1] = aVec2[0].x;
      bound_2d[2] = aVec2[0].y;
      bound_2d[3] = aVec2[0].y;
      for(int ipoin=1;ipoin<(int)aPo2D.size();ipoin++){
        if( aVec2[ipoin].x < bound_2d[0] ){ bound_2d[0] = aVec2[ipoin].x; }
        if( aVec2[ipoin].x > bound_2d[1] ){ bound_2d[1] = aVec2[ipoin].x; }
        if( aVec2[ipoin].y < bound_2d[2] ){ bound_2d[2] = aVec2[ipoin].y; }
        if( aVec2[ipoin].y > bound_2d[3] ){ bound_2d[3] = aVec2[ipoin].y; }
      }
      max_len = (bound_2d[1]-bound_2d[0]>bound_2d[3]-bound_2d[2]) ? bound_2d[1]-bound_2d[0] : bound_2d[3]-bound_2d[2];
      center[0] = (bound_2d[1]+bound_2d[0])*0.5;
      center[1] = (bound_2d[3]+bound_2d[2])*0.5;
    }
    
    const double tri_len = max_len * 4.0;
    const double tmp_len = tri_len * sqrt(3.0) / 6.0;
    
    const int npo = (int)aPo2D.size();
    aPoDel.push_back( (int)aPo2D.size()+0 );
    aPoDel.push_back( (int)aPo2D.size()+1 );
    aPoDel.push_back( (int)aPo2D.size()+2 );
    aPo2D.resize(npo+3);
    aVec2.resize(npo+3);
    aVec2[npo+0] = CVector2(center[0], center[1]+2.0*tmp_len);
    aPo2D[npo+0].e = 0;	aPo2D[npo+0].d = 0;
    aVec2[npo+1] = CVector2(center[0]-0.5*tri_len, center[1]-tmp_len);
    aPo2D[npo+1].e = 0;	aPo2D[npo+1].d = 1;
    aVec2[npo+2] = CVector2(center[0]+0.5*tri_len, center[1]-tmp_len);
    aPo2D[npo+2].e = 0;	aPo2D[npo+2].d = 2;
    
    aTri.resize(1);
    aTri[0].v[0] = npo+0;
    aTri[0].v[1] = npo+1;
    aTri[0].v[2] = npo+2;
    aTri[0].s2[0] =  -1;
    aTri[0].s2[1] =  -1;
    aTri[0].s2[2] =  -1;
    aTri[0].r2[0] =  0;
    aTri[0].r2[1] =  0;
    aTri[0].r2[2] =  0;
  }
  
  const double MIN_TRI_AREA = 1.0e-10;
  // Make Delaunay Division
  for(int ipoin=0;ipoin<(int)aPo2D.size();ipoin++){
    if( aPo2D[ipoin].e >= 0 ) continue;	// already added
    const CVector2& po_add = aVec2[ipoin];
    int itri_in = -1;
    int iedge = -1;
    int iflg1 = 0, iflg2 = 0;
    for(int itri=0;itri<(int)aTri.size();itri++){
      iflg1 = 0; iflg2 = 0;
      const ETri& ref_tri = aTri[itri];
      if( TriArea(po_add, aVec2[ref_tri.v[1]], aVec2[ref_tri.v[2]] ) > MIN_TRI_AREA ){
        iflg1++; iflg2 += 0;
      }
      if( TriArea(po_add, aVec2[ref_tri.v[2]], aVec2[ref_tri.v[0]] ) > MIN_TRI_AREA ){
        iflg1++; iflg2 += 1;
      }
      if( TriArea(po_add, aVec2[ref_tri.v[0]], aVec2[ref_tri.v[1]] ) > MIN_TRI_AREA ){
        iflg1++; iflg2 += 2;
      }
      if( iflg1 == 3 ){
        itri_in = itri;
        break;
      }
      else if( iflg1 == 2 ){
        const int ied0 = 3-iflg2;
        const int ipo_e0 = ref_tri.v[ (ied0+1)%3 ];
        const int ipo_e1 = ref_tri.v[ (ied0+2)%3 ];
        const unsigned int* rel = relTriTri[ ref_tri.r2[ied0] ];
        const int itri_s = ref_tri.s2[ied0];
        assert( aTri[itri_s].v[ rel[ (ied0+1)%3 ] ] == ipo_e0 );
        assert( aTri[itri_s].v[ rel[ (ied0+2)%3 ] ] == ipo_e1 );
        const int inoel_d = rel[ied0];
        assert( aTri[itri_s].s2[inoel_d] == itri );
        const int ipo_d = aTri[itri_s].v[inoel_d];
        assert( TriArea( po_add, aVec2[ipo_e1], aVec2[ aTri[itri].v[ied0] ] ) > MIN_TRI_AREA );
        assert( TriArea( po_add, aVec2[ aTri[itri].v[ied0] ], aVec2[ipo_e0] ) > MIN_TRI_AREA );
        if( TriArea( po_add, aVec2[ipo_e0], aVec2[ipo_d ] ) < MIN_TRI_AREA ){ continue;	}
        if( TriArea( po_add, aVec2[ipo_d ], aVec2[ipo_e1] ) < MIN_TRI_AREA ){ continue; }
        const int det_d =  DetDelaunay(po_add,aVec2[ipo_e0],aVec2[ipo_e1],aVec2[ipo_d]);
        if( det_d == 2 || det_d == 1 ) continue;
        itri_in = itri;
        iedge = ied0;
        break;
      }
    }
    if( itri_in == -1 ){
      //			std::cout << "Super Triangle Failure " << ipoin << " " << po_add.x << " " << po_add.y << std::endl;
      //			std::cout << aTri.size() << std::endl;
      assert(0);
      return false;
    }
    if( iedge == -1 ){
      InsertPoint_Elem(ipoin,itri_in,aPo2D,aTri);
    }
    else{
      InsertPoint_ElemEdge(ipoin,itri_in,iedge,aPo2D,aTri);
    }
    DelaunayAroundPoint2(ipoin,aPo2D,aTri,aVec2);
  }
  
  int itri0_ker  = (int)aTri.size(); // one internal triangle
  { // enforce edge
    const int nloop = (int)aPtrVtxInd.size()-1;
    for(int iloop=0;iloop<nloop;iloop++){
      const int nbar = aPtrVtxInd[iloop+1]-aPtrVtxInd[iloop];
      for(int ibar=0;ibar<nbar;ibar++){
        for(;;){
          int ipoi0, ipoi1;
          {
            int ind0 = aPtrVtxInd[iloop];
            ipoi0 = aVtxInd[ind0+ibar];
            if( ibar != nbar-1 ){ ipoi1 = aVtxInd[ind0+ibar+1]; }
            else{ ipoi1 = aVtxInd[ind0]; }
          }
          assert( ipoi0 < (int)aPo2D.size() ); assert( ipoi1 < (int)aPo2D.size() );
          int itri0;
          int inotri0,inotri1;
          if( FindEdge(itri0,inotri0,inotri1,ipoi0,ipoi1,aPo2D,aTri) ){ // this edge divide outside and inside
            assert( inotri0 != inotri1 );
            assert( inotri0 < 3 );
            assert( inotri1 < 3 );
            assert( aTri[itri0].v[ inotri0 ] == ipoi0 );
            assert( aTri[itri0].v[ inotri1 ] == ipoi1 );
            const int ied0 = 3 - inotri0 - inotri1;
            {
              const int itri1 = aTri[itri0].s2[ied0];
              const int ied1 = relTriTri[ aTri[itri0].r2[ied0] ][ied0];
              assert( aTri[itri1].s2[ied1] == itri0 );
              aTri[itri1].s2[ied1] = -1;
              aTri[itri0].s2[ied0] = -1;
            }
            itri0_ker = itri0;
            break;
          }
          else{ // this edge is devided from connection outer triangle
            double ratio;
            if( !FindEdgePoint_AcrossEdge2(itri0,inotri0,inotri1,ratio,
                                           ipoi0,ipoi1,
                                           aPo2D,aTri,aVec2) ){ assert(0); }
            assert( ratio > -1.0e-20 && ratio < 1.0+1.0e-20 );
            assert( TriArea( aVec2[ipoi0], aVec2[ aTri[itri0].v[inotri0] ], aVec2[ipoi1] ) > 1.0e-20 );
            assert( TriArea( aVec2[ipoi0], aVec2[ipoi1], aVec2[ aTri[itri0].v[inotri1] ] ) > 1.0e-20 );
            //						std::cout << ratio << std::endl;
            if( ratio < 1.0e-20 ){
              assert(0);
              return false;
            }
            else if( ratio > 1.0 - 1.0e-10 ){
              assert(0);
              return false;
            }
            else{
              const int ied0 = 3 - inotri0 - inotri1;
              assert( aTri[itri0].s2[ied0] >= 0 );
              const int itri1 = aTri[itri0].s2[ied0];
              const int ied1 = relTriTri[ aTri[itri0].r2[ied0] ][ied0];
              assert( aTri[itri1].s2[ied1] >= itri0 );
              FlipEdge(itri0,ied0,aPo2D,aTri);
              continue;
            }
          }
        }
      }
    }
  }
  
  {
    aTri_in.clear();
    int ntri_in;
    std::vector<int> inout_flg;
    {
      inout_flg.resize(aTri.size(),-1);
      inout_flg[itri0_ker] = 0;
      ntri_in = 1;
      std::stack<int> ind_stack;
      ind_stack.push(itri0_ker);
      for(;;){
        if( ind_stack.empty() ) break;
        const int itri_cur = ind_stack.top();
        ind_stack.pop();
        for(int inotri=0;inotri<3;inotri++){
          if( aTri[itri_cur].s2[inotri] == -1 ) continue;
          const int itri_s = aTri[itri_cur].s2[inotri];
          if( inout_flg[itri_s] == -1 ){
            inout_flg[itri_s] = ntri_in;
            ntri_in++;
            ind_stack.push(itri_s);
          }
        }
      }
    }
    aTri_in.resize( ntri_in );
    for(int itri=0;itri<(int)aTri.size();itri++){
      if( inout_flg[itri] != -1 ){
        int itri_in = inout_flg[itri];
        assert( itri_in >= 0 && (int)itri_in < ntri_in );
        aTri_in[itri_in] = aTri[itri];
      }
    }
    for(int itri=0;itri<(int)aTri_in.size();itri++){
      for(int ifatri=0;ifatri<3;ifatri++){
        if( aTri_in[itri].s2[ifatri] == -1 ) continue;
        int itri_s0 = aTri_in[itri].s2[ifatri];
        assert( itri_s0 >= 0 && (int)itri_s0 < (int)aTri.size() );
        int itri_in_s0 = inout_flg[itri_s0];
        assert( itri_in_s0 >= 0 && (int)itri_in_s0 < (int)aTri_in.size() );
        aTri_in[itri].s2[ifatri] = itri_in_s0;
      }
    }
  }
  ////////////////
  std::vector<int> map_po_del;
  int npo_pos;
  {
    map_po_del.resize( aPo2D.size(), -1 );
    for(int ipo=0;ipo<(int)aPoDel.size();ipo++){
      map_po_del[ aPoDel[ipo] ] = -2;
    }
    npo_pos = 0;
    for(int ipo=0;ipo<(int)aPo2D.size();ipo++){
      if( map_po_del[ipo] == -2 ) continue;
      map_po_del[ipo] = npo_pos;
      npo_pos++;
    }
  }
  {
    std::vector<CEPo2> aPo_tmp = aPo2D;
    aPo2D.resize( npo_pos );
    for(int ipo=0;ipo<(int)map_po_del.size();ipo++){
      if( map_po_del[ipo] == -2 ) continue;
      int ipo1 = map_po_del[ipo];
      aPo2D[ipo1] = aPo_tmp[ipo];
    }
  }
  for(int itri=0;itri<(int)aTri_in.size();itri++){
    for(int ifatri=0;ifatri<3;ifatri++){
      assert( aTri_in[itri].v[ifatri] != -2 );
      const int ipo = aTri_in[itri].v[ifatri];
      aTri_in[itri].v[ifatri] = map_po_del[ipo];
      aPo2D[ipo].e = itri;
      aPo2D[ipo].d = ifatri;
    }
  }
  
  return true;
}


bool IsInclude_Loop
(const double co[],
 const int ixy_stt, const int ixy_end,
 const std::vector<double>& aXY)
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
      const double p0[2] = {aXY[ipo0*2+0],aXY[ipo0*2+1]};
      const double p1[2] = {aXY[ipo1*2+0],aXY[ipo1*2+1]};
      const double area0 = TriArea2D(co,codir,p0);
      const double area1 = TriArea2D(co,p1,codir);
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
      double po2[2] = { r0*aXY[ipo0*2  ]+r1*aXY[ipo1*2  ],  r0*aXY[ipo0*2+1]+r1*aXY[ipo1*2+1] };
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

bool isOuterShapeOK
(int nloop,
 const std::vector<double>& aXY,
 const std::vector<int>& aIndXYs)
{
  ////////////////////////////////
  // enter Input check section
  
  { // make sure every loop has at least 3 points
    for(int iloop=0;iloop<nloop;iloop++){
      if( aIndXYs[iloop+1]-aIndXYs[iloop] < 3 ) return false;
    }
  }
  {
    ////////////////////////////////
    // check inclusion of loops
    for(int iloop=1;iloop<nloop;iloop++){
      for(int ipo=aIndXYs[iloop];ipo<aIndXYs[iloop+1];ipo++){
        const double pi[2] = {aXY[ipo*2+0],aXY[ipo*2+1]};
        if( !IsInclude_Loop(pi,
                            aIndXYs[0],aIndXYs[1],
                            aXY) ) return false;
      }
    }
    // check inclusion
    for(int iloop=1;iloop<nloop;iloop++){
      for(int jloop=0;jloop<nloop;jloop++){
        if( iloop == jloop ) continue;
        for(int jpo=aIndXYs[jloop];jpo<aIndXYs[jloop+1];jpo++){
          const double pj[2] = {aXY[jpo*2+0],aXY[jpo*2+1]};
          if( IsInclude_Loop(pj,
                             aIndXYs[iloop],aIndXYs[iloop+1],
                             aXY) ) return false;
        }
      }
    }
  }
  { // check intersection
    bool is_intersect = false;
    for(int iloop=0;iloop<nloop;iloop++){
      const int nbar_i = aIndXYs[iloop+1]-aIndXYs[iloop];
      for(int ibar=0;ibar<nbar_i;ibar++){
        const int ipo0 = aIndXYs[iloop] + ibar;
        int ipo1 = aIndXYs[iloop] + ibar+1;
        if( ibar == nbar_i-1 ){ ipo1 = aIndXYs[iloop]; }
        const double pi0[2] = {aXY[ipo0*2+0],aXY[ipo0*2+1]};
        const double pi1[2] = {aXY[ipo1*2+0],aXY[ipo1*2+1]};
        const double xmax_i = ( pi1[0] > pi0[0] ) ? pi1[0] : pi0[0];
        const double xmin_i = ( pi1[0] < pi0[0] ) ? pi1[0] : pi0[0];
        const double ymax_i = ( pi1[1] > pi0[1] ) ? pi1[1] : pi0[1];
        const double ymin_i = ( pi1[1] < pi0[1] ) ? pi1[1] : pi0[1];
        for(int jbar=ibar+2;jbar<nbar_i;jbar++){
          const int jpo0 = aIndXYs[iloop] + jbar;
          int jpo1 = aIndXYs[iloop] + jbar+1;
          if( jbar == nbar_i-1 ){
            if( ibar == 0 ) continue;
            jpo1 = aIndXYs[iloop];
          }
          const double pj0[2] = {aXY[jpo0*2+0],aXY[jpo0*2+1]};
          const double pj1[2] = {aXY[jpo1*2+0],aXY[jpo1*2+1]};
//          const double xmax_j = ( xys[jpo1*2+0] > xys[jpo0*2+0] ) ? xys[jpo1*2+0] : xys[jpo0*2+0];
//          const double xmin_j = ( xys[jpo1*2+0] < xys[jpo0*2+0] ) ? xys[jpo1*2+0] : xys[jpo0*2+0];
//          const double ymax_j = ( xys[jpo1*2+1] > xys[jpo0*2+1] ) ? xys[jpo1*2+1] : xys[jpo0*2+1];
//          const double ymin_j = ( xys[jpo1*2+1] < xys[jpo0*2+1] ) ? xys[jpo1*2+1] : xys[jpo0*2+1];
          const double xmax_j = ( pj1[0] > pj0[0] ) ? pj1[0] : pj0[0];
          const double xmin_j = ( pj1[0] < pj0[0] ) ? pj1[0] : pj0[0];
          const double ymax_j = ( pj1[1] > pj0[1] ) ? pj1[1] : pj0[1];
          const double ymin_j = ( pj1[1] < pj0[1] ) ? pj1[1] : pj0[1];
          if( xmin_j > xmax_i || xmax_j < xmin_i ) continue;	// åçˆÇ™Ç†ÇËÇ¶Ç»Ç¢ÉpÉ^Å[ÉìÇèúäO
          if( ymin_j > ymax_i || ymax_j < ymin_i ) continue;	// è„Ç…ìØÇ∂
          if( IsCrossLines(pi0,pi1,  pj0,pj1) ){
            is_intersect = true;
            break;
          }
        }
        if( is_intersect ) break;
        for(int jloop=iloop+1;jloop<nloop;jloop++){
          const int nbar_j = aIndXYs[jloop+1]-aIndXYs[jloop];
          for(int jbar=0;jbar<nbar_j;jbar++){
            const int jpo0 = aIndXYs[jloop] + jbar;
            int jpo1 = aIndXYs[jloop] + jbar+1;
            if( jbar == nbar_j-1 ){ jpo1 = aIndXYs[jloop]; }
            const double pj0[2] = {aXY[jpo0*2+0],aXY[jpo0*2+1]};
            const double pj1[2] = {aXY[jpo1*2+0],aXY[jpo1*2+1]};
            const double xmax_j = ( pj1[0] > pj0[0] ) ? pj1[0] : pj0[0];
            const double xmin_j = ( pj1[0] < pj0[0] ) ? pj1[0] : pj0[0];
            const double ymax_j = ( pj1[1] > pj0[1] ) ? pj1[1] : pj0[1];
            const double ymin_j = ( pj1[1] < pj0[1] ) ? pj1[1] : pj0[1];
            //            const double xmax_j = ( xys[jpo1*2+0] > xys[jpo0*2+0] ) ? xys[jpo1*2+0] : xys[jpo0*2+0];
            //            const double xmin_j = ( xys[jpo1*2+0] < xys[jpo0*2+0] ) ? xys[jpo1*2+0] : xys[jpo0*2+0];
            //            const double ymax_j = ( xys[jpo1*2+1] > xys[jpo0*2+1] ) ? xys[jpo1*2+1] : xys[jpo0*2+1];
            //            const double ymin_j = ( xys[jpo1*2+1] < xys[jpo0*2+1] ) ? xys[jpo1*2+1] : xys[jpo0*2+1];
            if( xmin_j > xmax_i || xmax_j < xmin_i ) continue;	// åçˆÇ™Ç†ÇËÇ¶Ç»Ç¢ÉpÉ^Å[ÉìÇèúäO
            if( ymin_j > ymax_i || ymax_j < ymin_i ) continue;	// è„Ç…ìØÇ∂
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
  ////////////////////////////////////////////////
  return true;
}

int delaunay_triangulation2
(std::vector<int>& aTri_out,		// out
 std::vector<double>& aXY_out, // out
 std::vector<int>& aPtrVtxInd, // out
 std::vector<int>& aVtxInd, // out
 ////
 const bool is_add_point_boundary,
 const std::vector<int>& aIndXYs, // in
 const std::vector<double>& aXY_in, // ind
 const double max_edge_length, // ind
 const CMeshDensity& mesh_density) // ind
{
  const int nloop = (int)aIndXYs.size()-1;
 
  if( !isOuterShapeOK(nloop,aXY_in,aIndXYs) ) return false;
  
  int nxys_presum = aIndXYs[nloop];
  std::vector<CEPo2> aPo2D;
  std::vector<CVector2> aVec2;
  aPo2D.resize(nxys_presum);
  aVec2.resize(nxys_presum);
  for(int ixys=0;ixys<nxys_presum;ixys++){
    aVec2[ixys] = CVector2(aXY_in[ixys*2+0], aXY_in[ixys*2+1]);
    aPo2D[ixys].e = -1;
    aPo2D[ixys].d = -1;
  }
  
  ////////////////////////////////
  // resampling
  // no resampling edge if(max_edge_length < 0)
  std::vector< std::vector<int> > aPoInEd;
  aPoInEd.resize(nxys_presum);
  if( max_edge_length > 0 ){
    for(int iloop=0;iloop<nloop;++iloop){
      int nadd = 0;
      const int nbar = aIndXYs[iloop+1]-aIndXYs[iloop];
      for(int ibar=0;ibar<nbar;ibar++){
        int ipo0 = aIndXYs[iloop]+ibar;
        int ipo1 = aIndXYs[iloop]+ibar+1;
        if( ibar == nbar-1 ){ ipo1 = aIndXYs[iloop]; }
        const double len = Distance( aVec2[ipo0], aVec2[ipo1] );
        nadd = (int)(len / max_edge_length);
        if( nadd == 0 || !is_add_point_boundary ) continue;
        const int ndiv = nadd+1;
        const double delx = (aVec2[ipo1].x - aVec2[ipo0].x)/ndiv;
        const double dely = (aVec2[ipo1].y - aVec2[ipo0].y)/ndiv;
        for(int iadd=0;iadd<nadd;++iadd){
          const unsigned int ipo = (int)aPo2D.size();
          CVector2 v2;
          v2.x = aVec2[ipo0].x + delx*(iadd+1);
          v2.y = aVec2[ipo0].y + dely*(iadd+1);
          CEPo2 po;
          po.e = -1;
          po.d = -1;
          aPo2D.push_back(po);
          aVec2.push_back(v2);
          aPoInEd[ aIndXYs[iloop]+ibar ].push_back(ipo);
        }
      }
    }
  }
  
  ////////////////////////////////
  {
    aPtrVtxInd.resize(nloop+1);
    aPtrVtxInd[0] = 0;
    for(int iloop=0;iloop<nloop;++iloop){
      const int nbar0 = aIndXYs[iloop+1]-aIndXYs[iloop];
      int nbar1 = nbar0;
      for(int ibar=0;ibar<nbar0;ibar++){
        nbar1 += aPoInEd[ aIndXYs[iloop]+ibar].size();
      }
      aPtrVtxInd[iloop+1] = aPtrVtxInd[iloop] + nbar1;
    }
    // adding new vertices on the outline
    aVtxInd.resize(aPtrVtxInd[nloop]);
    {
      int ivtx0 = 0;
      for(int iloop=0;iloop<nloop;iloop++){
        double area_loop = 0;
        { // area of this loop
          CVector2 vtmp(0,0);
          const int nbar = aPtrVtxInd[iloop+1]-aPtrVtxInd[iloop];
          for(int ibar=0;ibar<nbar;ibar++){
            int ipo0 = aPtrVtxInd[iloop]+ibar;
            int ipo1 = aPtrVtxInd[iloop]+ibar+1;
            if( ibar == nbar-1 ){ ipo1 = aPtrVtxInd[iloop]; }
            area_loop += TriArea( vtmp, aVec2[ipo0], aVec2[ipo1] );
          }
        }
        const int nbar0 = aIndXYs[iloop+1]-aIndXYs[iloop];
        if( (area_loop > 0) == (iloop == 0) ){ // outer loop
          for(int ibar=0;ibar<nbar0;ibar++){
            int ie = aIndXYs[iloop] + ibar;
            const std::vector<int>& add = aPoInEd[ie];
            aVtxInd[ivtx0] = ie;	ivtx0++;
            for(unsigned int iadd=0;iadd<add.size();iadd++){
              aVtxInd[ivtx0] = add[iadd];	ivtx0++;
            }
          }
        }
        else{
          for(int ibar=0;ibar<nbar0;ibar++){ // inner loop
            int ie = aIndXYs[iloop+1] - 1 - ibar;
            const std::vector<int>& add = aPoInEd[ie];
            const int nadd = (int)add.size();
            for(int iadd=0;iadd<nadd;iadd++){
              aVtxInd[ivtx0] = add[nadd-1-iadd];	ivtx0++;
            }
            aVtxInd[ivtx0] = ie;	ivtx0++;
          }
        }
      }
    }
  }
  ////////////////////////////////
  std::vector<ETri> aTri_in;
  if( !TriangulateOuterLoop2(aPo2D,aTri_in,aVec2,    aPtrVtxInd, aVtxInd) ){
    return true;
  }
  if( max_edge_length > 0 ){
    MeshingInside2(aPo2D,aTri_in,aVec2, aVtxInd,max_edge_length,mesh_density);
  }
  
  ////////////////////////////////
  // pushing back to STL vector
  const int ntri = (int)aTri_in.size();
  aTri_out.resize(ntri*3);
  for(int itri=0;itri<ntri;itri++){
    aTri_out[itri*3+0] = aTri_in[itri].v[0];
    aTri_out[itri*3+1] = aTri_in[itri].v[1];
    aTri_out[itri*3+2] = aTri_in[itri].v[2];
  }
  const int nxy_out = (int)aPo2D.size();
  aXY_out.resize(nxy_out*2);
  for(int ixy=0;ixy<nxy_out;ixy++){
    aXY_out[ixy*2+0] = aVec2[ixy].x;
    aXY_out[ixy*2+1] = aVec2[ixy].y;
  }
  return true;
}


bool GenerateTesselation2
(std::vector<int>& aTri_out, // out
 std::vector<double>& aXY_out, // out
 std::vector<int>& aInd_InVtxLoop,
 std::vector<int>& aIndVtxLoop,
 double elen,
 const CMeshDensity& mesh_density, 
 bool is_resample,
 const std::vector< std::vector<double> >& aVecAry0) // in
{
  if( aVecAry0.empty() ) return false;
  if( aVecAry0[0].size() < 3 ) return false;
  
  ////////////////////////////////
  std::vector<double> aXY;
  std::vector<int> aIndXYs;
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
    aXY.resize(npo_sum*2);
    npo_sum = 0;
    for(int iloop=0;iloop<(int)nloop;iloop++){
      const int nxys = (int)aVecAry0[iloop].size()/2;
      for(int ixys=0;ixys<nxys;ixys++){
        aXY[npo_sum*2+0] = aVecAry0[iloop][ixys*2+0];
        aXY[npo_sum*2+1] = aVecAry0[iloop][ixys*2+1];
        npo_sum++;
      }
    }
  }
  aTri_out.clear();
  aXY_out.clear();
  aInd_InVtxLoop.clear();
  aIndVtxLoop.clear();
  int res = delaunay_triangulation2(aTri_out,aXY_out,
                                     aInd_InVtxLoop,aIndVtxLoop,
                                     is_resample,
                                     aIndXYs,aXY,
                                     elen,mesh_density);

  return res;
}

bool GenerateTesselation2
(std::vector<int>& aTri_out, // out
 std::vector<double>& aXY_out, // out
 std::vector<int>& aPtrVtxInd,
 std::vector<int>& aVtxInd,
 ////
 double elen,
 bool is_uniform_resample_loop, // good for polyline curve
 const std::vector< std::vector<double> >& aVecAry0) // in
{
  class CMeshDensity_Uni : public CMeshDensity {
      virtual double edgeLengthRatio(double px, double py) const {
        return 1.0;
      }
  } mdu;
  return GenerateTesselation2(aTri_out,aXY_out,
                              aPtrVtxInd,aVtxInd,
          ///
                              elen,
                              mdu,
                              is_uniform_resample_loop,
                              aVecAry0);
}


