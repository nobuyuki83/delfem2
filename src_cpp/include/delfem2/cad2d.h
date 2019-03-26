#ifndef CAD2D_h
#define CAD2D_h

#include "delfem2/vec2.h"

#include "delfem2/funcs_gl.h"

class CCadTopo
{
public:
  void Clear(){
    nVertex = 0;
    aEdge.clear();
    aFace.clear();
  }
  void AddPolygon(int np){
    const int iv0 = nVertex;
    nVertex += np;
    const int ie0 = aEdge.size();
    for(int iie=0;iie<np;++iie){
      CEdge edge0;
      edge0.iv0 = iv0 + (iie+0)%np;
      edge0.iv1 = iv0 + (iie+1)%np;
      aEdge.push_back(edge0);
    }
    CFace face0;
    for(int iie=0;iie<np;++iie){
      face0.aIE.push_back( std::make_pair(ie0+iie,true ) );
    }
    aFace.push_back(face0);
  }
public:
  class CEdge{
  public:
    int iv0,iv1;
  };
  class CFace{
  public:
    std::vector< std::pair<int,bool> > aIE; // index of edge, is this edge ccw?
  };
public:
  int nVertex;
  std::vector<CEdge> aEdge;
  std::vector<CFace> aFace;
};



class CCad2D_VtxGeo{
public:
  CCad2D_VtxGeo(const CVector2& p) : pos(p){}
public:
  CVector2 pos;
};
class CCad2D_EdgeGeo{
public:
  void GenMesh(unsigned int iedge, const CCadTopo& topo,
               std::vector<CCad2D_VtxGeo>& aVtxGeo);
  double Distance(double x, double y) const;
public:
  CVector2 p0,p1;
  std::vector<CVector2> aP;
};
class CCad2D_FaceGeo{
public:
  //    std::vector<CVector2> aP;
  std::vector<int> aTri;
  std::vector<double> aXY;
public:
  void GenMesh(unsigned int iface0, const CCadTopo& topo, 
               std::vector<CCad2D_EdgeGeo>& aEdgeGeo);
};

//////////////////

class CCad2D
{
public:
  CCad2D(){
    ivtx_picked = -1;
  }
  void Clear(){
    aVtx.clear();
    aEdge.clear();
    aFace.clear();
    topo.Clear();
  }
  void Draw() const;
  void Mouse(int btn,int action,int mods,
             const std::vector<double>& src,
             const std::vector<double>& dir,
             double view_height);
  void Motion(const std::vector<double>& src0,
              const std::vector<double>& src1,
              const std::vector<double>& dir);
  std::vector<double> MinMaxXYZ() const;
  void AddPolygon(const std::vector<double>& aXY);
  void Meshing(std::vector<double>& aXY,
               std::vector<int>& aTri,
               double len) const;
  void setBCFlagEdge(int* pBC,
                     const double* pXY, int nXY,
                     const std::vector<int>& aIE,
                     int iflag, double torelance) const;
public:
public:
public:
  CCadTopo topo;
  /////
  std::vector<CCad2D_VtxGeo> aVtx;
  std::vector<CCad2D_EdgeGeo> aEdge;
  std::vector<CCad2D_FaceGeo> aFace;
  int ivtx_picked;
};



#endif
