#ifndef cadtopo_h
#define cadtopo_h

class CCadTopo
{
public:
  CCadTopo(){
    nVertex = 0;
  }
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
  bool AddPoint_Edge(int ieo){
    if( ieo < 0 || ieo >= (int)aEdge.size() ){ return false; }
    const int ivn = nVertex;
    nVertex += 1;
    const int iv0 = aEdge[ieo].iv0;
    const int iv1 = aEdge[ieo].iv1;
    const int ien = aEdge.size();
    aEdge.resize(aEdge.size()+1);
    aEdge[ieo].iv0 = iv0;
    aEdge[ieo].iv1 = ivn;
    aEdge[ien].iv0 = ivn;
    aEdge[ien].iv1 = iv1;
    for(unsigned int ifc=0;ifc<aFace.size();++ifc){
      const int ne = aFace[ifc].aIE.size();
      int iie = 0;
      for(;iie<ne;++iie){
        if( aFace[ifc].aIE[iie].first == ieo ){ break; }
      }
      if( iie == ne ){ continue; }
      if( aFace[ifc].aIE[iie].second ){
        aFace[ifc].aIE.insert(aFace[ifc].aIE.begin()+iie+1,std::make_pair(ien,true));
      }
      else{
        std::cout << "TODO: implement this" << std::endl;
      }
    }
    return true;
  }
  bool Check() const{
    for(unsigned int ifc=0;ifc<aFace.size();++ifc){
      const int ne = aFace[ifc].aIE.size();
      for(int iie=0;iie<ne;++iie){
        int ie0 = aFace[ifc].aIE[(iie+0)%ne].first;
        int ie1 = aFace[ifc].aIE[(iie+1)%ne].first;
        bool flg0 = aFace[ifc].aIE[(iie+0)%ne].second;
        bool flg1 = aFace[ifc].aIE[(iie+1)%ne].second;
        int iv0a = (flg0) ? aEdge[ie0].iv1 : aEdge[ie0].iv0;
        int iv0b = (flg1) ? aEdge[ie1].iv0 : aEdge[ie1].iv1;
        assert( iv0a == iv0b );
        if( iv0a != iv0b ) return false;
      }
    }
    return true;
  }
public:
  class CEdge{
  public:
    int iv0,iv1;
  };
  class CFace{
  public:
    std::vector<int> GetArray_IdVertex(const std::vector<CEdge>& aEdge) const
    {
      std::vector<int> res;
      for(unsigned int ie=0;ie<aIE.size();++ie){
        const int ie0 = aIE[ie].first;
        const bool dir = aIE[ie].second;
        const int iv0 = dir ? aEdge[ie0].iv0 : aEdge[ie0].iv1;
        res.push_back(iv0);
      }
      return res;
    }
  public:
    std::vector< std::pair<int,bool> > aIE; // index of edge, is this edge ccw?
  };
public:
  int nVertex;
  std::vector<CEdge> aEdge;
  std::vector<CFace> aFace;
};


#endif /* cadtopo_h */
