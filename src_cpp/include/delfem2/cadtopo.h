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
