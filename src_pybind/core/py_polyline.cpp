#include <stdio.h>
#include <vector>
#include <deque>
#include <map>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/vec2.h"

namespace py = pybind11;

int getLevelSet_AddEdge
(int ip0, int ip1,
 double v0, double v1,
 std::map<std::pair<int,int>, std::pair<int,double> >& aNode)
{
  assert(v0*v1<0);
  auto key = std::make_pair(ip0,ip1);
  if( ip0 > ip1 ){
    key = std::make_pair(ip1,ip0);
  }
  auto res = aNode.find(key);
  if( res == aNode.end() ){
    int ino0 = aNode.size();
    auto val = std::make_pair(ino0,-v1/(v0-v1));
    aNode.insert(std::make_pair(key,val));
    return ino0;
  }
  else{
    return res->second.first;
  }
  return -1;
}

class CLevSetNode{
public:
  int ip0;
  int ip1;
  double ratio;
  int in0, in1;
public:
  CLevSetNode(){
    in0 = -1;
    in1 = -1;
  }
  void setNeighbour(int ino_in){
    if( in0 == - 1 ){ in0 = ino_in; }
    else{
      assert( in1 == -1 );
      in1 = ino_in;
    }
  }
  void pos(double& x0, double& y0, int nw){
    assert(ip0<ip1);
    int ih0 = ip0/nw;
    int iw0 = ip0-ih0*nw;
    int ih1 = ip1/nw;
    int iw1 = ip1-ih1*nw;
    x0 = ratio * iw0 + (1-ratio)*iw1;
    y0 = ratio * ih0 + (1-ratio)*ih1;
  }
};


void getLevelSet_GetEdges
(std::map<std::pair<int,int>, std::pair<int,double> >& mapNode,
std::vector< std::pair<int,int> >& aEdge,
 int nh, int nw, float thres, const float* pVal)
{
//  const int mh = nh+1;
  const int mw = nw+1;
  for(int ih=0;ih<nh;++ih){ // x
    for(int iw=0;iw<nw;++iw){ // y
      const int ip0 = (ih+0)*mw+(iw+0);
      const int ip1 = (ih+0)*mw+(iw+1);
      const int ip2 = (ih+1)*mw+(iw+0);
      const int ip3 = (ih+1)*mw+(iw+1);
      const float v0 = pVal[ip0]-thres;
      const float v1 = pVal[ip1]-thres;
      const float v2 = pVal[ip2]-thres;
      const float v3 = pVal[ip3]-thres;
      if( (v0>0 && v1<0 && v2<0 && v3<0) ||  (v0<0 && v1>0 && v2>0 && v3>0) ){
        int in0 = getLevelSet_AddEdge(ip0,ip1,v0,v1,mapNode);
        int in1 = getLevelSet_AddEdge(ip0,ip2,v0,v2,mapNode);
        aEdge.push_back(std::make_pair(in0,in1));
      }
      if( (v0<0 && v1>0 && v2<0 && v3<0) ||  (v0>0 && v1<0 && v2>0 && v3>0) ){
        int in0 = getLevelSet_AddEdge(ip0,ip1,v0,v1,mapNode);
        int in1 = getLevelSet_AddEdge(ip1,ip3,v1,v3,mapNode);
        aEdge.push_back(std::make_pair(in0,in1));
      }
      if( (v0<0 && v1<0 && v2>0 && v3<0) ||  (v0>0 && v1>0 && v2<0 && v3>0) ){
        int in0 = getLevelSet_AddEdge(ip0,ip2,v0,v2,mapNode);
        int in1 = getLevelSet_AddEdge(ip2,ip3,v2,v3,mapNode);
        aEdge.push_back(std::make_pair(in0,in1));
      }
      if( (v0<0 && v1<0 && v2<0 && v3>0) ||  (v0>0 && v1>0 && v2>0 && v3<0) ){
        int in0 = getLevelSet_AddEdge(ip1,ip3,v1,v3,mapNode);
        int in1 = getLevelSet_AddEdge(ip2,ip3,v2,v3,mapNode);
        aEdge.push_back(std::make_pair(in0,in1));
      }
      if( (v0>0 && v1>0 && v2<0 && v3<0) ||  (v0<0 && v1<0 && v2>0 && v3>0) ){
        int in0 = getLevelSet_AddEdge(ip0,ip2,v0,v2,mapNode);
        int in1 = getLevelSet_AddEdge(ip1,ip3,v1,v3,mapNode);
        aEdge.push_back(std::make_pair(in0,in1));
      }
      if( (v0>0 && v1<0 && v2>0 && v3<0) ||  (v0<0 && v1>0 && v2<0 && v3>0) ){
        int in0 = getLevelSet_AddEdge(ip0,ip1,v0,v1,mapNode);
        int in1 = getLevelSet_AddEdge(ip2,ip3,v2,v3,mapNode);
        aEdge.push_back(std::make_pair(in0,in1));
      }
      if( (v0>0 && v1<0 && v2<0 && v3>0) ||  (v0<0 && v1>0 && v2>0 && v3<0) ){
        int in0 = getLevelSet_AddEdge(ip0,ip1,v0,v1,mapNode);
        int in1 = getLevelSet_AddEdge(ip0,ip2,v0,v2,mapNode);
        aEdge.push_back(std::make_pair(in0,in1));
        int in2 = getLevelSet_AddEdge(ip1,ip3,v1,v3,mapNode);
        int in3 = getLevelSet_AddEdge(ip2,ip3,v2,v3,mapNode);
        aEdge.push_back(std::make_pair(in2,in3));
      }
    }
  }
}

void getLevelSet_MakeLoop
(std::vector< std::deque<int> >& aLoop,
 const std::vector<CLevSetNode>& aNode)
{
  unsigned int ino_ker = 0;
  for(;;){
    for(;ino_ker<aNode.size();++ino_ker){
      bool is_included = false;
      for(unsigned int iloop=0;iloop<aLoop.size();++iloop){
        for(unsigned int jino=0;jino<aLoop[iloop].size();++jino){
          const int jno0 = aLoop[iloop][jino];
          if( jno0 == (int)ino_ker ){ is_included = true; break; }
        }
        if( is_included ){ break; }
      }
      if( !is_included ){ break; }
    }
    if( ino_ker == aNode.size() ){ break; }
    //    std::cout << "inoker: " << ino_ker << " " << aNode.size() << std::endl;
    const int iloop0 = aLoop.size();
    aLoop.resize(aLoop.size()+1);
    aLoop[iloop0].push_back(ino_ker);
    if( aNode[ino_ker].in0 != -1 ){
      aLoop[iloop0].push_front(aNode[ino_ker].in0);
      for(;;){
        int in1 = aLoop[iloop0][0];
        int in2 = aLoop[iloop0][1];
        int jn0a = aNode[in1].in0;
        int jn0b = aNode[in1].in1;
        int in0 = -1;
        if( jn0a == in2 ){ in0 = jn0b; }
        if( jn0b == in2 ){ in0 = jn0a; }
        if( in0 == -1 ){ break; }
        if( in0 == aLoop[iloop0].back() ){ break; }
        aLoop[iloop0].push_front(in0);
      }
    }
    if( aLoop[iloop0].front() != aNode[ino_ker].in1 && aNode[ino_ker].in1 != -1 ){
      aLoop[iloop0].push_back(aNode[ino_ker].in1);
      for(;;){
        int in1 = aLoop[iloop0][aLoop[iloop0].size()-1];
        int in0 = aLoop[iloop0][aLoop[iloop0].size()-2];
        int jn2a = aNode[in1].in0;
        int jn2b = aNode[in1].in1;
        int in2 = -1;
        if( jn2a == in0 ){ in2 = jn2b; }
        if( jn2b == in0 ){ in2 = jn2a; }
        if( in2 == -1 ){ break; }
        if( in2 == aLoop[iloop0].front() ){ break; }
        aLoop[iloop0].push_back(in2);
      }
    }
  }
}

std::vector<std::vector<double> > getLevelSet
(int mh, int mw,
 const float* pVal, double thres)
{
  const int nh = mh-1;
  const int nw = mw-1;
  std::vector<CLevSetNode> aNode;
  {
    std::vector< std::pair<int,int> > aEdge;
    std::map<std::pair<int,int>, std::pair<int,double> > mapNode;
    getLevelSet_GetEdges(mapNode,aEdge,
                         nh,nw,thres,pVal);
    const int nNode = mapNode.size();
    aNode.resize(nNode);
    for(auto itr=mapNode.begin();itr!=mapNode.end();++itr){
      CLevSetNode node;
      node.ip0 = itr->first.first;
      node.ip1 = itr->first.second;
      node.ratio = itr->second.second;
      int in0 = itr->second.first;
      aNode[in0] = node;
    }
    for(unsigned int ie=0;ie<aEdge.size();++ie){
      const int in0 = aEdge[ie].first;
      const int in1 = aEdge[ie].second;
      aNode[in0].setNeighbour(in1);
      aNode[in1].setNeighbour(in0);
    }
  }
  ///////////////////////////////////
  std::vector< std::deque<int> > aLoop;
  getLevelSet_MakeLoop(aLoop,aNode);
  /////////
  std::vector< std::vector<double> > aRes;
  aRes.resize(aLoop.size());
  for(unsigned int iloop=0;iloop<aLoop.size();++iloop){
    for(unsigned int iin=0;iin<aLoop[iloop].size();++iin){
      int in0 = aLoop[iloop][iin];
      assert( in0>=0 && in0<(int)aNode.size() );
      double x0, y0; aNode[in0].pos(x0,y0,mw);
      aRes[iloop].push_back(x0);
      aRes[iloop].push_back(y0);
    }
  }
  return aRes;
}


std::vector<std::vector<int> > pyGetLevelSet(const py::array_t<float>& a, double mag)
{
  assert(a.ndim()==2);
  const int mh0 = a.shape()[0];
  const int mw0 = a.shape()[1];
  const int mh1 = mh0+2;
  const int mw1 = mw0+2;
  std::vector<float> a1(mh1*mw1,0.0);
  for(int ih=0;ih<mh0;++ih){
  for(int iw=0;iw<mw0;++iw){
    const int ip0 = (ih+0)*mw0+(iw+0);
    const int ip1 = (ih+1)*mw1+(iw+1);
    a1[ip1] = a.data()[ip0];
  }
  }
  std::vector<std::vector<double> > aLoop = getLevelSet(mh1, mw1, a1.data(), 0.5);
  std::vector<std::vector<int> > aLoopInt(aLoop.size());
  for(unsigned int iloop=0;iloop<aLoop.size();++iloop){
    std::vector<double> loop = aLoop[iloop];
    const int np = loop.size()/2;
    for(int ip=0;ip<np;++ip){
      double x0 = loop[ip*2+0];
      double y0 = loop[ip*2+1];
      double ix0 = int((x0+0.5-1.0)*mag);
      double iy0 = int((y0+0.5-1.0)*mag);
      aLoopInt[iloop].push_back(ix0);
      aLoopInt[iloop].push_back(iy0);
    }
  }
  return aLoopInt;
}


double distance
(int i0, int i1, int i2,
 const std::vector<int>& aXY)
{
  const int np = aXY.size()/2;
  assert(i0>=0&&i0<np);
  assert(i1>=0&&i1<np);
  assert(i2>=0&&i2<np);
  CVector2 v0(aXY[i0*2+0], aXY[i0*2+1]);
  CVector2 v1(aXY[i1*2+0], aXY[i1*2+1]);
  CVector2 v2(aXY[i2*2+0], aXY[i2*2+1]);
  CVector2 vn = GetNearest_LineSeg_Point(v2,v0,v1);
  return Distance(vn,v2);
}

void set_flag_douglas_peucker
(std::vector<bool>& aFlg,
 const std::vector<int>& aXY,
 int i0, int i1, double eps)
{
  assert(i0<i1);
  assert(i0>=0&&i0<(int)aFlg.size()&&aFlg[i0]);
  assert(i1>=0&&i1<(int)aFlg.size()&&aFlg[i1]);
  int ifar = -1;
  double max_dist = 0.0;
  for(int i2=i0+1;i2<i1;++i2){
    assert( aFlg[i2] );
    double dist = distance(i0,i1,i2,aXY);
    if( ifar == -1 || dist > max_dist ){
      max_dist = dist;
      ifar = i2;
    }
  }
  if( max_dist < eps ){
    for(int i2=i0+1;i2<i1;++i2){ aFlg[i2] = false; }
  }
  else{
    if( ifar-i0 >= 2 ){ set_flag_douglas_peucker(aFlg, aXY, i0, ifar, eps); }
    if( i1-ifar >= 2 ){ set_flag_douglas_peucker(aFlg, aXY, ifar, i1, eps); }
  }
}


std::vector<int> simplifyPolyloop(std::vector<int>& aXY_in, double eps)
{
  std::vector<int> aXY_out;
  int np = aXY_in.size()/2;
  if( np <=3 ){ return aXY_out; }
  std::vector<bool> aFlg(np,true);
  set_flag_douglas_peucker(aFlg, aXY_in, 0, np/2, eps);
  set_flag_douglas_peucker(aFlg, aXY_in, np/2,np-1, eps);
  for(int ip=0;ip<np;++ip){
    if( !aFlg[ip] ) continue;
    aXY_out.push_back(aXY_in[ip*2+0]);
    aXY_out.push_back(aXY_in[ip*2+1]);
  }
  return aXY_out;
}


void init_polyline(py::module &m){
  m.def("get_level_set", pyGetLevelSet);
  m.def("simplify_polyloop",simplifyPolyloop);
}
