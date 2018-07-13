#include <stdio.h>
#include <set>
#include <math.h>

#include "delfem2/bvh.h"

double DetSide(const double p[3], const double org[3], const double n[3]){
  return (p[0]-org[0])*n[0] + (p[1]-org[1])*n[1] + (p[2]-org[2])*n[2];
}

void DevideElemAryConnex
(int iroot_node,
 std::vector<int>& aElem2Node,
 std::vector<CNodeBVH>& aNodeBVH,
 ////
 const std::vector<int>& list,
 const std::vector<int>& aElemSurInd,
 const std::vector<int>& aElemSur,
 const std::vector<double>& aElemCenter)
{
  assert( list.size() > 1 );
  double eps = 1.0e-10;
  double x_min,x_max, y_min,y_max, z_min,z_max;
  {
    {
      assert( list.size() > 0 );
      int itri = list[0];
      assert( aElem2Node[itri] == iroot_node );
      double cgx = aElemCenter[itri*3+0];
      double cgy = aElemCenter[itri*3+1];
      double cgz = aElemCenter[itri*3+2];
      x_min = cgx-eps;
      x_max = cgx+eps;
      y_min = cgy-eps;
      y_max = cgy+eps;
      z_min = cgz-eps;
      z_max = cgz+eps;
    }
    for(int il=1;il<list.size();il++){
      int itri = list[il];
      assert( itri < aElemCenter.size() );
      assert( aElem2Node[itri] == iroot_node );
      double cgx = aElemCenter[itri*3+0];
      double cgy = aElemCenter[itri*3+1];
      double cgz = aElemCenter[itri*3+2];
      x_min = ( cgx-eps < x_min ) ? cgx-eps : x_min;
      x_max = ( cgx+eps > x_max ) ? cgx+eps : x_max;
      y_min = ( cgy-eps < y_min ) ? cgy-eps : y_min;
      y_max = ( cgy+eps > y_max ) ? cgy+eps : y_max;
      z_min = ( cgz-eps < z_min ) ? cgz-eps : z_min;
      z_max = ( cgz+eps > z_max ) ? cgz+eps : z_max;
    }
  }
  double dir[3] = {0,0,0}; // longest direction of AABB
  {
    double lenx = x_max - x_min;
    double leny = y_max - y_min;
    double lenz = z_max - z_min;
    if( lenx > leny && lenx > lenz ){ dir[0] = 1; }
    if( leny > lenz && leny > lenx ){ dir[1] = 1; }
    if( lenz > lenx && lenz > leny ){ dir[2] = 1; }
  }
  double org[3] = {(x_min+x_max)*0.5,  (y_min+y_max)*0.5,  (z_min+z_max)*0.5};
  int itri_ker = -1;
  for(int il=0;il<list.size();il++){
    int itri0 = list[il];
    const double det0 = DetSide(aElemCenter.data()+itri0*3,org,dir);
    if( fabs(det0) < 1.0e-10 ) continue;
    if( det0 < 0 ){ dir[0]*=-1; dir[1]*=-1; dir[2]*=-1; }
    itri_ker = itri0;
    break;
  }
  if( itri_ker == -1 ){
    org[0]=0; org[1]=0; org[2]=0;
    for(int il=0;il<list.size();il++){ // center of the gravity of list
      int itri0 = list[il];
      org[0] += aElemCenter[itri0*3+0];
      org[1] += aElemCenter[itri0*3+1];
      org[2] += aElemCenter[itri0*3+2];
    }
    org[0] = org[0]/list.size();
    org[1] = org[1]/list.size();
    org[2] = org[2]/list.size();
    double mat[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
    for(int il=0;il<list.size();il++){
      int itri0 = list[il];
      const double vcg[3] = {
        aElemCenter[itri0*3+0]-org[0],
        aElemCenter[itri0*3+1]-org[1],
        aElemCenter[itri0*3+2]-org[2] };
      for(int i=0;i<3;i++){ for(int j=0;j<3;j++){ mat[i][j] += vcg[i]*vcg[j]; } }
    }
    dir[0] = 1; dir[1] = 1; dir[2] = 1;
    for(int i=0;i<10;i++){ // power method to find the max eigen value/vector
      double tmp[3] = {
        mat[0][0]*dir[0] + mat[0][1]*dir[1] + mat[0][2]*dir[2],
        mat[1][0]*dir[0] + mat[1][1]*dir[1] + mat[1][2]*dir[2],
        mat[2][0]*dir[0] + mat[2][1]*dir[1] + mat[2][2]*dir[2] };
      double len = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
      dir[0] = tmp[0]/len;
      dir[1] = tmp[1]/len;
      dir[2] = tmp[2]/len;
    }
    for(int il=0;il<list.size();il++){
      int itri0 = list[il];
      double det = DetSide(aElemCenter.data()+itri0*3,org,dir);
      if( fabs(det) < 1.0e-10 ) continue;
      if( det < 0 ){ dir[0]*=-1; dir[1]*=-1; dir[2]*=-1; }
      itri_ker = itri0;
      break;
    }
  }
  int inode_ch0 = (int)aNodeBVH.size();
  int inode_ch1 = (int)aNodeBVH.size()+1;
  aNodeBVH.resize(aNodeBVH.size()+2);
  aNodeBVH[inode_ch0].iroot = iroot_node;
  aNodeBVH[inode_ch1].iroot = iroot_node;
  aNodeBVH[iroot_node].ichild[0] = inode_ch0;
  aNodeBVH[iroot_node].ichild[1] = inode_ch1;
  std::vector<int> list_ch0;
  { // 子ノード０に含まれる三角形を抽出（itri_kerと接続していて，dirベクトルの正方向にある三角形）
    aElem2Node[itri_ker] = inode_ch0;
    list_ch0.push_back(itri_ker);
    std::stack<int> stack;
    stack.push(itri_ker);
    while(!stack.empty()){
      int itri0 = stack.top();
      stack.pop();
      for(int jje=aElemSurInd[itri0];jje<aElemSurInd[itri0+1];jje++){
        int jtri = aElemSur[jje*2+0];
        if( jtri == -1 ) continue;
        if( aElem2Node[jtri] != iroot_node ) continue;
        assert( jtri < aElemCenter.size() );
        double det = DetSide(aElemCenter.data()+jtri*3,org,dir);
        if( det < 0 ) continue;
        stack.push(jtri);
        aElem2Node[jtri] = inode_ch0;
        list_ch0.push_back(jtri);
      }
    }
    assert( list_ch0.size() > 0 );
  }
  // 子ノード１に含まれる三角形を抽出(入力リストから子ノード0に含まれる物を除外)
  std::vector<int> list_ch1;
  for(int il=0;il<list.size();il++){
    int itri = list[il];
    if( aElem2Node[itri] == inode_ch0 ) continue;
    assert( aElem2Node[itri] == iroot_node );
    aElem2Node[itri] = inode_ch1;
    list_ch1.push_back(itri);
  }
  assert( list_ch1.size() > 0 );
  
  /////
  if( list_ch0.size() == 1 ){
    aNodeBVH[inode_ch0].ichild[0] = list_ch0[0];
    aNodeBVH[inode_ch0].ichild[1] = -1;
  }
  else{ // 子ノード0にある三角形を再度分割
    DevideElemAryConnex(inode_ch0,aElem2Node,aNodeBVH,
                       list_ch0,aElemSurInd,aElemSur,aElemCenter);
  }
  list_ch0.clear();
  
  //////
  if( list_ch1.size() == 1 ){
    aNodeBVH[inode_ch1].ichild[0] = list_ch1[0];
    aNodeBVH[inode_ch1].ichild[1] = -1;
  }
  else{ // 子ノード1にある三角形を再度分割
    DevideElemAryConnex(inode_ch1,aElem2Node,aNodeBVH,
                       list_ch1,aElemSurInd,aElemSur,aElemCenter);
  }
}

int MakeTreeTopologyBVH_TopDown
(std::vector<CNodeBVH>& aNodeBVH,
 const std::vector<int>& aElemSurInd,
 const std::vector<int>& aElemSur,
 const std::vector<double>& aElemCenter)
{
  aNodeBVH.clear();
  const int nelem = (int)aElemCenter.size()/3;
  assert( aElemSurInd.size() == nelem+1 );
  std::vector<int> list(nelem);
  for(int ielem=0;ielem<nelem;ielem++){ list[ielem] = ielem; }
  std::vector<int> aElem2Node;
  aElem2Node.resize(nelem,0);
  aNodeBVH.resize(1);
  aNodeBVH[0].iroot = -1;
  DevideElemAryConnex(0,aElem2Node,aNodeBVH,
                      list,aElemSurInd,aElemSur,aElemCenter);
  return 0;
}


// end independent from type of bounding volume
///////////////////////////////////////////////////////////////////////////////////////////

