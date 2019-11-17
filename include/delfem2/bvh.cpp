/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <stdio.h>
#include <set>
#include <math.h>

#include "delfem2/bvh.h"

namespace dfm2 = delfem2;

// ------------------------------------

double DetSide(const double p[3], const double org[3], const double n[3]){
  return (p[0]-org[0])*n[0] + (p[1]-org[1])*n[1] + (p[2]-org[2])*n[2];
}

void DevideElemAryConnex
(int iroot_node,
 std::vector<int>& aElem2Node,
 std::vector<dfm2::CNodeBVH2>& aNodeBVH,
 ////
 const std::vector<int>& list,
 const int nfael,
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
    for(unsigned int il=1;il<list.size();il++){
      int itri = list[il];
      assert( itri < (int)aElemCenter.size() );
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
  for(unsigned int il=0;il<list.size();il++){
    int itri0 = list[il];
    const double det0 = DetSide(aElemCenter.data()+itri0*3,org,dir);
    if( fabs(det0) < 1.0e-10 ) continue;
    if( det0 < 0 ){ dir[0]*=-1; dir[1]*=-1; dir[2]*=-1; }
    itri_ker = itri0;
    break;
  }
  if( itri_ker == -1 ){
    org[0]=0; org[1]=0; org[2]=0;
    for(unsigned int il=0;il<list.size();il++){ // center of the gravity of list
      int itri0 = list[il];
      org[0] += aElemCenter[itri0*3+0];
      org[1] += aElemCenter[itri0*3+1];
      org[2] += aElemCenter[itri0*3+2];
    }
    org[0] = org[0]/list.size();
    org[1] = org[1]/list.size();
    org[2] = org[2]/list.size();
    double mat[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
    for(unsigned int il=0;il<list.size();il++){
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
    for(unsigned int il=0;il<list.size();il++){
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
      for(int ifael=0;ifael<nfael;ifael++){
        int jtri = aElemSur[itri0*6+nfael*2+0];
        if( jtri == -1 ) continue;
        if( aElem2Node[jtri] != iroot_node ) continue;
        assert( jtri < (int)aElemCenter.size() );
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
  for(unsigned int il=0;il<list.size();il++){
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
                       list_ch0,nfael,aElemSur,aElemCenter);
  }
  list_ch0.clear();
  
  //////
  if( list_ch1.size() == 1 ){
    aNodeBVH[inode_ch1].ichild[0] = list_ch1[0];
    aNodeBVH[inode_ch1].ichild[1] = -1;
  }
  else{ // 子ノード1にある三角形を再度分割
    DevideElemAryConnex(inode_ch1,aElem2Node,aNodeBVH,
                       list_ch1,nfael,aElemSur,aElemCenter);
  }
}

int dfm2::BVH_MakeTreeTopology
(std::vector<dfm2::CNodeBVH2>& aNodeBVH,
 const int nfael,
 const std::vector<int>& aElemSur,
 const std::vector<double>& aElemCenter)
{
  aNodeBVH.clear();
  const unsigned int nelem = aElemCenter.size()/3;
  std::vector<int> list(nelem);
  for(unsigned int ielem=0;ielem<nelem;ielem++){ list[ielem] = ielem; }
  std::vector<int> aElem2Node;
  aElem2Node.resize(nelem,0);
  aNodeBVH.resize(1);
  aNodeBVH[0].iroot = -1;
  DevideElemAryConnex(0,aElem2Node,aNodeBVH,
                      list,nfael,aElemSur,aElemCenter);
  return 0;
}



  // -------------------------------------------

  // Expands a 10-bit integer into 30 bits
  // by puting two zeros before each bit
  // "1011011111" -> "001000001001000001001001001001"
std::uint32_t expandBits(std::uint32_t v)
{
  v = (v * 0x00010001u) & 0xFF0000FFu;
  v = (v * 0x00000101u) & 0x0F00F00Fu;
  v = (v * 0x00000011u) & 0xC30C30C3u;
  v = (v * 0x00000005u) & 0x49249249u;
  return v;
}

std::uint32_t delfem2::MortonCode(double x, double y, double z)
{
  std::uint32_t ix = (unsigned int)fmin(fmax(x * 1024.0f, 0.0f), 1023.0f);
  std::uint32_t iy = (unsigned int)fmin(fmax(y * 1024.0f, 0.0f), 1023.0f);
  std::uint32_t iz = (unsigned int)fmin(fmax(z * 1024.0f, 0.0f), 1023.0f);
    //  std::cout << std::bitset<10>(ix) << " " << std::bitset<10>(iy) << " " << std::bitset<10>(iz) << std::endl;
  ix = expandBits(ix);
  iy = expandBits(iy);
  iz = expandBits(iz);
    //  std::cout << std::bitset<30>(ix) << " " << std::bitset<30>(iy) << " " << std::bitset<30>(iz) << std::endl;
  std::uint32_t ixyz = ix * 4 + iy * 2 + iz;
  return ixyz;
}

inline unsigned int clz(uint32_t x){
#ifdef __GNUC__ // GCC compiler
  return __builtin_clz(x);
#else
  int y = x;
  unsigned int n = 0;
  if (y == 0) return sizeof(y) * 8;
  while (1) {
    if (y < 0) break;
    n ++;
    y <<= 1;
  }
  return n;
#endif
}

int delta(int i, int j, const unsigned int* sorted_morton_code, int length)
{
  if (j<0 || j >= length){
    return -1;
  }
  else{
    return clz(sorted_morton_code[i] ^ sorted_morton_code[j]);
  }
}

std::pair<int,int> delfem2::determineRange
(const unsigned int* sorted_morton_code, int numInternalNode, int i)
{
  int size = numInternalNode + 1;
  int d = delta(i, i + 1, sorted_morton_code, size) - delta(i, i - 1, sorted_morton_code, size);
  d = d > 0 ? 1 : -1;
  
    //compute the upper bound for the length of the range
  int delta_min = delta(i, i - d, sorted_morton_code, size);
  int lmax = 2;
  while (delta(i, i + lmax*d, sorted_morton_code, size)>delta_min)
  {
    lmax = lmax * 2;
  }
  
    //find the other end using binary search
  int l = 0;
  for (int t = lmax / 2; t >= 1; t /= 2)
  {
    if (delta(i, i + (l + t)*d, sorted_morton_code, size)>delta_min)
    {
      l = l + t;
    }
  }
  int j = i + l*d;
  
  std::pair<int,int> range;
  if (i <= j) { range.first = i; range.second = j; }
  else { range.first = j; range.second = i; }
  return range;
}

bool is_diff_at_bit(unsigned int val1, unsigned int val2, int n)
{
  return val1 >> (31 - n) != val2 >> (31 - n);
}

int delfem2::findSplit(const unsigned int* sorted_morton_code, int start, int last)
{
    //return -1 if there is only
    //one primitive under this node.
  if (start == last)
  {
    return -1;
  }
  else
  {
    int common_prefix = clz(sorted_morton_code[start] ^ sorted_morton_code[last]);
    
      //handle duplicated morton code separately
    if (common_prefix == 32)
    {
      return (start + last) / 2;
    }
    
      // Use binary search to find where the next bit differs.
      // Specifically, we are looking for the highest object that
      // shares more than commonPrefix bits with the first one.
    
    int split = start; // initial guess
    int step = last - start;
    do
    {
      step = (step + 1) >> 1; // exponential decrease
      int newSplit = split + step; // proposed new position
      
      if (newSplit < last)
      {
        bool is_diff = is_diff_at_bit(sorted_morton_code[start],
                                      sorted_morton_code[newSplit],
                                      common_prefix);
        if (!is_diff)
        {
          split = newSplit; // accept proposal
        }
      }
    } while (step > 1);
    
    return split;
  }
}


class CPairMtcInd{
public:
  std::uint32_t imtc;
  std::uint32_t iobj;
public:
  bool operator < (const CPairMtcInd& rhs) const {
    return this->imtc < rhs.imtc;
  }
};


void dfm2::GetSortedMortenCode
 (std::vector<unsigned int>& aSortedId,
  std::vector<unsigned int>& aSortedMc,
  const std::vector<double>& aXYZ,
  const double minmax_xyz[6])
{
  std::vector<CPairMtcInd> aNodeBVH; // array of BVH node
  const int np = aXYZ.size()/3;
  aNodeBVH.resize(np);
  const double x_min = minmax_xyz[0];
  const double x_max = minmax_xyz[1];
  const double y_min = minmax_xyz[2];
  const double y_max = minmax_xyz[3];
  const double z_min = minmax_xyz[4];
  const double z_max = minmax_xyz[5];
  for(int ip=0;ip<np;++ip){
    double x = (aXYZ[ip*3+0]-x_min)/(x_max-x_min);
    double y = (aXYZ[ip*3+1]-y_min)/(y_max-y_min);
    double z = (aXYZ[ip*3+2]-z_min)/(z_max-z_min);
    aNodeBVH[ip].imtc = dfm2::MortonCode(x,y, z);
    aNodeBVH[ip].iobj = ip;
  }
  std::sort(aNodeBVH.begin(), aNodeBVH.end());
  aSortedId.resize(aNodeBVH.size());
  aSortedMc.resize(aNodeBVH.size());
  for(size_t ino=0;ino<aNodeBVH.size();++ino){
    aSortedMc[ino] = aNodeBVH[ino].imtc;
    aSortedId[ino] = aNodeBVH[ino].iobj;
      //        std::cout << std::bitset<32>(aNodeBVH[ino].imtc) << "  " << clz(aNodeBVH[ino].imtc) << "   " << ino << std::endl;
  }
}

void dfm2::BVH_TreeTopology_Morton
(std::vector<dfm2::CNodeBVH2>& aNodeBVH,
 const std::vector<unsigned int>& aSortedId,
 const std::vector<unsigned int>& aSortedMc)
{
  aNodeBVH.resize(aSortedMc.size()*2-1);
  aNodeBVH[0].iroot = -1;
  const unsigned int nni = aSortedMc.size()-1; // number of internal node
  for(int ini=0;ini<nni;++ini){
    const std::pair<int,int> range = dfm2::determineRange(aSortedMc.data(), aSortedMc.size()-1, ini);
    int isplit = dfm2::findSplit(aSortedMc.data(), range.first, range.second);
    assert( isplit != -1 );
    if( range.first == isplit ){
      const int inlA = nni+isplit;
      aNodeBVH[ini].ichild[0] = inlA;
      aNodeBVH[inlA].iroot = ini;
      aNodeBVH[inlA].ichild[0] = aSortedId[ini];
      aNodeBVH[inlA].ichild[1] = -1;
    }
    else{
      const int iniA = isplit;
      aNodeBVH[ini].ichild[0] = iniA;
      aNodeBVH[iniA].iroot = ini;
    }
      // ----
    if( range.second == isplit+1 ){
      const int inlB = nni+isplit+1;
      aNodeBVH[ini].ichild[1] = inlB;
      aNodeBVH[inlB].iroot = ini;
      aNodeBVH[inlB].ichild[0] = aSortedId[ini];
      aNodeBVH[inlB].ichild[1] = -1;
    }
    else{
      const int iniB = isplit+1;
      aNodeBVH[ini].ichild[1] = iniB;
      aNodeBVH[iniB].iroot = ini;
    }
  }
}
