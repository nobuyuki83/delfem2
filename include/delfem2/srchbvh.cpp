/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <algorithm>
#include <climits> // UNINT_MAX

#include "delfem2/srchbvh.h"

// ------------------------------------

namespace delfem2 {
namespace bvh {

DFM2_INLINE double DetSide(const double p[3], const double org[3], const double n[3]){
  return (p[0]-org[0])*n[0] + (p[1]-org[1])*n[1] + (p[2]-org[2])*n[2];
}

DFM2_INLINE void DevideElemAryConnex
(unsigned int iroot_node,
 std::vector<unsigned int>& aElem2Node,
 std::vector<CNodeBVH2>& aNodeBVH,
 //
 const std::vector<int>& list,
 const int nfael,
 const std::vector<unsigned int>& aElSuEl,
 const std::vector<double>& aElemCenter)
{
  assert( list.size() > 1 );
  double eps = 1.0e-10;
  double x_min,x_max, y_min,y_max, z_min,z_max;
  {
    {
      assert( !list.empty() );
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
    for(std::size_t il=1;il<list.size();il++){
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
  for(int itri0 : list){
    const double det0 = DetSide(aElemCenter.data()+itri0*3,org,dir);
    if( fabs(det0) < 1.0e-10 ) continue;
    if( det0 < 0 ){ dir[0]*=-1; dir[1]*=-1; dir[2]*=-1; }
    itri_ker = itri0;
    break;
  }
  if( itri_ker == -1 ){
    org[0]=0; org[1]=0; org[2]=0;
    for(int itri0 : list){ // center of the gravity of list
      org[0] += aElemCenter[itri0*3+0];
      org[1] += aElemCenter[itri0*3+1];
      org[2] += aElemCenter[itri0*3+2];
    }
    org[0] = org[0]/list.size();
    org[1] = org[1]/list.size();
    org[2] = org[2]/list.size();
    double mat[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
    for(int itri0 : list){
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
    for(int itri0 : list){
      double det = DetSide(aElemCenter.data()+itri0*3,org,dir);
      if( fabs(det) < 1.0e-10 ) continue;
      if( det < 0 ){ dir[0]*=-1; dir[1]*=-1; dir[2]*=-1; }
      itri_ker = itri0;
      break;
    }
  }
  const unsigned int inode_ch0 = static_cast<unsigned int>(aNodeBVH.size());
  const unsigned int inode_ch1 = static_cast<unsigned int>(aNodeBVH.size()+1);
  aNodeBVH.resize(aNodeBVH.size()+2);
  aNodeBVH[inode_ch0].iparent = iroot_node;
  aNodeBVH[inode_ch1].iparent = iroot_node;
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
        const unsigned int jtri = aElSuEl[itri0*nfael+ifael];
        if( jtri == UINT_MAX ) continue;
        if( aElem2Node[jtri] != iroot_node ) continue;
        assert( jtri < aElemCenter.size() );
        double det = DetSide(aElemCenter.data()+jtri*3,org,dir);
        if( det < 0 ) continue;
        stack.push(jtri);
        aElem2Node[jtri] = inode_ch0;
        list_ch0.push_back(jtri);
      }
    }
    assert( !list_ch0.empty() );
  }
  // 子ノード１に含まれる三角形を抽出(入力リストから子ノード0に含まれる物を除外)
  std::vector<int> list_ch1;
  for(int itri : list){
    if( aElem2Node[itri] == inode_ch0 ) continue;
    assert( aElem2Node[itri] == iroot_node );
    aElem2Node[itri] = inode_ch1;
    list_ch1.push_back(itri);
  }
  assert( !list_ch1.empty() );
  // ---------------------------
  if( list_ch0.size() == 1 ){
    aNodeBVH[inode_ch0].ichild[0] = list_ch0[0];
    aNodeBVH[inode_ch0].ichild[1] = UINT_MAX;
  }
  else{ // 子ノード0にある三角形を再度分割
    DevideElemAryConnex(inode_ch0,aElem2Node,aNodeBVH,
                       list_ch0,nfael,aElSuEl,aElemCenter);
  }
  list_ch0.clear();
  // -----------------------------
  if( list_ch1.size() == 1 ){
    aNodeBVH[inode_ch1].ichild[0] = list_ch1[0];
    aNodeBVH[inode_ch1].ichild[1] = UINT_MAX;
  }
  else{ // 子ノード1にある三角形を再度分割
    DevideElemAryConnex(inode_ch1,aElem2Node,aNodeBVH,
                        list_ch1,nfael,aElSuEl,aElemCenter);
  }
}

DFM2_INLINE int delta(
	int i,
	int j,
	const unsigned int* sorted_morton_code, 
	size_t length)
{
  if (j<0 || j >= (int)length){
    return -1;
  }
  else{
    return nbits_leading_zero(sorted_morton_code[i] ^ sorted_morton_code[j]);
  }
}

// Expands a 10-bit integer into 30 bits
// by puting two zeros before each bit
// "1011011111" -> "001000001001000001001001001001"
DFM2_INLINE std::uint32_t expandBits(std::uint32_t v)
{
  v = (v * 0x00010001u) & 0xFF0000FFu;
  v = (v * 0x00000101u) & 0x0F00F00Fu;
  v = (v * 0x00000011u) & 0xC30C30C3u;
  v = (v * 0x00000005u) & 0x49249249u;
  return v;
}

class CPairMtcInd{
public:
  std::uint32_t imtc;
  unsigned int iobj;
public:
  bool operator < (const CPairMtcInd& rhs) const {
    return (this->imtc < rhs.imtc);
  }
};


DFM2_INLINE void mark_child(
	std::vector<int>& aFlg,
	unsigned int inode0,                      
	const std::vector<CNodeBVH2>& aNode)
{
  assert( inode0 < aNode.size() );
  if( aNode[inode0].ichild[1] == UINT_MAX ){ // leaf
    const unsigned int in0 = aNode[inode0].ichild[0];
    assert( in0 < aFlg.size() );
    aFlg[in0] += 1;
    return;
  }
  const unsigned int in0 = aNode[inode0].ichild[0];
  const unsigned int in1 = aNode[inode0].ichild[1];
  mark_child(aFlg, in0, aNode);
  mark_child(aFlg, in1, aNode);
}

}
}

// ===========================================================

/**
 * @brief compute number of leading zeros
 * @function compute number of leading zeros
 * @param x input
 * @details clz(0) needs to be 32 to run BVH
 */
DFM2_INLINE unsigned int delfem2::nbits_leading_zero(
	uint32_t x){
  // avoid using buit-in functions such as "__builtin_clz(x)",
  // becuase application to 0 is typically undefiend.
  int32_t y = x;
  unsigned int n = 0;
  if (y == 0){ return sizeof(y) * 8; }
  while (true) {
    if (y < 0){ break; }
    n ++;
    y <<= 1;
  }
  return n;
}

DFM2_INLINE int delfem2::BVHTopology_TopDown_MeshElem(
    std::vector<CNodeBVH2>& aNodeBVH,
    const unsigned int nfael,
    const std::vector<unsigned int>& aElSuEl,
    const std::vector<double>& aElemCenter)
{
  aNodeBVH.clear();
  const size_t nelem = aElemCenter.size()/3;
  std::vector<int> list(nelem);
  for(unsigned int ielem=0;ielem<nelem;ielem++){ list[ielem] = ielem; }
  std::vector<unsigned int> aElem2Node;
  aElem2Node.resize(nelem,0);
  aNodeBVH.resize(1);
  aNodeBVH[0].iparent = UINT_MAX;
  bvh::DevideElemAryConnex(0,aElem2Node,aNodeBVH,
                           list,nfael,aElSuEl,aElemCenter);
  return 0;
}


template <typename REAL>
DFM2_INLINE std::uint32_t delfem2::MortonCode(REAL x, REAL y, REAL z)
{
  auto ix = (std::uint32_t)fmin(fmax(x * 1024.0f, 0.0f), 1023.0f);
  auto iy = (std::uint32_t)fmin(fmax(y * 1024.0f, 0.0f), 1023.0f);
  auto iz = (std::uint32_t)fmin(fmax(z * 1024.0f, 0.0f), 1023.0f);
    //  std::cout << std::bitset<10>(ix) << " " << std::bitset<10>(iy) << " " << std::bitset<10>(iz) << std::endl;
  ix = bvh::expandBits(ix);
  iy = bvh::expandBits(iy);
  iz = bvh::expandBits(iz);
    //  std::cout << std::bitset<30>(ix) << " " << std::bitset<30>(iy) << " " << std::bitset<30>(iz) << std::endl;
  std::uint32_t ixyz = ix * 4 + iy * 2 + iz;
  return ixyz;
}
#ifdef DFM2_STATIC_LIBRARY
template std::uint32_t delfem2::MortonCode(float x, float y, float z);
template std::uint32_t delfem2::MortonCode(double x, double y, double z);
#endif



DFM2_INLINE std::pair<unsigned int,unsigned int> delfem2::MortonCode_DeterminRange(
    const std::uint32_t* sortedMC,
    size_t nMC,
    unsigned int imc)
{
  assert( nMC > 0 );
  if( imc == 0 ){ 
	  return std::make_pair(
		  0,
		  static_cast<unsigned int>(nMC-1) ); 
  }
  if( imc == nMC-1 ){ 
	  return std::make_pair(
		  static_cast<unsigned int>(nMC-1),
		  static_cast<unsigned int>(nMC-1) ); 
  }
  // ----------------------
  const std::uint32_t mc0 = sortedMC[imc-1];
  const std::uint32_t mc1 = sortedMC[imc+0];
  const std::uint32_t mc2 = sortedMC[imc+1];
  if( mc0 == mc1 && mc1 == mc2 ){ // for hash value collision
    unsigned int jmc=imc+1;
    for(;jmc<nMC;++jmc){
      if( sortedMC[jmc] != mc1 ) break;
    }
    return std::make_pair(imc,jmc-1);
  }
  // get direction
  // (d==+1) -> imc is left-end, move forward
  // (d==-1) -> imc is right-end, move backward
  int d = bvh::delta(imc, imc + 1, sortedMC, nMC) - bvh::delta(imc, imc - 1, sortedMC, nMC);
  d = d > 0 ? 1 : -1;
  
  //compute the upper bound for the length of the range
  const int delta_min = bvh::delta(imc, imc - d, sortedMC, nMC);
  int lmax = 2;
  while ( bvh::delta(imc, imc + lmax*d, sortedMC, nMC)>delta_min )
  {
    lmax = lmax * 2;
  }
  
  //find the other end using binary search
  int l = 0;
  for (int t = lmax / 2; t >= 1; t /= 2)
  {
    if (bvh::delta(imc, imc + (l + t)*d, sortedMC, nMC)>delta_min)
    {
      l = l + t;
    }
  }
  unsigned int j = imc + l*d;
  
  std::pair<unsigned int, unsigned int> range;
  if (imc <= j) { range.first = imc; range.second = j; }
  else { range.first = j; range.second = imc; }
  return range;
}

DFM2_INLINE unsigned int delfem2::MortonCode_FindSplit(
    const std::uint32_t* aMC,
    unsigned int iMC_start,
    unsigned int iMC_last)
{
  if (iMC_start == iMC_last) { return UINT_MAX; }

  const std::uint32_t mcStart = aMC[iMC_start];

  const unsigned int nbitcommon0 = nbits_leading_zero(mcStart ^ aMC[iMC_last]);
  
  // handle duplicated morton code
  if ( nbitcommon0 == 32 ){ return iMC_start; } // sizeof(std::uint32_t)*8
  
  // Use binary search to find where the next bit differs.
  // Specifically, we are looking for the highest object that
  // shares more than commonPrefix bits with the first one.
  unsigned int iMC_split = iMC_start; // initial guess
  unsigned int step = iMC_last - iMC_start;
  while (step > 1){
    step = (step + 1) / 2; // half step
    const unsigned int iMC_new = iMC_split + step; // proposed new position
    if ( iMC_new >= iMC_last ) { continue; }
    const unsigned int nbitcommon1 = nbits_leading_zero(mcStart ^ aMC[iMC_new] );
    if ( nbitcommon1 > nbitcommon0 ){
      iMC_split = iMC_new; // accept proposal
    }
  }
  return iMC_split;
}

template <typename REAL>
DFM2_INLINE void delfem2::SortedMortenCode_Points3(
    std::vector<unsigned int> &aSortedId,
    std::vector<std::uint32_t> &aSortedMc,
    const std::vector<REAL> &aXYZ,
    const REAL min_xyz[3],
    const REAL max_xyz[3])
{
  std::vector<bvh::CPairMtcInd> aNodeBVH0; // array of BVH node
  const std::size_t np = aXYZ.size()/3;
  aNodeBVH0.resize(np);
  const REAL x_min = min_xyz[0];
  const REAL y_min = min_xyz[1];
  const REAL z_min = min_xyz[2];
  const REAL x_max = max_xyz[0];
  const REAL y_max = max_xyz[1];
  const REAL z_max = max_xyz[2];
  for(unsigned ip=0;ip<np;++ip){
    const REAL x = (aXYZ[ip*3+0]-x_min)/(x_max-x_min);
    const REAL y = (aXYZ[ip*3+1]-y_min)/(y_max-y_min);
    const REAL z = (aXYZ[ip*3+2]-z_min)/(z_max-z_min);
    aNodeBVH0[ip].imtc = MortonCode(x,y,z);
    aNodeBVH0[ip].iobj = ip;
  }
  std::sort(aNodeBVH0.begin(), aNodeBVH0.end());
  aSortedId.resize(aNodeBVH0.size());
  aSortedMc.resize(aNodeBVH0.size());
  for(size_t ino=0;ino<aNodeBVH0.size();++ino){
    aSortedMc[ino] = aNodeBVH0[ino].imtc;
    aSortedId[ino] = aNodeBVH0[ino].iobj;
      //        std::cout << std::bitset<32>(aNodeBVH[ino].imtc) << "  " << clz(aNodeBVH[ino].imtc) << "   " << ino << std::endl;
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::SortedMortenCode_Points3(
    std::vector<unsigned int>& aSortedId,
    std::vector<std::uint32_t>& aSortedMc,
    const std::vector<float>& aXYZ,
    const float min_xyz[3],
    const float max_xyz[3]);
template void delfem2::SortedMortenCode_Points3(
    std::vector<unsigned int>& aSortedId,
    std::vector<std::uint32_t>& aSortedMc,
    const std::vector<double>& aXYZ,
    const double min_xyz[3],
    const double max_xyz[3]);
#endif

// ----------------------------------

DFM2_INLINE void delfem2::BVHTopology_Morton(
    std::vector<CNodeBVH2>& aNodeBVH,
    const std::vector<unsigned int>& aSortedId,
    const std::vector<std::uint32_t>& aSortedMc)
{
  assert( aSortedId.size() == aSortedMc.size() );
  assert( !aSortedMc.empty() );
  aNodeBVH.resize(aSortedMc.size()*2-1);
  aNodeBVH[0].iparent = UINT_MAX;
  const unsigned int nni = static_cast<unsigned int>(aSortedMc.size()-1); // number of internal node
  for(unsigned int ini=0;ini<nni;++ini){
    const std::pair<unsigned int, unsigned int> range = MortonCode_DeterminRange(aSortedMc.data(), aSortedMc.size(), ini);
    unsigned int isplit = MortonCode_FindSplit(aSortedMc.data(), range.first, range.second);
    assert( isplit != UINT_MAX );
    if( range.first == isplit ){
      const unsigned int inlA = nni+isplit;
      aNodeBVH[ini].ichild[0] = inlA;
      aNodeBVH[inlA].iparent = ini;
      aNodeBVH[inlA].ichild[0] = aSortedId[isplit];
      aNodeBVH[inlA].ichild[1] = UINT_MAX;
    }
    else{
      const unsigned int iniA = isplit;
      aNodeBVH[ini].ichild[0] = iniA;
      aNodeBVH[iniA].iparent = ini;
    }
    // ----
    if( range.second == isplit+1 ){
      const unsigned int inlB = nni+isplit+1;
      aNodeBVH[ini].ichild[1] = inlB;
      aNodeBVH[inlB].iparent = ini;
      aNodeBVH[inlB].ichild[0] = aSortedId[isplit+1];
      aNodeBVH[inlB].ichild[1] = UINT_MAX;
    }
    else{
      const unsigned int iniB = isplit+1;
      aNodeBVH[ini].ichild[1] = iniB;
      aNodeBVH[iniB].iparent = ini;
    }
  }
}


DFM2_INLINE void delfem2::Check_MortonCode_Sort(
  [[maybe_unused]] const std::vector<unsigned int>& aSortedId,
  [[maybe_unused]] const std::vector<std::uint32_t>& aSortedMc,
  [[maybe_unused]] const std::vector<double>& aXYZ,
  [[maybe_unused]] const double bbmin[3],
  [[maybe_unused]] const double bbmax[3])
{
#ifdef NDEBUG
  return;
#else
  for(unsigned int imc=1;imc<aSortedMc.size();++imc){
    std::uint32_t mc0 = aSortedMc[imc-1];
    std::uint32_t mc1 = aSortedMc[imc+0];
    assert( mc0 <= mc1 );
  }
  for(unsigned int imc=0;imc<aSortedMc.size();++imc){
    std::uint32_t mc0 = aSortedMc[imc];
    unsigned int ip = aSortedId[imc];
    double x0 = aXYZ[ip*3+0];
    double y0 = aXYZ[ip*3+1];
    double z0 = aXYZ[ip*3+2];
    double x1 = (x0-bbmin[0])/(bbmax[0]-bbmin[0]);
    double y1 = (y0-bbmin[1])/(bbmax[1]-bbmin[1]);
    double z1 = (z0-bbmin[2])/(bbmax[2]-bbmin[2]);
    std::uint32_t mc1 = MortonCode(x1,y1,z1);
    assert( mc0 == mc1 );
  }
  /*
  for(int ip=0;ip<aSortedId.size();++ip){
    std::uint32_t mc0 = aSortedMc[ip];
    std::cout << std::bitset<32>(mc0) << " " << ip << " " << mc0 << std::endl;
  }
   */
#endif
}

DFM2_INLINE void delfem2::Check_MortonCode_RangeSplit(
    [[maybe_unused]] const std::vector<std::uint32_t>& aSortedMc)
{
#ifdef NDEBUG
  return;
#else
  assert(!aSortedMc.empty());
  for(unsigned int ini=0;ini<aSortedMc.size()-1;++ini){
    const std::pair<unsigned int,unsigned int> range = MortonCode_DeterminRange(aSortedMc.data(), aSortedMc.size(), ini);
    const unsigned int isplit = MortonCode_FindSplit(aSortedMc.data(), range.first, range.second);
    const std::pair<unsigned int,unsigned int> rangeA = MortonCode_DeterminRange(aSortedMc.data(), aSortedMc.size(), isplit);
    const std::pair<unsigned int,unsigned int> rangeB = MortonCode_DeterminRange(aSortedMc.data(), aSortedMc.size(), isplit+1);
    assert( range.first == rangeA.first );
    assert( range.second == rangeB.second );
    {
      const unsigned int last1 = ( isplit == range.first ) ? isplit : rangeA.second;
      const unsigned int first1 = ( isplit+1 == range.second ) ? isplit+1 : rangeB.first;
      assert( last1+1 == first1 );
    }
  }
#endif
}

DFM2_INLINE void delfem2::Check_BVH(
	const std::vector<CNodeBVH2>& aNodeBVH,
	size_t N)
{
  std::vector<int> aFlg(N,0);
  bvh::mark_child(aFlg, 0, aNodeBVH);
  for(unsigned int i=0;i<N;++i){
    assert(aFlg[i]==1);
  }
}
