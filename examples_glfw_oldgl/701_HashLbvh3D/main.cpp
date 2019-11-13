#include <iostream>
#include <vector>
#include <bitset>
#include <stdint.h> // std::uint32_t
#include "delfem2/bv.h"
#include "delfem2/bvh.h"
#include "delfem2/primitive.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshtopo.h"
#include "delfem2/srchbi_v3bvh.h"


// ---------------------------------

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/gl2_funcs.h"

namespace dfm2 = delfem2;

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

std::uint32_t MortonCode(double x, double y, double z){
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

class CNodeBVH2_LBV{
public:
  std::uint32_t imc;
  std::uint32_t iobj;
  std::uint32_t ino_first;
  std::uint32_t ino_split;
  std::uint32_t ino_last;
public:
};


bool comp0(const CNodeBVH2_LBV& lhs, const CNodeBVH2_LBV& rhs){
  return lhs.imc < rhs.imc;
}

int findSplit
 (const std::vector<CNodeBVH2_LBV>& aNode,
  int inode)
{
  const CNodeBVH2_LBV& no = aNode[inode];
  unsigned int fist = no.ino_first;
  unsigned int last = no.ino_last;
  unsigned int firstCode = aNode[fist].imc;
  unsigned int lastCode = aNode[last].imc;
  
  if (firstCode == lastCode)
    return (fist + last) >> 1;
  
  const int commonPrefix = clz(firstCode ^ lastCode);
  
  int split = fist; // initial guess
  int step = last - fist;
  
  do
  {
    step = (step + 1) >> 1; // exponential decrease
    int newSplit = split + step; // proposed new position
    
    if (newSplit < last) {
      unsigned int splitCode = aNode[newSplit].imc;
      int splitPrefix = clz(firstCode ^ splitCode);
      if (splitPrefix > commonPrefix)
        split = newSplit; // accept proposal
    }
  }
  while (step > 1);
  
  return split;
}



// ------------------------------------
// input parameter for simulation
std::vector<double> aXYZ; // 3d points
std::vector<CNodeBVH2_LBV> aNodeBVH; // array of BVH node

// data for camera
double cur_time = 0;

// ----------------------------------------

void myGlutDisplay(void)
{
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  ::glPointSize(3);
  ::glBegin(GL_POINTS);
  for(int ip=0;ip<aXYZ.size()/3;++ip){
    ::glVertex3d(aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]);
  }
  ::glEnd();
}

int main(int argc,char* argv[])
{
  {
    dfm2::CBV3D_AABB bb(-1,+1, -1,+1, -1,+1);
    const unsigned int N = 1000;
    aXYZ.resize(N*3);
    for(int i=0;i<N;++i){
      aXYZ[i*3+0] = (bb.x_max -  bb.x_min) * rand()/(RAND_MAX+1.0) + bb.x_min;
      aXYZ[i*3+1] = (bb.y_max -  bb.y_min) * rand()/(RAND_MAX+1.0) + bb.y_min;
      aXYZ[i*3+2] = (bb.z_max -  bb.z_min) * rand()/(RAND_MAX+1.0) + bb.z_min;
    }
    aNodeBVH.resize(N);
    for(int ip=0;ip<aXYZ.size()/3;++ip){
      double x = (aXYZ[ip*3+0]-bb.x_min)/(bb.x_max-bb.x_min);
      double y = (aXYZ[ip*3+1]-bb.y_min)/(bb.y_max-bb.y_min);
      double z = (aXYZ[ip*3+2]-bb.z_min)/(bb.z_max-bb.z_min);
      aNodeBVH[ip].imc = MortonCode(x,y, z);
      aNodeBVH[ip].iobj = ip;
    }
    std::sort(aNodeBVH.begin(), aNodeBVH.end(), comp0);
    for(int ino=0;ino<aNodeBVH.size();++ino){
      std::cout << std::bitset<32>(aNodeBVH[ino].imc) << "  " << clz(aNodeBVH[ino].imc) << "   " << ino << std::endl;
    }
    aNodeBVH[0].ino_first = 0;
    aNodeBVH[0].ino_last = aNodeBVH.size()-1;
    std::cout << "split " << findSplit(aNodeBVH, 0) << std::endl;
  }
  
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.5;
  delfem2::opengl::setSomeLighting();
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


