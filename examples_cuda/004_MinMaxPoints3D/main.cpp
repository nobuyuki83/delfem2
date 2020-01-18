#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>

#include "delfem2/cuda/cu_matvec.h"

namespace dfm2 = delfem2;

int main() {
  std::vector<float> aXYZ;
  {
    const unsigned int np = 16*4000+1;
    aXYZ.resize(np*3);
    std::mt19937 engine(0);
    std::uniform_real_distribution<> dist1(+1.0, 2.0);
    for(int ip=0;ip<np;++ip) {
      aXYZ[ip * 3 + 0] = dist1(engine)*1.0;
      aXYZ[ip * 3 + 1] = dist1(engine)*0.5;
      aXYZ[ip * 3 + 2] = dist1(engine)*0.3;
    }
  }

  clock_t time0 = clock();

  float min3[3],max3[3];
  dfm2::cuda::cuda_Min3Max3_Points3F(min3,max3,
      aXYZ.data(), aXYZ.size()/3);
  printf("min xyz %f %f %f\n",min3[0],min3[1],min3[2]);
  printf("max xyz %f %f %f\n",max3[0],max3[1],max3[2]);

  clock_t time1 = clock();
  printf("on gpu it takes %.2f sec \n",(double)(time1-time0)/CLOCKS_PER_SEC);

}
