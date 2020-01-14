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

  float minmax[6];
  dfm2::cuda::cuda_MinMax_Point3D(minmax,
      aXYZ.data(), aXYZ.size()/3);
  printf("min xyz %f %f %f\n",minmax[0],minmax[1],minmax[2]);
  printf("max xyz %f %f %f\n",minmax[3],minmax[4],minmax[5]);

  clock_t time1 = clock();
  printf("on gpu it takes %.2f sec \n",(double)(time1-time0)/CLOCKS_PER_SEC);

}
