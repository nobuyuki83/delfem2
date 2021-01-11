#include "thrust/host_vector.h"
#include "thrust/device_vector.h"


__global__ void kernel(int *out, int *in, const int n)
{
  unsigned int i = threadIdx.x;
  if (i < n) {
    out[i] = in[i] * 2;
  }
}

int main(){
  thrust::host_vector<int> hVectorIn(1024);
  for(int i=0;i<1024;i++){
    hVectorIn[i]=i;
  }

  thrust::device_vector<int> dVectorIn(1024);
  thrust::device_vector<int> dVectorOut(1024);
  dVectorIn=hVectorIn;

  kernel<<<1,1024>>>(
      thrust::raw_pointer_cast(dVectorOut.data()),
      thrust::raw_pointer_cast(dVectorIn.data()),
      1024);

  thrust::host_vector<int> hVectorOut = dVectorOut;
  for(int i=0;i<1024;i++){
    std::cout << i << " " << hVectorOut[i] << " " << hVectorIn[i] << std::endl;
  }
}