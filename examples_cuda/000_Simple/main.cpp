#include <cstdio>
#include <cstdlib>
#include <vector>
#include "delfem2/cuda/cu_matvec.h"

int main()
{
  printf("Hello\n");

  const int n = 4000;
  std::vector<float> in(n);
  std::vector<float> out(n);
  std::vector<float> answer(n);

  float scale = 2.0;
  for (int i = 0; i < n; i++) in[i] = rand() % 100;
  for (int i = 0; i < n; i++) answer[i] = in[i] * scale;

  cuda_VecScale(out.data(), in.data(), scale, n);

  for (int i = 0; i < n; i++) {
      if (answer[i] != out[i]) {
          printf("error at index = %d\n", i);
          break;
      }
  }
  printf("OK\n");

  return 0;
}

