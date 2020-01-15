namespace delfem2 {
namespace cuda{

void cuda_VecScale(
    float *hOut,
    const float *hIn,
    float scale,
    const int n);

// -------------------------------------------------------------

float cuda_Dot(
    const float *h_A,
    const float *h_B,
    unsigned int n);

// -------------------------------------------------------------


void cuda_MatMat(
    float *h_C_gpu,
    const float *h_A,
    const float *h_B,
    unsigned int WIDTH);

// -------------------------------------------------------------\

void cuda_MinMax_Point3D(
    float *p_minmax,
    const float *p_XYZ,
    unsigned int np);

} // cuda
} // delfem2