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

// -------------------------------------------------------------

void cuda_MinMax_Point3D(
    float *p_minmax,
    const float *p_XYZ,
    unsigned int np);

// -------------------------------------------------------------

void cuda_CentRad_MeshTri3D(
    float* pXYZ_cent,
    float* pRad,
    const float *aXYZ,
    const unsigned int nXYZ,
    const unsigned int *aTri,
    const unsigned int nTri);

} // cuda
} // delfem2