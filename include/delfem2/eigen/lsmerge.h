

#include <Eigen/Core>

template <int nrow, int ncol, int ndimrow, int ndimcol, typename T>
bool Merge(
    Eigen::MatrixXd& A,
    const unsigned int* aIpRow,
    const unsigned int* aIpCol,
    const T emat[nrow][ncol][ndimrow][ndimcol],
    std::vector<unsigned int>& tmp_buffer)
{
  for(unsigned int in=0;in<nrow;++in){
    for(unsigned int jn=0;jn<ncol;++jn) {
      const unsigned int ip = aIpRow[in];
      const unsigned int jp = aIpCol[jn];
      for(unsigned int idim=0;idim<ndimrow;++idim){
        for(unsigned int jdim=0;jdim<ndimcol;++jdim) {
          A(ip*ndimrow+idim,jp*ndimcol+jdim) += emat[in][jn][idim][jdim];
        }
      }
    }
  }
  return true;
}