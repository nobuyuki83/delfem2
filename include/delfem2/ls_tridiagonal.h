#ifndef DFM2_LS_TRIDIAGONAL_H
#define DFM2_LS_TRIDIAGONAL_H


/**
 * 0 1 # # #
 * 2 3 4 # #
 * # 5 6 7 #
 * # # 8 9 10
 */
class BlockTriDiagonalMatrix {
 public:
  // initialize with block size n
  BlockTriDiagonalMatrix() { n_ = 0; }
  ~BlockTriDiagonalMatrix() {}

  void Initialize(size_t n) {
    this->n_ = n;
    v.resize((n_ * 3 - 2) * 9);
  }
  void Clear() {
    n_ = 0;
    v.clear();
  }

  //! clear value
  void SetZero() { v.assign((n_ * 3 - 2) * 9, 0.0); }

  //! marge element stiffness matrix to the position (idiv,idiv+1)
  void Merge(unsigned int idiv, double eM[][2][3][3]) {
    for (unsigned int i = 0; i < 36; i++) { v[idiv * 27 + i] += (&eM[0][0][0][0])[i]; }
  }

  //! define fixed boudnary condition
  void FixBC(unsigned int ino, unsigned int idof);

  // execute LU factorization
  void Decompose();

  // solve matrix
  void Solve(std::vector<double> &res);

 private:
  static inline void Inverse_Mat3(double a[], double t[]) {
    const double det =
        +a[0] * a[4] * a[8] + a[3] * a[7] * a[2] + a[6] * a[1] * a[5]
            - a[0] * a[7] * a[5] - a[6] * a[4] * a[2] - a[3] * a[1] * a[8];
    const double inv_det = 1.0 / det;
    for (int i = 0; i < 9; i++) { t[i] = a[i]; }
    a[0] = inv_det * (t[4] * t[8] - t[5] * t[7]);
    a[1] = inv_det * (t[2] * t[7] - t[1] * t[8]);
    a[2] = inv_det * (t[1] * t[5] - t[2] * t[4]);
    a[3] = inv_det * (t[5] * t[6] - t[3] * t[8]);
    a[4] = inv_det * (t[0] * t[8] - t[2] * t[6]);
    a[5] = inv_det * (t[2] * t[3] - t[0] * t[5]);
    a[6] = inv_det * (t[3] * t[7] - t[4] * t[6]);
    a[7] = inv_det * (t[1] * t[6] - t[0] * t[7]);
    a[8] = inv_det * (t[0] * t[4] - t[1] * t[3]);
  }
  
  void Sub_MatMat3(
      double pVal_ij[9],
      const double pVal_ik[9],
      const double pVal_kj[9]){
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; ++j) {
        pVal_ij[i * 3 + j] -=
            pVal_ik[i * 3 + 0] * pVal_kj[0 * 3 + j] +
                pVal_ik[i * 3 + 1] * pVal_kj[1 * 3 + j] +
                pVal_ik[i * 3 + 2] * pVal_kj[2 * 3 + j];
      }
    }
  }
  
  void MatMat3(double pVal_ij[9],
           const double pVal_ii[9],
           const double tmpBlk[9]){
    for (unsigned int i = 0; i < 3; i++) {
      for(unsigned int j=0;j<3;++j){
        pVal_ij[i * 3 + j]
        = pVal_ii[i * 3 + 0] * tmpBlk[0*3+j]
        + pVal_ii[i * 3 + 1] * tmpBlk[1*3+j]
        + pVal_ii[i * 3 + 2] * tmpBlk[2*3+j];
      }
    }
  }
  void Sub_MatVec3(double pTmpVec[3],
                   const double pVal_ij[9],
                   const double valj[3]){
    for(unsigned int i=0;i<3;++i){
      pTmpVec[i] -= pVal_ij[i*3+0] * valj[0] + pVal_ij[i*3+1] * valj[1] + pVal_ij[i*3+2] * valj[2];
    }
  }
  void MatVec3(double res[3],
               const double pVal_ii[9],
               const double pTmpVec[3]){
    for(unsigned int i=0;i<3;++i){
      res[i] = pVal_ii[i*3+0] * pTmpVec[0] + pVal_ii[i*3+1] * pTmpVec[1] + pVal_ii[i*3+2] * pTmpVec[2];
    }

  }
 private:
  unsigned int n_;
  std::vector<double> v;
};

// ------------------------------------------------------

// define fixed boudnary condition
void BlockTriDiagonalMatrix::FixBC(unsigned int ino, unsigned int idof) {
  assert(idof < 3 && ino < n_);
  if (ino != 0) {
    double *pvu = v.data() + ino * 27 - 18;
    double *pvl = v.data() + ino * 27 - 9;
    pvu[0 * 3 + idof] = 0;
    pvu[1 * 3 + idof] = 0;
    pvu[2 * 3 + idof] = 0;
    pvl[idof * 3 + 0] = 0;
    pvl[idof * 3 + 1] = 0;
    pvl[idof * 3 + 2] = 0;
  }
  if (ino != n_ - 1) {
    double *pvu = v.data() + ino * 27 + 18;
    double *pvl = v.data() + ino * 27 + 9;
    pvu[0 * 3 + idof] = 0;
    pvu[1 * 3 + idof] = 0;
    pvu[2 * 3 + idof] = 0;
    pvl[idof * 3 + 0] = 0;
    pvl[idof * 3 + 1] = 0;
    pvl[idof * 3 + 2] = 0;
  }
  double *pvc = v.data() + ino * 27;
  pvc[0 * 3 + idof] = 0;
  pvc[1 * 3 + idof] = 0;
  pvc[2 * 3 + idof] = 0;
  pvc[idof * 3 + 0] = 0;
  pvc[idof * 3 + 1] = 0;
  pvc[idof * 3 + 2] = 0;
  pvc[idof * 3 + idof] = 1;
}

// execute ILU factorization
void BlockTriDiagonalMatrix::Decompose() {
  double tmpBlk[9];
  constexpr int blksize = 9;
  for (unsigned int iblk = 0; iblk < n_; iblk++) {
    if (iblk != 0) {
      const double *pVal_ik = v.data() + blksize * 3 * iblk - blksize;
      const double *pVal_kj = v.data() + blksize * 3 * iblk - blksize * 2;
      double *pVal_ij = v.data() + blksize * 3 * iblk;
      Sub_MatMat3(pVal_ij,pVal_ik,pVal_kj);
    }
    {   // calc inverse of diagonal
      double *pVal_ii = v.data() + blksize * 3 * iblk;
      Inverse_Mat3(pVal_ii, tmpBlk);
    }
    // [U] = [1/D][U]
    if (iblk != n_ - 1) {
      double *pVal_ij = v.data() + blksize * 3 * iblk + blksize;
      const double *pVal_ii = v.data() + blksize * 3 * iblk;
      for (unsigned int i = 0; i < 9; i++) { tmpBlk[i] = pVal_ij[i]; }
      MatMat3(pVal_ij, pVal_ii, tmpBlk);
    }
  } // end iblk
}

// solve matrix
void BlockTriDiagonalMatrix::Solve(std::vector<double> &res) {
  double pTmpVec[3];
  for (unsigned int iblk = 0; iblk < n_; iblk++) {
    pTmpVec[0] = res[iblk * 3 + 0];
    pTmpVec[1] = res[iblk * 3 + 1];
    pTmpVec[2] = res[iblk * 3 + 2];
    if (iblk != 0) {
      const double *pVal_ij = v.data() + iblk * 27 - 9;
      const double *valj = res.data() + (iblk - 1) * 3;
      Sub_MatVec3(pTmpVec,pVal_ij,valj);
    }
    const double *pVal_ii = v.data() + 27 * iblk;
    MatVec3(res.data()+iblk*3,pVal_ii,pTmpVec);
  }
  for (int iblk = int(n_) - 1; iblk >= 0; iblk--) {
    pTmpVec[0] = res[iblk * 3 + 0];
    pTmpVec[1] = res[iblk * 3 + 1];
    pTmpVec[2] = res[iblk * 3 + 2];
    if (iblk != (int) n_ - 1) {
      const double *pVal_ij = v.data() + 27 * iblk + 9;
      const double *valj = res.data() + (iblk + 1) * 3;
      Sub_MatVec3(pTmpVec, pVal_ij, valj);
    }
    res[iblk * 3 + 0] = pTmpVec[0];
    res[iblk * 3 + 1] = pTmpVec[1];
    res[iblk * 3 + 2] = pTmpVec[2];
  }
}

#endif 
