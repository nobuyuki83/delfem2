
#include "delfem2/matb.h"

static inline void CalcInvMat3(double a[], double t[] ){
  const double det = a[0]*a[4]*a[8] + a[3]*a[7]*a[2] + a[6]*a[1]*a[5]
  - a[0]*a[7]*a[5] - a[6]*a[4]*a[2] - a[3]*a[1]*a[8];
  const double inv_det = 1.0/det;
  for(int i=0;i<9;i++){ t[i] = a[i]; }
  a[0] = inv_det*(t[4]*t[8]-t[5]*t[7]);
  a[1] = inv_det*(t[2]*t[7]-t[1]*t[8]);
  a[2] = inv_det*(t[1]*t[5]-t[2]*t[4]);
  a[3] = inv_det*(t[5]*t[6]-t[3]*t[8]);
  a[4] = inv_det*(t[0]*t[8]-t[2]*t[6]);
  a[5] = inv_det*(t[2]*t[3]-t[0]*t[5]);
  a[6] = inv_det*(t[3]*t[7]-t[4]*t[6]);
  a[7] = inv_det*(t[1]*t[6]-t[0]*t[7]);
  a[8] = inv_det*(t[0]*t[4]-t[1]*t[3]);
}

// define fixed boudnary condition
void CTriDiaMat3::FixBC(unsigned int ino, unsigned int idof){
  assert( idof < 3 && ino < n );
  if( ino != 0 ){
    double* pvu = v+ino*27-18;
    double* pvl = v+ino*27-9;
    pvu[0*3+idof] = 0;	pvu[1*3+idof] = 0;	pvu[2*3+idof] = 0;
    pvl[idof*3+0] = 0;	pvl[idof*3+1] = 0;	pvl[idof*3+2] = 0;
  }
  if( ino != n-1 ){
    double* pvu = v+ino*27+18;
    double* pvl = v+ino*27+9;
    pvu[0*3+idof] = 0;	pvu[1*3+idof] = 0;	pvu[2*3+idof] = 0;
    pvl[idof*3+0] = 0;	pvl[idof*3+1] = 0;	pvl[idof*3+2] = 0;
  }
  double* pvc = v+ino*27;
  pvc[0*3+idof] = 0;	pvc[1*3+idof] = 0;	pvc[2*3+idof] = 0;
  pvc[idof*3+0] = 0;	pvc[idof*3+1] = 0;	pvc[idof*3+2] = 0;
  pvc[idof*3+idof] = 1;
}
// execute ILU factorization
void CTriDiaMat3::ILU_Frac()
{
  double tmpBlk[9];
  for(unsigned int iblk=0;iblk<n;iblk++){
    if( iblk != 0 ){
      const double* pVal_ik = v+27*iblk-9;
      const double* pVal_kj = v+27*iblk-18;
      double* pVal_ij = v+27*iblk;
      for(unsigned int i=0;i<3;i++){
        pVal_ij[i*3+0] -= pVal_ik[i*3+0]*pVal_kj[0] + pVal_ik[i*3+1]*pVal_kj[3] + pVal_ik[i*3+2]*pVal_kj[6];
        pVal_ij[i*3+1] -= pVal_ik[i*3+0]*pVal_kj[1] + pVal_ik[i*3+1]*pVal_kj[4] + pVal_ik[i*3+2]*pVal_kj[7];
        pVal_ij[i*3+2] -= pVal_ik[i*3+0]*pVal_kj[2] + pVal_ik[i*3+1]*pVal_kj[5] + pVal_ik[i*3+2]*pVal_kj[8];
      }
    }
    {   // calc inverse of diagonal
      double* pVal_ii = v+27*iblk;
      CalcInvMat3(pVal_ii,tmpBlk);
    }
    // [U] = [1/D][U]
    if( iblk !=  n-1 ){
      double* pVal_ij = v+27*iblk+9;
      const double* pVal_ii = v+27*iblk;
      for(unsigned int i=0;i<9;i++){ tmpBlk[i] = pVal_ij[i]; }
      for(unsigned int i=0;i<3;i++){
        pVal_ij[i*3+0] = pVal_ii[i*3+0]*tmpBlk[0] + pVal_ii[i*3+1]*tmpBlk[3] + pVal_ii[i*3+2]*tmpBlk[6];
        pVal_ij[i*3+1] = pVal_ii[i*3+0]*tmpBlk[1] + pVal_ii[i*3+1]*tmpBlk[4] + pVal_ii[i*3+2]*tmpBlk[7];
        pVal_ij[i*3+2] = pVal_ii[i*3+0]*tmpBlk[2] + pVal_ii[i*3+1]*tmpBlk[5] + pVal_ii[i*3+2]*tmpBlk[8];
      }
    }
  }	// end iblk
}

// solve matrix
void CTriDiaMat3::Solve(std::vector<double>& res){
  double pTmpVec[3];
  for(unsigned int iblk=0;iblk<n;iblk++){
    pTmpVec[0] = res[iblk*3+0];
    pTmpVec[1] = res[iblk*3+1];
    pTmpVec[2] = res[iblk*3+2];
    if( iblk != 0 ){
      const double* pVal_ij = v+iblk*27-9;
      const double valj0 = res[(iblk-1)*3+0];
      const double valj1 = res[(iblk-1)*3+1];
      const double valj2 = res[(iblk-1)*3+2];
      pTmpVec[0] -= pVal_ij[0]*valj0+pVal_ij[1]*valj1+pVal_ij[2]*valj2;
      pTmpVec[1] -= pVal_ij[3]*valj0+pVal_ij[4]*valj1+pVal_ij[5]*valj2;
      pTmpVec[2] -= pVal_ij[6]*valj0+pVal_ij[7]*valj1+pVal_ij[8]*valj2;
    }
    const double* pVal_ii = v+27*iblk;
    res[iblk*3+0] = pVal_ii[0]*pTmpVec[0]+pVal_ii[1]*pTmpVec[1]+pVal_ii[2]*pTmpVec[2];
    res[iblk*3+1] = pVal_ii[3]*pTmpVec[0]+pVal_ii[4]*pTmpVec[1]+pVal_ii[5]*pTmpVec[2];
    res[iblk*3+2] = pVal_ii[6]*pTmpVec[0]+pVal_ii[7]*pTmpVec[1]+pVal_ii[8]*pTmpVec[2];
  }
  for(int iblk=n-1;iblk>=0;iblk--){
    pTmpVec[0] = res[iblk*3+0];
    pTmpVec[1] = res[iblk*3+1];
    pTmpVec[2] = res[iblk*3+2];
    if( iblk != (int)n-1 ){
      const double* pVal_ij = v+27*iblk+9;
      const double valj0 = res[(iblk+1)*3+0];
      const double valj1 = res[(iblk+1)*3+1];
      const double valj2 = res[(iblk+1)*3+2];
      pTmpVec[0] -= pVal_ij[0]*valj0+pVal_ij[1]*valj1+pVal_ij[2]*valj2;
      pTmpVec[1] -= pVal_ij[3]*valj0+pVal_ij[4]*valj1+pVal_ij[5]*valj2;
      pTmpVec[2] -= pVal_ij[6]*valj0+pVal_ij[7]*valj1+pVal_ij[8]*valj2;
    }
    res[iblk*3+0] = pTmpVec[0];
    res[iblk*3+1] = pTmpVec[1];
    res[iblk*3+2] = pTmpVec[2];
  }
}
