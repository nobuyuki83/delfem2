/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef LP_H
#define LP_H

void LinPro_Solve
(std::vector<double>& solution,
 int& nitr,
 int nvar, int nineq,
 std::vector<double>& A,
 std::vector<int>& flg_row,
 std::vector<int>& map_col2row)
{
  const int ncol = (nineq+1);
  const int nrow = (nvar+nineq+2);
  for(int itr=0;itr<nitr;++itr){
    int jrow_min = -1;
    double val_min = 0.0;
    {
      jrow_min = -1;
      for(int j=0;j<nrow;++j){
        if( flg_row[j] != 0 ){ continue; }
        double d0 = A[(ncol-1)*nrow+j];
        if( jrow_min == -1 || d0 < val_min ){
          jrow_min = j;
          val_min = d0;
          continue;
        }
      }
    }
    assert( jrow_min >=0 && jrow_min < nrow );
    if( val_min > 0 ){
      solution.assign(nrow,0.0);
      for(int icol=0;icol<ncol;++icol){
        int jrow = map_col2row[icol];
        solution[jrow]  = A[icol*nrow+(nrow-1)];
      }
      nitr = itr;
      break;
    }
    int icol_min = -1;
    { // detecting positive minimum bound
      double val_min = 0;
      for(int icol=0;icol<ncol-1;++icol){
        double v0 = A[icol*nrow+jrow_min];
        double v1 = A[icol*nrow+(nrow-1)];
        double v2 = v1/v0;
        if( icol_min == -1 || (v2>0 && v2<val_min ) ){
          val_min = v2;
          icol_min = icol;
        }
      }
    }
    //    std::cout << icol_min << " " << jrow_min << std::endl;
    {
      const double vpiv = A[icol_min*nrow+jrow_min];
      for(int jrow=0;jrow<nrow;++jrow){ A[icol_min*nrow+jrow] /= vpiv; }
    }
    for(int icol=0;icol<ncol;++icol){
      if( icol == icol_min ) continue;
      const double vpiv = A[icol*nrow+jrow_min];
      for(int jrow=0;jrow<nrow;jrow++){
        A[icol*nrow+jrow] -= vpiv*A[icol_min*nrow+jrow];
      }
    }
    /*
     for(int icol=0;icol<ncol;++icol){
     for(int jrow=0;jrow<nrow;++jrow){
     std::cout << "   " << icol << " " << jrow << " " <<  A[icol*nrow+jrow] << std::endl;
     }
     }
     */
    flg_row[jrow_min] = 1;
    int jrow_new = map_col2row[icol_min];
    flg_row[jrow_new] = 0;
    map_col2row[icol_min] = jrow_min;
  }
}


void LinPro_SetTarget
(std::vector<double>& A,
 std::vector<int>& flg_row,
 std::vector<int>& map_col2row,
 int nvar, int nineq,
 const std::vector<double>& aW)
{
  const int ncol = (nineq+1);
  const int nrow = (nvar+nineq+2);
  for(int iw=0;iw<aW.size();++iw){
    A[(ncol-1)*nrow+1+iw] = -aW[iw];
  }
}

void LinPro_SetIneqLe
(std::vector<double>& A,
 std::vector<int>& flg_row,
 std::vector<int>& map_col2row,
 int nvar, int nineq, int icol,
 const std::vector<double>& aW)
{
  const int ncol = (nineq+1);
  const int nrow = (nvar+nineq+2);
  assert( aW.size()-1 == nvar );
  for(int iw=0;iw<nvar;++iw){
    A[icol*nrow+1+iw] = aW[iw];
  }
  A[icol*nrow+(nrow-1)] = aW[aW.size()-1];
//  A[icol*nrow+1+nvar+icol] = 1.0; // will be -1 for Ge
}

void LinPro_Init
(std::vector<double>& A,
 std::vector<int>& flg_row,
 std::vector<int>& map_col2row,
 int nvar, int nineq)
{
  const int ncol = (nineq+1);
  const int nrow = (nvar+nineq+2);
  A.assign(ncol*nrow,0.0);
  for(int ieq=0;ieq<nineq;ieq++){
    A[ieq*nrow+1+nvar+ieq] = 1.0;
  }
  A[(ncol-1)*nrow+0] = 1.0;
  ///////
  flg_row.assign(nrow,0); // 0:base 1:non_base 2:trg
  flg_row[0] = 1;
  for(int ieq=0;ieq<nineq;ieq++){
    flg_row[1+nvar+ieq] = 1;
  }
  flg_row[nrow-1] = 2;
  ///
  map_col2row.assign(ncol,0); // 0:base 1:non_base 2:trg
  map_col2row[ncol-1] = 0;
  for(int ieq=0;ieq<nineq;ieq++){
    map_col2row[ieq] = 1+nvar+ieq;
  }
}


#endif
