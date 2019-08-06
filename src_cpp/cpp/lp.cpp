/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cassert>

#include "delfem2/lp.h"


void Print(const std::vector<double>& A, int ncol, int nrow,
           const std::vector<int>& map_col2row)
{
  for(int icol=0;icol<ncol;++icol){
    std::cout << icol << " " << map_col2row[icol] << " --> ";
    for(int irow=0;irow<nrow;++irow){
      std::cout << A[icol*nrow+irow] << "   ";
    }
    std::cout << std::endl;
  }
}

void Print(const std::vector<double>& A, int ncol, int nrow)
{
  for(int icol=0;icol<ncol;++icol){
    std::cout << icol << " --> ";
    for(int irow=0;irow<nrow;++irow){
      std::cout << A[icol*nrow+irow] << "   ";
    }
    std::cout << std::endl;
  }
}


// return value
// 0: converged
// 1: input value wrong
// 2: no bound
int LinPro_SolveTable
(int& nitr,
 std::vector<double>& B,
 std::vector<int>& map_col2row,
 int ncol, int nrow)
{
#ifndef NDEBUG
  {
    assert(B.size()==ncol*nrow);
    assert(map_col2row.size()==ncol);
    std::vector<int> aFlgRow(nrow,0);
    for(int jcol=0;jcol<ncol;++jcol){
      const int jrow = map_col2row[jcol];
      assert( jrow >= 0 && jrow < nrow );
      assert(aFlgRow[jrow] == 0 );
      aFlgRow[jrow] = 1;
    }
  }
#endif
  
  for(int icol=0;icol<ncol-1;++icol){
    int jrow = map_col2row[icol];
    double b0 = B[icol*nrow+jrow];
    double b1 = B[(ncol-1)*nrow+jrow];
    if( fabs(b1)<1.0e-30 ) continue;
    assert( fabs(b0) > 1.0e-30 );
    double ratio = b1/b0;
    for(int jrow=0;jrow<nrow;++jrow){
      B[(ncol-1)*nrow+jrow] -= ratio*B[icol*nrow+jrow];
    }
  }
//  ::Print(B, ncol, nrow, map_col2row);
  ////
  std::vector<int> flg_row(nrow,0); // 0:base 1:non_base 2:trg
  for(int ieq=0;ieq<ncol;++ieq){
    int jrow = map_col2row[ieq];
    flg_row[jrow] = 1;
  }
  flg_row[0] = 2;
  ////
#ifndef NDEBUG
  for(int jcol=0;jcol<ncol;++jcol){
    int jrow = map_col2row[jcol];
    for(int icol=0;icol<ncol;++icol){
      double b = B[icol*nrow+jrow];
      if( icol == jcol ){
        assert(fabs(b)>1.0e-30); // typically +1 or -1
        double b2 = B[jcol*nrow]; // the origin (0,0,...) should be valid
        if( jcol != ncol-1 ){
          if( b*b2 < 0 ){
            std::cout << "wrong entry" << icol << " " << jrow << " " << b2 << " " << b << std::endl;
            return 1;
          }
          assert(b*b2>=0); // the initial target value can be negtive
        }
      }
      else{
        assert(fabs(b)<1.0e-30);
      }
    }
  }
#endif
  ////
  for(int itr=0;itr<nitr;++itr){
//    std::cout << "iteration: " << itr << std::endl;
//    ::Print(B,ncol,nrow,map_col2row);
    int icol_min = -1, jrow_min = -1;
    { // find minimum row
      bool is_optim = true;
      for(int jrow=0;jrow<nrow;++jrow){
        if( flg_row[jrow] != 0 ){ continue; }
        double d0 = B[(ncol-1)*nrow+jrow];
        if( d0 >= 0 ){ continue; }
        is_optim = false;
        double min_bound = -1;
        for(int icol=0;icol<ncol-1;++icol){
          double v0 = B[icol*nrow+jrow];
          double v1 = B[icol*nrow];
//          std::cout << "       hogehoge" << v0 << " " << v1 << std::endl;
          if( v0*v1<=0 ){ continue; }
          double bound0 = v1/v0;
//          std::cout << "        val improv" << icol << " " << jrow << " " << bound0 << std::endl;
          if( min_bound < 0 || bound0 < min_bound ){
            min_bound = bound0;
            jrow_min = jrow;
            icol_min = icol;
          }
        }
        if( min_bound < 0 ){ // no bound for improvemnet
          continue;
        }
        else{
          assert( jrow_min != -1 && icol_min != -1 );
          break;
        }
      }
//      std::cout << "mininum index " << icol_min << " " << jrow_min << std::endl;
      if( is_optim ){ // no way to furthre increase the value
        nitr = itr;
        return 0; // converged
      }
    }
    if( icol_min == -1 || jrow_min == -1 ){
      nitr = itr;
      return 0;
    }
//    std::cout << "minimum bottom" << itr << " " << B[(ncol-1)*nrow+jrow_min] << "     " << B[(ncol-1)*nrow]  << " " << jrow_min << " " << std::endl;

    assert( icol_min >= 0 );
//    std::cout << "itr: " << itr << " " << icol_min << " " << jrow_min << " " << B[icol_min*nrow+jrow_min] << std::endl;
    if( fabs(B[icol_min*nrow+jrow_min])<1.0e-30 ){
//      std::cout << "invalid: " << itr << " " << B[icol_min*nrow+jrow_min] << std::endl;
      return 2;
    }
//    std::cout << "icolmin" << icol_min << " " << jrow_min << std::endl;
    { // Gauss's sweep
      const double vpiv = B[icol_min*nrow+jrow_min];
      for(int jrow=0;jrow<nrow;++jrow){ B[icol_min*nrow+jrow] /= vpiv; }
      for(int icol=0;icol<ncol;++icol){
        if( icol == icol_min ) continue;
        const double vpiv = B[icol*nrow+jrow_min];
        for(int jrow=0;jrow<nrow;jrow++){
          B[icol*nrow+jrow] -= vpiv*B[icol_min*nrow+jrow];
        }
      }
    }
    // swap basic-nonbasic index
    flg_row[jrow_min] = 1;
    int jrow_new = map_col2row[icol_min];
    flg_row[jrow_new] = 0;
    map_col2row[icol_min] = jrow_min;
  }
  return false;
}




////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////


bool CLinPro::Solve
(std::vector<double>& solution,
 double& opt_val,
 int& nitr,
 const std::vector<double>& aCoeffTrg) const
{
  const int ncol = (neq+1);
  const int nrow = (1+1+nvar+nslk);
  assert(this->A.size()==ncol*nrow);
  std::vector<int> map_col2rowB = map_col2row;
  std::vector<double> B = this->A;
  B[(ncol-1)*nrow+1] = 1;
  for(int ic=0;ic<aCoeffTrg.size();++ic){
    B[(ncol-1)*nrow+2+ic] = -aCoeffTrg[ic];
  }
//  std::cout << "before solve" << std::endl;
// ::Print(B, ncol, nrow, map_col2rowB);
  LinPro_SolveTable(nitr, B, map_col2rowB,
                    ncol, nrow);
//  std::cout << "after solve" << std::endl;
//  ::Print(B, ncol, nrow, map_col2rowB);
  std::vector<double> buff(nrow,0.0);
  for(int icol=0;icol<ncol;++icol){
    int jrow = map_col2rowB[icol];
    buff[jrow] = B[icol*nrow];
  }
  solution.resize(nvar,0.0);
  for(int ivar=0;ivar<nvar;++ivar){
    solution[ivar] = buff[ivar+2];
//    std::cout << " sol" << ivar << " " << solution[ivar] << std::endl;
  }
  opt_val = buff[1];
//  std::cout << " opt" << opt_val << std::endl;
  return 0;
}


void CLinPro::AddEqn
(const std::vector<double>& aW,
 double rhs,
 EQ_TYPE type)
{
  CEq eq;
  eq.aCoeff.assign(aW.begin(),aW.end());
  eq.rhs = rhs;
  eq.itype = type;
  aEq.push_back(eq);
}


// 0 converged
// 2 nobound
int CLinPro::Precomp(int& nitr)
{
  neq = aEq.size();
  nvar = 0;
  nslk = 0;
  nart = 0;
  std::vector<int> mapEq2Slk(neq,-1);
  std::vector<int> mapEq2Art(neq,-1);
  for(int ieq=0;ieq<aEq.size();++ieq){
    const int nv = aEq[ieq].aCoeff.size();
    if( nvar < nv ){ nvar = nv; }
    if( aEq[ieq].itype == EQ ){ mapEq2Slk[ieq] = nslk; nslk++; }
    if( aEq[ieq].itype == LE || aEq[ieq].itype == GE ){ mapEq2Slk[ieq] = nslk; nslk++; }
    if( aEq[ieq].itype == LE && aEq[ieq].rhs < 0  ){ mapEq2Art[ieq] = nart; nart++; }
    if( aEq[ieq].itype == GE && aEq[ieq].rhs > 0  ){ mapEq2Art[ieq] = nart; nart++; }
  }
  const int ncol = (neq+1); // neq,trg
  const int nrow = (1+1+nvar+nslk+nart); // rhs,trg,var,slk,art
//  std::cout << neq << " " << nvar << " " << nslk << " " << nart << "  " << ncol << " " << nrow << std::endl;
  A.resize(ncol*nrow,0.0);
  for(int ieq=0;ieq<neq;++ieq){
    A[ieq*nrow] = aEq[ieq].rhs;
    for(int ic=0;ic<aEq[ieq].aCoeff.size();++ic){
      A[ieq*nrow+2+ic] = aEq[ieq].aCoeff[ic];
    }
    if( aEq[ieq].itype == LE ){
      const int islk = mapEq2Slk[ieq]; assert(islk!=-1);
      A[ieq*nrow+2+nvar+islk] = +1;
    }
    if( aEq[ieq].itype == GE ){
      const int islk = mapEq2Slk[ieq]; assert(islk!=-1);
      A[ieq*nrow+2+nvar+islk] = -1;
    }
    if( aEq[ieq].itype == LE && aEq[ieq].rhs < 0 ){
      const int iart = mapEq2Art[ieq]; assert(iart!=-1);
      A[ieq*nrow+2+nvar+nslk+iart] = -1;
    }
    if( aEq[ieq].itype == GE && aEq[ieq].rhs > 0 ){
      const int iart = mapEq2Art[ieq]; assert(iart!=-1);
      A[ieq*nrow+2+nvar+nslk+iart] = +1;
    }
    if( aEq[ieq].itype == EQ ){
      const int islk = mapEq2Slk[ieq]; assert(islk!=-1);
       if( aEq[ieq].rhs >= 0 ){
         A[ieq*nrow+2+nvar+islk] = +1;
       }
       else{
         A[ieq*nrow+2+nvar+islk] = -1;
       }
    }
  }
  map_col2row.assign(ncol,0); // 0:base 1:non_base 2:trg
  map_col2row[ncol-1] = 1;
  for(int ieq=0;ieq<neq;ieq++){
    map_col2row[ieq] = 2+nvar+ieq;
  }
//  std::cout << "precomp" << std::endl;
//  ::Print(A, ncol, nrow, map_col2row);
  if( nart == 0 ){ nitr=0; return 0; }
  ////
  A[(ncol-1)*nrow+1] = 1.0;
  for(int ieq=0;ieq<aEq.size();++ieq){
    if( mapEq2Art[ieq] == -1 ) continue;
    int iart = mapEq2Art[ieq];
    int jrow1 = 1+1+nvar+nslk+iart;
    A[(ncol-1)*nrow+jrow1] = 1;
//    std::cout << "new entry: " << ieq << " " << jrow1 << std::endl;
    /*
    if( A[ieq*nrow] < 0 ){
      for(int jrow=0;jrow<nrow;++jrow){
        A[(ncol-1)*nrow+jrow] += A[ieq*nrow+jrow];
      }
    }
    else{
      for(int jrow=0;jrow<nrow;++jrow){
        A[(ncol-1)*nrow+jrow] -= A[ieq*nrow+jrow];
      }
    }
     */
    map_col2row[ieq] = jrow1;
  }
//  std::cout << "  precomp art" << std::endl;
//  ::Print(A, ncol, nrow, map_col2row);
  int res = LinPro_SolveTable(nitr, A, map_col2row, ncol,nrow);
  if( res == 2 ){
    std::cout << "no bound in solution" << std::endl;
    return 2;
  }
//  std::cout << "  solved" << std::endl;
//  ::Print(A, ncol, nrow, map_col2row);
  if( fabs(A[(ncol-1)*nrow])> 1.0e-8 ){
    std::cout << "couldn't found solution" << std::endl;
    return 3;
  }
//  std::cout << "succesfully found a valid solution" << std::endl;
  const std::vector<double> B = A;
  const int nrow1 = 1+1+nvar+nslk;
//  std::cout << "hugahugat" << ncol << " " << nrow1 << std::endl;
  A.resize(ncol*nrow1);
  for(int icol=0;icol<ncol;++icol){
    for(int irow=0;irow<nrow1;++irow){
      A[icol*nrow1+irow] = B[icol*nrow+irow];
    }
  }
  return 0;
}



