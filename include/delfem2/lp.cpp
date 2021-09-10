/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <iostream>
#include <cassert>
#include <map>
#include <climits>

#include "delfem2/lp.h"


DFM2_INLINE void Print(
    const std::vector<double>& A,
    unsigned int ncol,
    unsigned int nrow,
    const std::vector<unsigned int>& map_col2row)
{
  for(unsigned int icol=0;icol<ncol;++icol){
    std::cout << icol << " " << map_col2row[icol] << " --> ";
    for(unsigned int irow=0;irow<nrow;++irow){
      std::cout << A[icol*nrow+irow] << "   ";
    }
    std::cout << std::endl;
  }
}

DFM2_INLINE void Print(
    const std::vector<double>& A,
    int ncol,
    int nrow)
{
  for(int icol=0;icol<ncol;++icol){
    std::cout << icol << " --> ";
    for(int irow=0;irow<nrow;++irow){
      std::cout << A[icol*nrow+irow] << "   ";
    }
    std::cout << std::endl;
  }
}


bool LinPro_CheckTable(
    std::vector<double>& B,
    std::vector<unsigned int>& map_col2row,
    int ncol,
    int nrow)
{
  { // check if entries of map_col2row are unique
    assert(B.size()==ncol*nrow);
    assert(map_col2row.size()==ncol);
    std::vector<int> aFlgRow(nrow,0);
    for(int jcol=0;jcol<ncol;++jcol){
      const unsigned int jrow = map_col2row[jcol];
      assert( jrow >= 0 && jrow < nrow );
      assert(aFlgRow[jrow] == 0 );
      aFlgRow[jrow] = 1;
    }
  }
  for(int jcol=0;jcol<ncol;++jcol){
    unsigned int jrow_dia = map_col2row[jcol];
    for(int icol=0;icol<ncol-1;++icol){
      double dia0 = B[icol*nrow+jrow_dia];
      if( icol == jcol ){
        assert(fabs(dia0-1)<1.0e-30); // typically +1 or -1
        double rhs0 = B[jcol*nrow]; // the origin (0,0,...) should be valid
        if( jcol != ncol-1 ){
          if( rhs0 < -1.0e-10 ){
            std::cout << "wrong entry" << icol << " " << jrow_dia << " " << rhs0 << " " << dia0 << std::endl;
            return false;
          }
          assert(rhs0>-1.0e-10); // the initial target value can be negtive
        }
      }
      else{
        assert(fabs(dia0)<1.0e-10);
      }
    }
  }
  return true;
}

// return value
// 0: converged
// 1: input value wrong
// 2: no bound
int LinPro_SolveTable(
    int& nitr,
    std::vector<double>& B,
    std::vector<unsigned int>& map_col2row,
    unsigned int ncol,
    unsigned int nrow)
{
  assert(ncol>0);
  for(unsigned int icol=0;icol<ncol-1;++icol){
    unsigned int jrow = map_col2row[icol];
    double b0 = B[icol*nrow+jrow];
    double b1 = B[(ncol-1)*nrow+jrow];
    if( fabs(b1)<1.0e-30 ) continue;
    assert( fabs(b0) > 1.0e-30 );
    double ratio = b1/b0;
    for(unsigned int krow=0;krow<nrow;++krow){
      B[(ncol-1)*nrow+krow] -= ratio*B[icol*nrow+krow];
    }
  }
//  ::Print(B, ncol, nrow, map_col2row);
  //
  std::vector<int> flg_row(nrow,0); // 0:base 1:non_base 2:trg
  for(unsigned int ieq=0;ieq<ncol;++ieq){
    unsigned int jrow = map_col2row[ieq];
    flg_row[jrow] = 1;
  }
  flg_row[0] = 2;
  // -----
  for(int itr=0;itr<nitr;++itr){
//    std::cout << "iteration: " << itr << " " << nitr << std::endl;
//    ::Print(B,ncol,nrow,map_col2row);
    assert(LinPro_CheckTable(B, map_col2row, ncol, nrow));
    { // check convegence
      bool is_optim = true;
      for(unsigned int jrow=0;jrow<nrow;++jrow){
        if( flg_row[jrow] != 0 ){ continue; }
        double d0 = B[(ncol-1)*nrow+jrow];
        if( d0 >= 0 ){ continue; }
        is_optim = false;
      }
      if( is_optim ){ // no way to furthre increase the value
        nitr = itr;
        return 0; // converged
      }
    }
	unsigned int icol_min, jrow_min;
    { // find minimum row
      std::map<double, std::pair<unsigned int,unsigned int> > mapBoundIndex;
      for(unsigned int jrow=0;jrow<nrow;++jrow){
        if( flg_row[jrow] != 0 ){ continue; }
        double d0 = B[(ncol-1)*nrow+jrow];
        if( d0 >= 0 ){ continue; }
        double min_bound = -1;
        unsigned int icol_min_cand = UINT_MAX;
        for(unsigned int icol=0;icol<ncol-1;++icol){
          const double dia0 = B[icol*nrow+jrow];
          const double rhs0 = B[icol*nrow];
//          std::cout << "hoge" << dia0 << std::endl;
          if( dia0 <= 0 ){ continue; } // negative dia makes the rhs netgative
          double bound0 = rhs0/dia0;
          if( bound0 > 1.0e+10 ) continue;
          if( icol_min_cand == UINT_MAX || bound0 < min_bound ){
            min_bound = bound0;
            icol_min_cand = icol;
          }
        }
        if( icol_min_cand == -1 ){ continue; }
        mapBoundIndex.insert( std::make_pair(min_bound*d0,std::make_pair(icol_min_cand,jrow)) );
      }
      if( mapBoundIndex.empty() ){
        nitr = itr;
        std::cout << "no boundary" << std::endl;
        return 2; // no boundary
      }
      icol_min = mapBoundIndex.begin()->second.first;
      jrow_min = mapBoundIndex.begin()->second.second;
      assert( icol_min >= 0 && icol_min < ncol );
      assert( jrow_min >= 0 && jrow_min < nrow );
    }
    assert( fabs(B[icol_min*nrow+jrow_min]) > 0 );
    { // Gauss's sweep
      const double vpiv = B[icol_min*nrow+jrow_min];
      for(unsigned int jrow=0;jrow<nrow;++jrow){ B[icol_min*nrow+jrow] /= vpiv; }
      for(unsigned int icol=0;icol<ncol;++icol){
        if( icol == icol_min ) continue;
        const double vpiv0 = B[icol*nrow+jrow_min];
        for(unsigned int jrow=0;jrow<nrow;jrow++){
          B[icol*nrow+jrow] -= vpiv0*B[icol_min*nrow+jrow];
        }
      }
    }
    // swap basic-nonbasic index
    flg_row[jrow_min] = 1;
    unsigned int jrow_new = map_col2row[icol_min];
    flg_row[jrow_new] = 0;
    map_col2row[icol_min] = jrow_min;
  }
  return false;
}

//
// ---------------------------------------------
//

std::vector<double> delfem2::CLinPro::GetValid() const
{
  std::vector<double> res(nvar,0.0);
  const unsigned int ncol = neq+1;
  const unsigned int nrow = 1+1+nvar+nslk;
//  Print(A, ncol, nrow, map_col2row);
  assert(A.size()==ncol*nrow);
  for(unsigned int icol=0;icol<neq;++icol){
    unsigned int jrow = map_col2row[icol];
    if( jrow >= 2 && jrow < 2+nvar ){
      double rhs = A[icol*nrow];
      double dia = A[icol*nrow+jrow];
      assert(fabs(dia)>0.0);
      res[jrow-2] = rhs/dia;
    }
  }
  return res;
}

void delfem2::CLinPro::Print() const{
  const unsigned int ncol = neq+1;
  const unsigned int nrow = 1+1+nvar+nslk;
  ::Print(A, ncol, nrow, map_col2row);
}

int delfem2::CLinPro::Solve
(std::vector<double>& solution,
 double& opt_val,
 int& nitr,
 const std::vector<double>& aCoeffTrg) const
{
  const unsigned int ncol = (neq+1);
  const unsigned int nrow = (1+1+nvar+nslk);
  assert(this->A.size()==ncol*nrow);
  std::vector<unsigned int> map_col2rowB = map_col2row;
  std::vector<double> B = this->A;
  B[(ncol-1)*nrow+1] = 1;
  for(int ic=0;ic<aCoeffTrg.size();++ic){
    B[(ncol-1)*nrow+2+ic] = -aCoeffTrg[ic];
  }
//  std::cout << "before solve" << std::endl;
// ::Print(B, ncol, nrow, map_col2rowB);
  assert(LinPro_CheckTable(B, map_col2rowB, ncol, nrow));
  int res = LinPro_SolveTable(nitr, B, map_col2rowB,
                              ncol, nrow);
//  if( res != 0 ){ return res; }
//  std::cout << "res" << res << std::endl;
//  std::cout << "after solve" << std::endl;
//  ::Print(B, ncol, nrow, map_col2rowB);
  std::vector<double> buff(nrow,0.0);
  for(unsigned int icol=0;icol<ncol;++icol){
    unsigned int jrow = map_col2rowB[icol];
    buff[jrow] = B[icol*nrow];
  }
  solution.resize(nvar,0.0);
  for(unsigned int ivar=0;ivar<nvar;++ivar){
    solution[ivar] = buff[ivar+2];
//    std::cout << " sol" << ivar << " " << solution[ivar] << std::endl;
  }
  opt_val = buff[1];
//  std::cout << " opt" << opt_val << std::endl;
  return res;
}


void delfem2::CLinPro::AddEqn(
    const std::vector<double>& aW,
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
// 3 no valid solution
int delfem2::CLinPro::Precomp(int& nitr)
{
  neq = static_cast<unsigned int>(aEq.size());
  nvar = 0;
  nslk = 0;
  nart = 0;
  std::vector<int> mapEq2Slk(neq,-1);
  std::vector<int> mapEq2Art(neq,-1);
  for(unsigned int ieq=0;ieq<aEq.size();++ieq){
    const size_t nv = aEq[ieq].aCoeff.size();
    if( nvar < nv ){ nvar = static_cast<unsigned int>(nv); }
    if( aEq[ieq].itype == LE || aEq[ieq].itype == GE ){ mapEq2Slk[ieq] = nslk; nslk++; }
    ////
    if( aEq[ieq].itype == LE && aEq[ieq].rhs < 0  ){ mapEq2Art[ieq] = nart; nart++; }
    if( aEq[ieq].itype == GE && aEq[ieq].rhs > 0  ){ mapEq2Art[ieq] = nart; nart++; }
    if( aEq[ieq].itype == EQ ){ mapEq2Art[ieq] = nart; nart++; }
  }
  const unsigned int ncol = (neq+1); // neq,trg
  const unsigned int nrow = (1+1+nvar+nslk+nart); // rhs,trg,var,slk,art
//  std::cout << neq << " " << nvar << " " << nslk << " " << nart << "  " << ncol << " " << nrow << std::endl;
  A.resize(ncol*nrow,0.0);
  for(unsigned int ieq=0;ieq<neq;++ieq){
    A[ieq*nrow] = aEq[ieq].rhs;
    for(size_t ic=0;ic<aEq[ieq].aCoeff.size();++ic){
      A[ieq*nrow+2+ic] = aEq[ieq].aCoeff[ic];
    }
    if( aEq[ieq].itype == LE ){
      const int islk = mapEq2Slk[ieq]; assert(islk!=-1);
      A[ieq*nrow+2+nvar+islk] = +1;
      const int iart = mapEq2Art[ieq];
      if( iart != -1 ){
        A[ieq*nrow+2+nvar+nslk+iart] = -1;
      }
    }
    if( aEq[ieq].itype == GE ){
      const int islk = mapEq2Slk[ieq]; assert(islk!=-1);
      A[ieq*nrow+2+nvar+islk] = -1;
      const int iart = mapEq2Art[ieq];
      if( iart != -1 ){
        A[ieq*nrow+2+nvar+nslk+iart] = +1;
      }
    }
    // equality
    if( aEq[ieq].itype == EQ ){
      const int iart = mapEq2Art[ieq];
      if( iart == -1 ){ continue; }
      if( aEq[ieq].rhs >= 0 ){
        A[ieq*nrow+2+nvar+nslk+iart] = +1;
        
      }
      else{
        A[ieq*nrow+2+nvar+nslk+iart] = -1;
      }
    }
  }
  map_col2row.assign(ncol,0); // 0:base 1:non_base 2:trg
  map_col2row[ncol-1] = 1;
  for(unsigned int ieq=0;ieq<neq;ieq++){
    if( mapEq2Slk[ieq] == -1 ) continue;
    int islk = mapEq2Slk[ieq];
    map_col2row[ieq] = 1+1+nvar+islk;
  }
//  std::cout << "precomp" << std::endl;
//  ::Print(A, ncol, nrow, map_col2row);
  if( nart == 0 ){ nitr=0; return 0; }
  ////
  A[(ncol-1)*nrow+1] = 1.0;
  for(unsigned int ieq=0;ieq<aEq.size();++ieq){
    if( mapEq2Art[ieq] == -1 ) continue;
    int iart = mapEq2Art[ieq];
    unsigned int jrow1 = 1+1+nvar+nslk+iart;
    A[(ncol-1)*nrow+jrow1] = 1;
//    std::cout << "new entry: " << ieq << " " << jrow1 << std::endl;
    map_col2row[ieq] = jrow1;
  }
  for(unsigned int icol=0;icol<ncol;++icol){
    unsigned int jrow_dia = map_col2row[icol];
    if( A[icol*nrow+jrow_dia] < 0 ){
      for(unsigned int jrow=0;jrow<nrow;++jrow){
        A[icol*nrow+jrow] = -A[icol*nrow+jrow];
      }
    }
  }
//  std::cout << "  precomp art" << std::endl;
//  ::Print(A, ncol, nrow, map_col2row);
  assert(LinPro_CheckTable(A, map_col2row, ncol, nrow));
  int res = LinPro_SolveTable(nitr, A, map_col2row, ncol,nrow);
//  std::cout << "  after solve precomp" << std::endl;
//  ::Print(A, ncol, nrow, map_col2row);
  assert(LinPro_CheckTable(A, map_col2row, ncol, nrow));
  if( res == 2 ){ // no bound in solution
//    std::cout << "no bound in solution" << std::endl;
    return 2;
  }
  for(unsigned int icol=0;icol<ncol;++icol){
    unsigned int jrow_dia = map_col2row[icol];
    if( jrow_dia < 1 || jrow_dia >= 2+nvar+nslk ){
//      std::cout << "precomputation is not solved: " << jrow_dia << " " << 2+nvar+nslk << std::endl;
      return 3;
    }
  }
  assert( fabs(A[(ncol-1)*nrow]) < 1.0e-8 );
//  std::cout << "  solved" << std::endl;
//  ::(A, ncol, nrow, map_col2row);
//    std::cout << "couldn't found solution" << std::endl;
//    return 3;
//  }
//  std::cout << "succesfully found a valid solution" << std::endl;
  const std::vector<double> B = A;
  const unsigned int nrow1 = 1+1+nvar+nslk;
//  std::cout << "hugahugat" << ncol << " " << nrow1 << std::endl;
  A.resize(ncol*nrow1);
  for(unsigned int icol=0;icol<ncol;++icol){
    for(unsigned int irow=0;irow<nrow1;++irow){
      A[icol*nrow1+irow] = B[icol*nrow+irow];
    }
  }
  return 0;
}



