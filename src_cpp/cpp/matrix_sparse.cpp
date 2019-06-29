/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cassert>
#include <math.h>
#include <vector>

#include "delfem2/matrix_sparse.h"


CMatrixSparse::CMatrixSparse()
{
  std::cout << "MatrixSquareSparse -- construct" << std::endl;
	nblk_col = 0;
	len_col = 0;
  
  nblk_row = 0;
  len_row = 0;
  
	ncrs = 0;
  
  is_dia = true;
	valCrs = 0;
	valDia = 0;
}

CMatrixSparse::~CMatrixSparse()
{
  colInd.clear();
  rowPtr.clear();
  
	if( valCrs != 0 ){ delete[] valCrs; valCrs = 0; }
	if( valDia != 0 ){ delete[] valDia; valDia = 0; }
}


void CMatrixSparse::Initialize(int nblk, int len, bool is_dia)
{
  this->is_dia = is_dia;
  this->nblk_col = nblk;
  this->len_col = len;
  this->nblk_row = nblk;
  this->len_row = len;
  const int blksize = len*len;
  ///
  if( valDia != 0 ){ delete[] valDia; valDia = 0; }
  if( is_dia ){
    valDia = new double [nblk*blksize];
    for(int i=0;i<nblk*blksize;i++){ valDia[i] = 0; }
  }
  ////
  colInd.assign(nblk+1,0);
  ////  
  ncrs = 0;
  rowPtr.clear();
	if( valCrs != 0 ){ delete[] valCrs; valCrs = 0; }  
}

void CMatrixSparse::operator = (const CMatrixSparse& m)
{
  std::cout << "CMatrixSquareSparse -- copy" << std::endl;
  if( is_dia ){
    assert( nblk_col == nblk_row );
    assert( len_col == len_row );
  }
  else{
    assert( m.valDia == 0 );
  }
  this->is_dia = m.is_dia;  
  this->nblk_col = m.nblk_col;
  this->len_col  = m.len_col;
  this->nblk_row = m.nblk_row;
  this->len_row  = m.len_row;
  const int blksize = len_col*len_row;
  this->ncrs = m.ncrs;
  colInd.clear();
  rowPtr.clear();
	if( valCrs != 0 ){ delete[] valCrs; valCrs = 0; }
  if( valDia != 0 ){ delete[] valDia; valDia = 0; }
  colInd.resize(nblk_col+1);
  rowPtr.resize(ncrs);
  valCrs = new double [ncrs*blksize];
  for(unsigned int i=0;i<nblk_col+1;  i++){ colInd[i] = m.colInd[i]; }
  for(unsigned int i=0;i<ncrs;        i++){ rowPtr[i] = m.rowPtr[i]; }
  for(unsigned int i=0;i<ncrs*blksize;i++){ valCrs[i] = m.valCrs[i]; }
  ///
  if( m.is_dia ){
    valDia = new double [nblk_col*blksize];
    for(unsigned int i=0;i<nblk_col*blksize;i++){ valDia[i] = m.valDia[i]; }
  }
}


bool CMatrixSparse::SetZero()
{
  if( is_dia ){
    assert( len_col == len_row );
    assert( nblk_col == nblk_row );
    for(unsigned int i=0;i<len_col*len_col*nblk_col;i++){ valDia[i] = 0.0; }
  }
  else{
    assert( valDia == 0);
  }
  for(unsigned int i=0;i<len_col*len_row*ncrs;i++){ valCrs[i] = 0.0; }
	return true;
}

bool CMatrixSparse::Mearge
(unsigned int nblkel_col, const unsigned int* blkel_col,
 unsigned int nblkel_row, const unsigned int* blkel_row,
 unsigned int blksize, const double* emat,
 std::vector<int>& marge_buffer)
{
//  assert( colInd != 0 );
	assert( valCrs != 0 );
	assert( valDia != 0 );

//	assert( nblkel_col == nblkel_row );
  assert( blksize == len_col*len_row );
  
  if( marge_buffer.size() < nblk_row ){
    marge_buffer.resize(nblk_row);
  }

	const unsigned int* colind = colInd.data();
	const unsigned int* rowptr = rowPtr.data();
	double* vcrs = valCrs;
	double* vdia = valDia;

	for(unsigned int iblkel=0;iblkel<nblkel_col;iblkel++){
		const unsigned int iblk1 = blkel_col[iblkel];
    assert( iblk1 < nblk_col );
		for(unsigned int jpsup=colind[iblk1];jpsup<colind[iblk1+1];jpsup++){
			assert( jpsup < ncrs );
			const int jblk1 = rowptr[jpsup];
			marge_buffer[jblk1] = jpsup;
		}
		for(unsigned int jblkel=0;jblkel<nblkel_row;jblkel++){
      const unsigned int jblk1 = blkel_row[jblkel];
      assert( jblk1 < nblk_row );
			if( iblk1 == jblk1 ){	// Marge Diagonal
				const double* pval_in = &emat[(iblkel*nblkel_row+iblkel)*blksize];
				double* pval_out = &vdia[iblk1*blksize];
				for(unsigned int i=0;i<blksize;i++){ pval_out[i] += pval_in[i]; }
			}
			else{	// Marge Non-Diagonal
				if( marge_buffer[jblk1] == -1 ) continue;
        assert( marge_buffer[jblk1] >= 0 && marge_buffer[jblk1] < (int)ncrs );
				const int jpsup1 = marge_buffer[jblk1];
				assert( rowPtr[jpsup1] == jblk1 );
				const double* pval_in = &emat[(iblkel*nblkel_row+jblkel)*blksize];
				double* pval_out = &vcrs[jpsup1*blksize];
				for(unsigned int i=0;i<blksize;i++){ pval_out[i] += pval_in[i]; }
			}
		}
		for(unsigned int jpsup=colind[iblk1];jpsup<colind[iblk1+1];jpsup++){
			assert( jpsup < ncrs );
			const int jblk1 = rowptr[jpsup];
			marge_buffer[jblk1] = -1;
		}
	}
	return true;
}

void CMatrixSparse::SetPattern
(const int* pColInd, unsigned int ncolind,
 const int* pRowPtr, unsigned int nrowptr)
{
//  assert( colInd != 0 );
	assert( ncrs == 0 );
//  assert( rowPtr == 0 );
  
  assert( ncolind == nblk_col+1 );
  for(unsigned int iblk=0;iblk<nblk_col+1;iblk++){
    colInd[iblk] = pColInd[iblk];
  }
  ncrs = pColInd[nblk_col];
  assert( ncrs == nrowptr );
  ////
//  if( rowPtr != 0 ){ delete[] rowPtr; rowPtr = 0; }
//  rowPtr = new int [ncrs];
  rowPtr.resize(ncrs);
  for(unsigned int icrs=0;icrs<ncrs;icrs++){
    rowPtr[icrs] = pRowPtr[icrs];
  }
  ////
  const int blksize = len_col*len_row;
  if( valCrs != 0 ){ delete[] valCrs; valCrs = 0; }
  valCrs = new double [ncrs*blksize];
}

// Calc Matrix Vector Product
// {y} = alpha*[A]{x} + beta*{y}
void CMatrixSparse::MatVec
(double alpha,
 const std::vector<double>& x,
 double beta,
 std::vector<double>& y) const
{
	const int blksize = len_col*len_col;

	if( len_col == 1 && len_row == 1 ){
		const double* vcrs  = valCrs;
		const double* vdia = valDia;
		const unsigned int* colind = colInd.data();
		const unsigned int* rowptr = rowPtr.data();
		////////////////
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
			double& vy = y[iblk];
			vy *= beta;
			const unsigned int colind0 = colind[iblk];
			const unsigned int colind1 = colind[iblk+1];
			for(unsigned int icrs=colind0;icrs<colind1;icrs++){
				assert( icrs < ncrs );
				const unsigned int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk_row );
				vy += alpha * vcrs[icrs] * x[jblk0];
			}
			vy += alpha * vdia[iblk] * x[iblk];
		}
	}
	else if( len_col == 2 && len_row == 2 ){
		const double* vcrs  = valCrs;
		const double* vdia = valDia;
		const unsigned int* colind = colInd.data();
		const unsigned int* rowptr = rowPtr.data();
		////////////////
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
			y[iblk*2+0] *= beta;
			y[iblk*2+1] *= beta;
			const unsigned int icrs0 = colind[iblk];
			const unsigned int icrs1 = colind[iblk+1];
			for(unsigned int icrs=icrs0;icrs<icrs1;icrs++){
				assert( icrs < ncrs );
				const unsigned int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk_row );
				y[iblk*2+0] += alpha * ( vcrs[icrs*4  ]*x[jblk0*2+0] + vcrs[icrs*4+1]*x[jblk0*2+1] );
				y[iblk*2+1] += alpha * ( vcrs[icrs*4+2]*x[jblk0*2+0] + vcrs[icrs*4+3]*x[jblk0*2+1] );
			}
			y[iblk*2+0] += alpha * ( vdia[iblk*4+0]*x[iblk*2+0] + vdia[iblk*4+1]*x[iblk*2+1] );
			y[iblk*2+1] += alpha * ( vdia[iblk*4+2]*x[iblk*2+0] + vdia[iblk*4+3]*x[iblk*2+1] );
		}
	}
	else if( len_col == 3 && len_row == 3 ){
		const double* vcrs  = valCrs;
		const double* vdia = valDia;
		const unsigned int* colind = colInd.data();
		const unsigned int* rowptr = rowPtr.data();
		////////////////
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
			y[iblk*3+0] *= beta;
			y[iblk*3+1] *= beta;
			y[iblk*3+2] *= beta;
			const unsigned int icrs0 = colind[iblk];
			const unsigned int icrs1 = colind[iblk+1];
			for(unsigned int icrs=icrs0;icrs<icrs1;icrs++){
				assert( icrs < ncrs );
				const unsigned int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk_row );
        const int i0 = iblk*3;
        const int j0 = jblk0*3;
        const int k0 = icrs*9;
				y[i0+0] += alpha*(vcrs[k0+0]*x[j0+0]+vcrs[k0+1]*x[j0+1]+vcrs[k0+2]*x[j0+2]);
				y[i0+1] += alpha*(vcrs[k0+3]*x[j0+0]+vcrs[k0+4]*x[j0+1]+vcrs[k0+5]*x[j0+2]);
				y[i0+2] += alpha*(vcrs[k0+6]*x[j0+0]+vcrs[k0+7]*x[j0+1]+vcrs[k0+8]*x[j0+2]);
			}
      {
        const int i0 = iblk*3;
        const int k0 = iblk*9;
        y[i0+0] += alpha*(vdia[k0+0]*x[i0+0]+vdia[k0+1]*x[i0+1]+vdia[k0+2]*x[i0+2]);
        y[i0+1] += alpha*(vdia[k0+3]*x[i0+0]+vdia[k0+4]*x[i0+1]+vdia[k0+5]*x[i0+2]);
        y[i0+2] += alpha*(vdia[k0+6]*x[i0+0]+vdia[k0+7]*x[i0+1]+vdia[k0+8]*x[i0+2]);
      }
		}
  }
  else if( len_col == 4 && len_row == 4 ){
    const double* vcrs  = valCrs;
    const double* vdia = valDia;
    const unsigned int* colind = colInd.data();
    const unsigned int* rowptr = rowPtr.data();
    ////////////////
    for(unsigned int iblk=0;iblk<nblk_col;iblk++){
      y[iblk*4+0] *= beta;
      y[iblk*4+1] *= beta;
      y[iblk*4+2] *= beta;
      y[iblk*4+3] *= beta;
      const unsigned int icrs0 = colind[iblk];
      const unsigned int icrs1 = colind[iblk+1];
      for(unsigned int icrs=icrs0;icrs<icrs1;icrs++){
        assert( icrs < ncrs );
        const unsigned int jblk0 = rowptr[icrs];
        assert( jblk0 < nblk_row );
        const int i0 = iblk*4;
        const int j0 = jblk0*4;
        const int k0 = icrs*16;
        y[i0+0] += alpha*(vcrs[k0+ 0]*x[j0+0]+vcrs[k0+ 1]*x[j0+1]+vcrs[k0+ 2]*x[j0+2]+vcrs[k0+ 3]*x[j0+3]);
        y[i0+1] += alpha*(vcrs[k0+ 4]*x[j0+0]+vcrs[k0+ 5]*x[j0+1]+vcrs[k0+ 6]*x[j0+2]+vcrs[k0+ 7]*x[j0+3]);
        y[i0+2] += alpha*(vcrs[k0+ 8]*x[j0+0]+vcrs[k0+ 9]*x[j0+1]+vcrs[k0+10]*x[j0+2]+vcrs[k0+11]*x[j0+3]);
        y[i0+3] += alpha*(vcrs[k0+12]*x[j0+0]+vcrs[k0+13]*x[j0+1]+vcrs[k0+14]*x[j0+2]+vcrs[k0+15]*x[j0+3]);
      }
      {
        const int i0 = iblk*4;
        const int k0 = iblk*16;
        y[i0+0] += alpha*(vdia[k0+ 0]*x[i0+0]+vdia[k0+ 1]*x[i0+1]+vdia[k0+ 2]*x[i0+2]+vdia[k0+ 3]*x[i0+3]);
        y[i0+1] += alpha*(vdia[k0+ 4]*x[i0+0]+vdia[k0+ 5]*x[i0+1]+vdia[k0+ 6]*x[i0+2]+vdia[k0+ 7]*x[i0+3]);
        y[i0+2] += alpha*(vdia[k0+ 8]*x[i0+0]+vdia[k0+ 9]*x[i0+1]+vdia[k0+10]*x[i0+2]+vdia[k0+11]*x[i0+3]);
        y[i0+3] += alpha*(vdia[k0+12]*x[i0+0]+vdia[k0+13]*x[i0+1]+vdia[k0+14]*x[i0+2]+vdia[k0+15]*x[i0+3]);
      }
    }
  }
	else{
		const double* vcrs  = valCrs;
		const double* vdia = valDia;
		const unsigned int* colind = colInd.data();
		const unsigned int* rowptr = rowPtr.data();
		////////////////
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
			for(unsigned int idof=0;idof<len_col;idof++){ y[iblk*len_col+idof] *= beta; }
			const unsigned int colind0 = colind[iblk];
			const unsigned int colind1 = colind[iblk+1];
			for(unsigned int icrs=colind0;icrs<colind1;icrs++){
				assert( icrs < ncrs );
				const unsigned int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk_row );
				for(unsigned int idof=0;idof<len_col;idof++){
				for(unsigned int jdof=0;jdof<len_row;jdof++){
					y[iblk*len_col+idof] += alpha * vcrs[icrs*blksize+idof*len_col+jdof] * x[jblk0*len_row+jdof];
				}
				}
			}
			for(unsigned int idof=0;idof<len_col;idof++){
			for(unsigned int jdof=0;jdof<len_row;jdof++){
				y[iblk*len_col+idof] += alpha * vdia[iblk*blksize+idof*len_col+jdof] * x[iblk*len_row+jdof];
			}
			}
		}
	}
}

void CMatrixSparse::SetBoundaryCondition
(const int* bc_flag, int np, int ndimval)
{
  assert( this->is_dia );
  assert( this->nblk_row == this->nblk_col );
  assert( this->len_row == this->len_col );
	const int blksize = len_col*len_row;
  assert( np == nblk_col );
  assert( ndimval == len_col );
	
	for(unsigned int iblk=0;iblk<nblk_col;iblk++){ // set diagonal
		for(unsigned int ilen=0;ilen<len_col;ilen++){
      if( bc_flag[iblk*len_col+ilen] == 0 ) continue;
			for(unsigned int jlen=0;jlen<len_row;jlen++){
        valDia[iblk*blksize+ilen*len_col+jlen] = 0.0;
        valDia[iblk*blksize+jlen*len_col+ilen] = 0.0;
      }
			valDia[iblk*blksize+ilen*len_col+ilen] = 1.0;
    }
  }
  /////
  for(unsigned int iblk=0;iblk<nblk_col;iblk++){ // set row
    for(unsigned int icrs=colInd[iblk];icrs<colInd[iblk+1];icrs++){
      for(unsigned int ilen=0;ilen<len_col;ilen++){
        if( bc_flag[iblk*len_col+ilen] == 0 ) continue;
        for(unsigned int jlen=0;jlen<len_row;jlen++){
          valCrs[icrs*blksize+ilen*len_col+jlen] = 0.0;
        }
			}
    }
  }
  /////
  for(unsigned int icrs=0;icrs<ncrs;icrs++){ // set column
		const int jblk1 = rowPtr[icrs];
    for(unsigned int jlen=0;jlen<len_row;jlen++){
      if( bc_flag[jblk1*len_row+jlen] == 0 ) continue;
      for(unsigned int ilen=0;ilen<len_col;ilen++){
        valCrs[icrs*blksize+ilen*len_col+jlen] = 0.0;
      }
		}
	}
}

void CMatrixSparse::SetMasterSlave
(const int* aMSFlag)
{
  assert( this->is_dia );
  assert( this->nblk_row == this->nblk_col );
  assert( this->len_row == this->len_col );
  const unsigned int blksize = len_col*len_row;
  const unsigned int ndof = nblk_col*len_col;
  /////
  std::vector<int> row2crs(nblk_row,-1);
  for(unsigned int idof1=0;idof1<ndof;++idof1){ // add row
    int idof0 = aMSFlag[idof1];
    if( idof0 == -1 ) continue;
    int ino0 = idof0 / len_col;
    int ilen0 = idof0 - ino0*len_col;
    assert( ilen0 >=0 && ilen0 < (int)len_col );
    assert( ino0 < (int)nblk_col && ilen0 < (int)len_col );
    int ino1 = idof1 / len_col;
    int ilen1 = idof1 - ino1*len_col;
    assert( ino1 < (int)nblk_col && ilen1 < (int)len_col );
    assert( ilen0 == ilen1 );
    for(unsigned int icrs0=colInd[ino0];icrs0<colInd[ino0+1];++icrs0){
      int jno0 = rowPtr[icrs0];
      assert( jno0 >= 0 && jno0 < (int)nblk_row );
      row2crs[jno0] = icrs0;
    }
    for(unsigned int icrs1=colInd[ino1];icrs1<colInd[ino1+1];++icrs1){
      int jno1 = rowPtr[icrs1];
      assert( jno1 >= 0 && jno1 < (int)nblk_row );
      assert( jno1 != ino1 );
      if( jno1 != ino0 ){ // add non-diagonal 1 to non-diagonal 0
        const int icrs0 = row2crs[jno1];
        assert( icrs0 >= 0 && icrs0 < (int)ncrs );
        for(unsigned int jdim=0;jdim<len_row;++jdim){
          valCrs[icrs0*blksize+ilen0*len_col+jdim] += valCrs[icrs1*blksize+ilen1*len_col+jdim];
        }
      }
      else{ // add non-diagonal 1 to diagonal 0
        for(unsigned int jdim=0;jdim<len_row;++jdim){
          valDia[ino0*blksize+ilen0*len_col+jdim] += valCrs[icrs1*blksize+ilen1*len_col+jdim];
        }
      }
    }
    { // add diagonal 1 to non-diagonal 0
      const int icrs0 = row2crs[ino1];
      assert( icrs0 >= 0 && icrs0 < (int)ncrs );
      for(unsigned int jdim=0;jdim<len_row;++jdim){
        valCrs[icrs0*blksize+ilen0*len_col+jdim] += valDia[ino1*blksize+ilen1*len_col+jdim];
      }
    }
    for(unsigned int icrs0=colInd[ino0];icrs0<colInd[ino0+1];++icrs0){
      int jno0 = rowPtr[icrs0];
      assert( jno0 >= 0 && jno0 < (int)nblk_row );
      row2crs[jno0] = -1;
    }
  }
  //////
  row2crs.assign(nblk_row,-1);
  for(unsigned int ino=0;ino<nblk_col;ino++){
    for(unsigned int icrs=colInd[ino];icrs<colInd[ino+1];++icrs){
      int jno0 = rowPtr[icrs];
      assert( jno0 >= 0 && jno0 < (int)nblk_row );
      row2crs[jno0] = icrs;
    }
    for(unsigned int jlen1=0;jlen1<len_row;jlen1++){
      int jdof0 = aMSFlag[ino*len_row+jlen1];
      if( jdof0 == -1 ) continue;
      int jno0 = jdof0/len_row;
      assert( jdof0 - jno0*len_row == jlen1 );
      const int icrs0 = row2crs[jno0];
      assert( icrs0 >= 0 && icrs0 < (int)ncrs );
      for(unsigned int ilen=0;ilen<len_col;ilen++){
        valCrs[icrs0*blksize+ilen*len_col+jlen1] +=valDia[ino*blksize+ilen*len_col+jlen1];
      }
    }
    for(unsigned int icrs1=colInd[ino];icrs1<colInd[ino+1];icrs1++){
      const unsigned int jno1 = rowPtr[icrs1];
      assert( jno1 >= 0 && jno1 < nblk_row );
      for(unsigned int jlen1=0;jlen1<len_row;jlen1++){
        int jdof0 = aMSFlag[jno1*len_row+jlen1];
        if( jdof0 == -1 ) continue;
        int jno0 = jdof0/len_row;
        assert( jno0 >= 0 && jno0 < (int)nblk_row );
        assert( jdof0 - jno0*len_row == jlen1 );
        if( (int)ino == jno0 ){
          for(unsigned int ilen=0;ilen<len_col;ilen++){
            valDia[jno0*blksize+ilen*len_col+jlen1] += valCrs[icrs1*blksize+ilen*len_col+jlen1];
          }
        }
        else{
          const int icrs0 = row2crs[jno0];
          assert( icrs0 >= 0 && icrs0 < (int)ncrs );
          for(unsigned int ilen=0;ilen<len_col;ilen++){
            valCrs[icrs0*blksize+ilen*len_col+jlen1] += valCrs[icrs1*blksize+ilen*len_col+jlen1];
          }
        }
      }
    }
    for(unsigned int icrs=colInd[ino];icrs<colInd[ino+1];++icrs){
      unsigned int jno0 = rowPtr[icrs];
      assert( jno0 >= 0 && jno0 < nblk_row );
      row2crs[jno0] = -1;
    }
  }
  //////
  for(unsigned int iblk=0;iblk<nblk_col;iblk++){
    for(unsigned int ilen=0;ilen<len_col;ilen++){
      if( aMSFlag[iblk*len_col+ilen] == -1 ) continue;
      for(unsigned int jlen=0;jlen<len_row;jlen++){
        valDia[iblk*blksize+ilen*len_col+jlen] = 0.0;
        valDia[iblk*blksize+jlen*len_col+ilen] = 0.0;
      }
      valDia[iblk*blksize+ilen*len_col+ilen] = 1.0;
    }
  }
  
  ////
  for(unsigned int iblk=0;iblk<nblk_col;iblk++){
    for(unsigned int icrs=colInd[iblk];icrs<colInd[iblk+1];icrs++){
      for(unsigned int idim=0;idim<len_col;idim++){
        int idof0 = aMSFlag[iblk*len_col+idim];
        if( idof0 == -1 ) continue;
        int jblk = rowPtr[icrs];
        for(unsigned int jdim=0;jdim<len_row;jdim++){
          int idof1 = jblk*len_row+jdim;
          if( idof0 != idof1 ){ valCrs[icrs*blksize+idim*len_col+jdim] = +0.0; }
          else{                 valCrs[icrs*blksize+idim*len_col+jdim] = -1.0; }
          valCrs[icrs*blksize+idim*len_col+jdim] = +0.0;
        }
      }
    }
  }
  /////
  for(unsigned int iblk=0;iblk<nblk_col;iblk++){
    for(unsigned int icrs=colInd[iblk];icrs<colInd[iblk+1];icrs++){
      const int jblk1 = rowPtr[icrs];
      for(unsigned int jdim=0;jdim<len_row;jdim++){
        int idof0 = aMSFlag[jblk1*len_row+jdim];
        if( idof0 == -1 ) continue;
        for(unsigned int idim=0;idim<len_col;idim++){
          int idof1 = iblk*len_col+idim;
          if( idof0 != idof1 ){ valCrs[icrs*blksize+idim*len_col+jdim] = +0.0; }
          else{                 valCrs[icrs*blksize+idim*len_col+jdim] = -1.0; }
          valCrs[icrs*blksize+idim*len_col+jdim] = +0.0;
        }
      }
    }
  }
}

double MatNorm_Assym
(const double* V0, int n0, int m0,
 const double* V1)
{
  double s = 0.0;
  for(int i=0;i<n0;++i){
    for(int j=0;j<m0;++j){
      double v0 = V0[i*m0+j];
      double v1 = V1[j*n0+i];
      s += (v0-v1)*(v0-v1);
    }
  }
  return s;
}

double MatNorm
(const double* V, int n, int m)
{
  double s = 0.0;
  for(int i=0;i<n;++i){
    for(int j=0;j<m;++j){
      double v = V[i*m+j];
      s += v*v;
    }
  }
  return s;
}

double MatNorm_Assym
(const double* V, int n)
{
  double s = 0.0;
  for(int i=0;i<n;++i){
    for(int j=0;j<n;++j){
      double v0 = V[i*n+j];
      double v1 = V[j*n+i];
      s += (v0-v1)*(v0-v1);
    }
  }
  return s;
}

double CMatrixSparse::CheckSymmetry() const
{
  assert( this->nblk_row == this->nblk_col );
  assert( this->len_row == this->len_col );
  const int blksize = len_col*len_row;
  const int nlen = len_col;
  ////
  double sum = 0;
  for(unsigned int ino=0;ino<nblk_col;++ino){
    for(unsigned int icrs0=colInd[ino];icrs0<colInd[ino+1];++icrs0){
      int jno = rowPtr[icrs0];
      unsigned int icrs1 = colInd[jno];
      for(;icrs1<colInd[jno+1];++icrs1){
        if( rowPtr[icrs1] == ino ){ break; }
      }
      if( icrs1 == colInd[jno+1] ){ // no counterpart
        sum += MatNorm(valCrs+blksize*icrs0,len_col,len_row);
      }
      else{
        sum += MatNorm_Assym(valCrs+blksize*icrs0,len_col,len_row,
                             valCrs+blksize*icrs1);
      }
    }
    sum += MatNorm_Assym(valDia+blksize*ino,nlen);
  }
  return sum;
}

void CMatrixSparse::ScaleLeftRight(const double* scale){
  assert( this->nblk_row == this->nblk_col );
  assert( this->len_row == this->len_col );
  const int blksize = len_col*len_row;
  for(unsigned int ino=0;ino<nblk_col;++ino){
    for(unsigned int icrs0=colInd[ino];icrs0<colInd[ino+1];++icrs0){
      const int jno = rowPtr[icrs0];
      const double s0 = scale[ino]*scale[jno];
      for(int i=0;i<blksize;++i){ valCrs[icrs0*blksize+i] *= s0; }
    }
  }
  if( is_dia ){
    for(unsigned int ino=0;ino<nblk_col;++ino){
      double s0 = scale[ino]*scale[ino];
      for(int i=0;i<blksize;++i){ valDia[ino*blksize+i] *= s0; }
    }
  }
}

void CMatrixSparse::AddDia(double eps){
  assert( this->nblk_row == this->nblk_col );
  assert( this->len_row == this->len_col );
  const int blksize = len_col*len_row;
  const int nlen = this->len_col;
  if( !is_dia ){ return; }
  for(unsigned int ino=0;ino<nblk_col;++ino){
    for(int ilen=0;ilen<nlen;++ilen){
      valDia[ino*blksize+ilen*nlen+ilen] += eps;
    }
  }
}

//////////////////////////////////////////////////////////////////////////

double InnerProduct
(const std::vector<double>& r_vec,
 const std::vector<double>& u_vec)
{
  const unsigned int n = r_vec.size();
  assert( u_vec.size() == n );
  double r = 0.0;
  for(unsigned int i=0;i<n;i++){
    r += r_vec[i]*u_vec[i];
  }
  return r;
}

double InnerProduct
(const double* r_vec,
 const double* u_vec,
 int n)
{
  double r = 0.0;
  for(int i=0;i<n;i++){
    r += r_vec[i]*u_vec[i];
  }
  return r;
}

// {y} = {y} + a * {x}
void AXPY
(double a,
 const std::vector<double>& x,
 std::vector<double>& y)
{
  const unsigned int n = x.size();
  assert( y.size() == n );
  for(unsigned int i=0;i<n;i++){
    y[i] += a*x[i];
  }
}

// {y} = {y} + a * {x}
void AXPY
(double a,
 const double* x,
 double* y,
 int n)
{
  for(int i=0;i<n;i++){
    y[i] += a*x[i];
  }
}

void Solve_CG
(double& conv_ratio,
 int& iteration,
 const CMatrixSparse& mat,
 std::vector<double>& r_vec,
 std::vector<double>& x_vec)
{
  assert( mat.is_dia );
  assert( mat.nblk_col == mat.nblk_row );
  assert( mat.len_col == mat.len_row );
  
	const double conv_ratio_tol = conv_ratio;
	const int mx_iter = iteration;
  
	const unsigned int nblk = mat.nblk_col;
  const unsigned int len = mat.len_col;
  assert(r_vec.size() == nblk*len);
  const int ndof = nblk*len;
  
  // {x} = 0
  x_vec.assign(ndof,0.0);
  
  double sqnorm_res = InnerProduct(r_vec,r_vec);
  if( sqnorm_res < 1.0e-30 ){
    conv_ratio = 0.0;
    iteration = 0;
    return;
  }
	double inv_sqnorm_res_ini = 1.0 / sqnorm_res;
  
  std::vector<double> Ap_vec(ndof);
  
	// Set Initial Serch Direction
	// {p} = {r}
  std::vector<double>  p_vec = r_vec;
  
	iteration = mx_iter;
	for(int iitr=1;iitr<mx_iter;iitr++){
    
		double alpha;
		{	// alpha = (r,r) / (p,Ap)
			mat.MatVec(1.0,p_vec,0.0,Ap_vec);
			const double pAp = InnerProduct(p_vec,Ap_vec);
			alpha = sqnorm_res / pAp;
		}
    
		// update x
		// {x} = +alpha*{ p} + {x}
		AXPY(alpha,p_vec,x_vec);
    
		// {r} = -alpha*{Ap} + {r}
		AXPY(-alpha,Ap_vec,r_vec);
    
		double sqnorm_res_new = InnerProduct(r_vec,r_vec);
    // Converge Judgement    
    
		if( sqnorm_res_new * inv_sqnorm_res_ini < conv_ratio_tol*conv_ratio_tol ){
      conv_ratio = sqrt( sqnorm_res * inv_sqnorm_res_ini );
      iteration = iitr;
      return;
		}
    
		// beta = (r1,r1) / (r0,r0)
		const double beta = sqnorm_res_new / sqnorm_res;
		sqnorm_res = sqnorm_res_new;
    
    // {p} = {r} + beta*{p}
    for(int i=0;i<ndof;i++){ p_vec[i] = r_vec[i] + beta*p_vec[i]; }
	}
  
	return;
}


bool Solve_BiCGSTAB
(double& conv_ratio,
 int& num_iter,
 const CMatrixSparse& mat,
 std::vector<double>& r_vec,
 std::vector<double>& x_vec)
{
  assert( mat.is_dia );
  assert( mat.nblk_col == mat.nblk_row );
  assert( mat.len_col == mat.len_row );
  
  const unsigned int nblk = mat.nblk_col;
  const unsigned int len = mat.len_col;
  assert(r_vec.size() == nblk*len);
  const int ndof = nblk*len;
  
  const unsigned int max_iter = num_iter;
  const double tolerance = conv_ratio;

  std::vector<double> s_vec(ndof);
  std::vector<double> As_vec(ndof);
  std::vector<double> p_vec(ndof);
  std::vector<double> Ap_vec(ndof);
  std::vector<double> r2_vec(ndof);
  
  x_vec.assign(ndof,0.0);
  
  double sq_inv_norm_res_ini;
  {
    const double sq_norm_res_ini = InnerProduct(r_vec, r_vec);
    std::cout << "initial norm " << sq_norm_res_ini << std::endl;
    if( sq_norm_res_ini < 1.0e-60 ){
      conv_ratio = sqrt( sq_norm_res_ini );
      num_iter = 0;
      return true;
    }
    sq_inv_norm_res_ini = 1.0 / sq_norm_res_ini;
  }
  
  // {r2} = {r}
  r2_vec = r_vec;
  
  // {p} = {r}
  p_vec = r_vec;
  
  // calc ({r},{r2})
  double r_r2 = InnerProduct(r_vec,r2_vec);
  
  num_iter = max_iter;
  for(unsigned int iitr=1;iitr<max_iter;iitr++){
    
    // calc {Ap} = [A]*{p}
    mat.MatVec(1.0,p_vec,0.0,Ap_vec);
//    ls.MATVEC(1.0,ip,0.0,iAp);
//    std::cout << " sq_norm iAp : " << ls.DOT(iAp,iAp) << std::endl;
//    std::cout << " sq_norm iAp : " << InnerProduct(Ap_vec,Ap_vec) << std::endl;
    
    // calc alpha
    // alhpa = ({r},{r2}) / ({Ap},{r2})
    double alpha;
    {
      const double denominator = InnerProduct(Ap_vec,r2_vec);
//      std::cout << " alpha deno : " << denominator << std::endl;
      alpha = r_r2 / denominator;
    }
    
    // {s} = {r} - alpha*{Ap}
    s_vec = r_vec;
    AXPY(-alpha, Ap_vec, s_vec);
    
    // calc {As} = [A]*{s}
    mat.MatVec(1.0,s_vec,0.0,As_vec);
//    ls.MATVEC(1.0,is,0.0,iAs);
    
    // calc omega
    // omega = ({As},{s}) / ({As},{As})
    double omega;
    {
      const double denominator = InnerProduct(As_vec,As_vec);
      const double numerator = InnerProduct(As_vec,s_vec);
      omega = numerator / denominator;
    }
    
    // update solution
    // ix += alpha*{p} + omega*{s}
    AXPY(alpha,p_vec,x_vec);
    AXPY(omega,s_vec,x_vec);
    
    // update residual
    // {r} = {s} - omega*{As}
    r_vec = s_vec;
    AXPY(-omega,As_vec,r_vec);
    
    {
      const double sq_norm_res = InnerProduct(r_vec,r_vec);
      const double sq_conv_ratio = sq_norm_res * sq_inv_norm_res_ini;
//      std::cout << iitr << " " << sq_norm_res << " " << sqrt(sq_conv_ratio) << " " << sqrt(sq_norm_res) << std::endl;
      if( sq_conv_ratio < tolerance*tolerance ){
        conv_ratio = sqrt( sq_conv_ratio );
        num_iter = iitr;
        break;
      }
    }
    
    // calc beta
    // beta = ({r},{r2})^new/({r},{r2})^old * alpha / omega
    double beta;
    {
      const double tmp1 = InnerProduct(r_vec,r2_vec);
      beta = (tmp1*alpha) / (r_r2*omega);
      r_r2 = tmp1;
    }
    
    // update p_vector
    // {p} = {r} + beta*({p}-omega*[A]*{p})
    for(int i=0;i<ndof;++i){ p_vec[i] *= beta; }
    AXPY(1.0,r_vec,p_vec);
    AXPY(-beta*omega,Ap_vec,p_vec);
//    ls.SCAL(beta,ip);
//    ls.AXPY(1.0,ir,ip);
//    ls.AXPY(-beta*omega,iAp,ip);
  }
  
  return true;
}

