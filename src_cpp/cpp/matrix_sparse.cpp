#include <iostream>
#include <cassert>
#include <math.h>
#include <vector>

#include "delfem2/matrix_sparse.h"


CMatrixSquareSparse::CMatrixSquareSparse()
{
  std::cout << "MatrixSquareSparse -- construct" << std::endl;
	m_nblk_col = 0;
	m_len_col = 0;
  
  m_nblk_row = 0;
  m_len_row = 0;
  
	m_ncrs = 0;
	m_colInd = 0;
	m_rowPtr = 0;
  
  is_dia = true;
	m_valCrs = 0;
	m_valDia = 0;
}

CMatrixSquareSparse::~CMatrixSquareSparse()
{
  if( m_colInd != 0 ){ delete[] m_colInd; m_colInd = 0; }
	if( m_rowPtr != 0 ){ delete[] m_rowPtr; m_rowPtr = 0; }
  
	if( m_valCrs != 0 ){ delete[] m_valCrs; m_valCrs = 0; }
	if( m_valDia != 0 ){ delete[] m_valDia; m_valDia = 0; }
}


void CMatrixSquareSparse::Initialize(int nblk, int len, bool is_dia)
{
  this->is_dia = is_dia;
  this->m_nblk_col = nblk;
  this->m_len_col = len;
  this->m_nblk_row = nblk;
  this->m_len_row = len;
  const int blksize = len*len;
  ///
  if( m_valDia != 0 ){ delete[] m_valDia; m_valDia = 0; }
  if( is_dia ){
    m_valDia = new double [nblk*blksize];
    for(int i=0;i<nblk*blksize;i++){ m_valDia[i] = 0; }
  }
  ////
  if( m_colInd != 0 ){ delete[] m_colInd; m_colInd = 0; }
  m_colInd = new int [nblk+1];
  for(int i=0;i<nblk+1;i++){ m_colInd[i] = 0; }
  ////  
  m_ncrs = 0;
	if( m_rowPtr != 0 ){ delete[] m_rowPtr; m_rowPtr = 0; }
	if( m_valCrs != 0 ){ delete[] m_valCrs; m_valCrs = 0; }  
}

void CMatrixSquareSparse::operator = (const CMatrixSquareSparse& m)
{
  std::cout << "CMatrixSquareSparse -- copy" << std::endl;
  if( is_dia ){
    assert( m_nblk_col == m_nblk_row );
    assert( m_len_col == m_len_row );
  }
  else{
    assert( m.m_valDia == 0 );
  }
  this->is_dia = m.is_dia;  
  this->m_nblk_col = m.m_nblk_col;
  this->m_len_col  = m.m_len_col;
  this->m_nblk_row = m.m_nblk_row;
  this->m_len_row  = m.m_len_row;
  const int blksize = m_len_col*m_len_row;
  this->m_ncrs = m.m_ncrs;  
  if( m_colInd != 0 ){ delete[] m_colInd; m_colInd = 0; }
  if( m_rowPtr != 0 ){ delete[] m_rowPtr; m_rowPtr = 0; }
	if( m_valCrs != 0 ){ delete[] m_valCrs; m_valCrs = 0; }
  if( m_valDia != 0 ){ delete[] m_valDia; m_valDia = 0; }
  m_colInd = new int    [m_nblk_col+1];
  m_rowPtr = new int    [m_ncrs];
  m_valCrs = new double [m_ncrs*blksize];
  for(int i=0;i<m_nblk_col+1;      i++){ m_colInd[i] = m.m_colInd[i]; }
  for(int i=0;i<m_ncrs;        i++){ m_rowPtr[i] = m.m_rowPtr[i]; }
  for(int i=0;i<m_ncrs*blksize;i++){ m_valCrs[i] = m.m_valCrs[i]; }
  ///
  if( m.is_dia ){
    m_valDia = new double [m_nblk_col*blksize];
    for(int i=0;i<m_nblk_col*blksize;i++){ m_valDia[i] = m.m_valDia[i]; }
  }
}


bool CMatrixSquareSparse::SetZero()
{
  if( is_dia ){
    assert( m_len_col == m_len_row );
    assert( m_nblk_col == m_nblk_row );
    for(int i=0;i<m_len_col*m_len_col*m_nblk_col;i++){ m_valDia[i] = 0.0; }
  }
  else{
    assert( m_valDia == 0);
  }
  for(int i=0;i<m_len_col*m_len_row*m_ncrs;i++){ m_valCrs[i] = 0.0; }
	return true;
}

bool CMatrixSquareSparse::Mearge
(int nblkel_col, const int* blkel_col,
 int nblkel_row, const int* blkel_row,
 int blksize, const double* emat,
 std::vector<int>& marge_buffer)
{
  assert( m_colInd != 0 );
	assert( m_valCrs != 0 );
	assert( m_valDia != 0 );

//	assert( nblkel_col == nblkel_row );
  assert( blksize == m_len_col*m_len_row );
  
  if( marge_buffer.size() < m_nblk_row ){
    marge_buffer.resize(m_nblk_row);
  }

	const int* colind = m_colInd;
	const int* rowptr = m_rowPtr;
	double* vcrs = m_valCrs;
	double* vdia = m_valDia;

	for(int iblkel=0;iblkel<nblkel_col;iblkel++){
		const int iblk1 = blkel_col[iblkel];
    assert( iblk1 < m_nblk_col );
		for(int jpsup=colind[iblk1];jpsup<colind[iblk1+1];jpsup++){
			assert( jpsup < m_ncrs );
			const int jblk1 = rowptr[jpsup];
			marge_buffer[jblk1] = jpsup;
		}
		for(int jblkel=0;jblkel<nblkel_row;jblkel++){
      const int jblk1 = blkel_row[jblkel];
      assert( jblk1 < m_nblk_row );
			if( iblk1 == jblk1 ){	// Marge Diagonal
				const double* pval_in = &emat[(iblkel*nblkel_row+iblkel)*blksize];
				double* pval_out = &vdia[iblk1*blksize];
				for(int i=0;i<blksize;i++){ pval_out[i] += pval_in[i]; }
			}
			else{	// Marge Non-Diagonal
				if( marge_buffer[jblk1] == -1 ) continue;
        assert( marge_buffer[jblk1] >= 0 && marge_buffer[jblk1] < (int)m_ncrs );
				const int jpsup1 = marge_buffer[jblk1];
				assert( m_rowPtr[jpsup1] == jblk1 );
				const double* pval_in = &emat[(iblkel*nblkel_row+jblkel)*blksize];
				double* pval_out = &vcrs[jpsup1*blksize];
				for(int i=0;i<blksize;i++){ pval_out[i] += pval_in[i]; }
			}
		}
		for(int jpsup=colind[iblk1];jpsup<colind[iblk1+1];jpsup++){
			assert( jpsup < m_ncrs );
			const int jblk1 = rowptr[jpsup];
			marge_buffer[jblk1] = -1;
		}
	}
	return true;
}

void CMatrixSquareSparse::SetPattern
(const int* pColInd, int ncolind,
 const int* pRowPtr, int nrowptr)
{
  assert( m_colInd != 0 );
	assert( m_ncrs == 0 );
  assert( m_rowPtr == 0 );
  
  assert( ncolind == m_nblk_col+1 );
  for(int iblk=0;iblk<m_nblk_col+1;iblk++){
    m_colInd[iblk] = pColInd[iblk];
  }
  m_ncrs = pColInd[m_nblk_col];
  assert( m_ncrs == nrowptr );
  ////
  if( m_rowPtr != 0 ){ delete[] m_rowPtr; m_rowPtr = 0; }
  m_rowPtr = new int [m_ncrs];
  for(int icrs=0;icrs<m_ncrs;icrs++){
    m_rowPtr[icrs] = pRowPtr[icrs];
  }
  ////
  const int blksize = m_len_col*m_len_row;
  if( m_valCrs != 0 ){ delete[] m_valCrs; m_valCrs = 0; }
  m_valCrs = new double [m_ncrs*blksize];
}

// Calc Matrix Vector Product
// {y} = alpha*[A]{x} + beta*{y}
void CMatrixSquareSparse::MatVec
(double alpha,
 const std::vector<double>& x,
 double beta,
 std::vector<double>& y) const
{
	const int blksize = m_len_col*m_len_col;

	if( m_len_col == 1 && m_len_row == 1 ){
		const double* vcrs  = m_valCrs;
		const double* vdia = m_valDia;
		const int* colind = m_colInd;
		const int* rowptr = m_rowPtr;
		////////////////
		for(int iblk=0;iblk<m_nblk_col;iblk++){
			double& vy = y[iblk];
			vy *= beta;
			const int colind0 = colind[iblk];
			const int colind1 = colind[iblk+1];
			for(int icrs=colind0;icrs<colind1;icrs++){
				assert( icrs < m_ncrs );
				const int jblk0 = rowptr[icrs];
				assert( jblk0 < m_nblk_row );
				vy += alpha * vcrs[icrs] * x[jblk0];
			}
			vy += alpha * vdia[iblk] * x[iblk];
		}
	}
	else if( m_len_col == 2 && m_len_row == 2 ){
		const double* vcrs  = m_valCrs;
		const double* vdia = m_valDia;
		const int* colind = m_colInd;
		const int* rowptr = m_rowPtr;
		////////////////
		for(int iblk=0;iblk<m_nblk_col;iblk++){
			y[iblk*2+0] *= beta;
			y[iblk*2+1] *= beta;
			const int icrs0 = colind[iblk];
			const int icrs1 = colind[iblk+1];
			for(int icrs=icrs0;icrs<icrs1;icrs++){
				assert( icrs < m_ncrs );
				const int jblk0 = rowptr[icrs];
				assert( jblk0 < m_nblk_row );
				y[iblk*2+0] += alpha * ( vcrs[icrs*4  ]*x[jblk0*2+0] + vcrs[icrs*4+1]*x[jblk0*2+1] );
				y[iblk*2+1] += alpha * ( vcrs[icrs*4+2]*x[jblk0*2+0] + vcrs[icrs*4+3]*x[jblk0*2+1] );
			}
			y[iblk*2+0] += alpha * ( vdia[iblk*4+0]*x[iblk*2+0] + vdia[iblk*4+1]*x[iblk*2+1] );
			y[iblk*2+1] += alpha * ( vdia[iblk*4+2]*x[iblk*2+0] + vdia[iblk*4+3]*x[iblk*2+1] );
		}
	}
	else if( m_len_col == 3 && m_len_row == 3 ){
		const double* vcrs  = m_valCrs;
		const double* vdia = m_valDia;
		const int* colind = m_colInd;
		const int* rowptr = m_rowPtr;
		////////////////
		for(int iblk=0;iblk<m_nblk_col;iblk++){
			y[iblk*3+0] *= beta;
			y[iblk*3+1] *= beta;
			y[iblk*3+2] *= beta;
			const int icrs0 = colind[iblk];
			const int icrs1 = colind[iblk+1];
			for(int icrs=icrs0;icrs<icrs1;icrs++){
				assert( icrs < m_ncrs );
				const int jblk0 = rowptr[icrs];
				assert( jblk0 < m_nblk_row );
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
  else if( m_len_col == 4 && m_len_row == 4 ){
    const double* vcrs  = m_valCrs;
    const double* vdia = m_valDia;
    const int* colind = m_colInd;
    const int* rowptr = m_rowPtr;
    ////////////////
    for(int iblk=0;iblk<m_nblk_col;iblk++){
      y[iblk*4+0] *= beta;
      y[iblk*4+1] *= beta;
      y[iblk*4+2] *= beta;
      y[iblk*4+3] *= beta;
      const int icrs0 = colind[iblk];
      const int icrs1 = colind[iblk+1];
      for(int icrs=icrs0;icrs<icrs1;icrs++){
        assert( icrs < m_ncrs );
        const int jblk0 = rowptr[icrs];
        assert( jblk0 < m_nblk_row );
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
		const double* vcrs  = m_valCrs;
		const double* vdia = m_valDia;
		const int* colind = m_colInd;
		const int* rowptr = m_rowPtr;
		////////////////
		for(int iblk=0;iblk<m_nblk_col;iblk++){
			for(int idof=0;idof<m_len_col;idof++){ y[iblk*m_len_col+idof] *= beta; }
			const int colind0 = colind[iblk];
			const int colind1 = colind[iblk+1];
			for(int icrs=colind0;icrs<colind1;icrs++){
				assert( icrs < m_ncrs );
				const int jblk0 = rowptr[icrs];
				assert( jblk0 < m_nblk_row );
				for(int idof=0;idof<m_len_col;idof++){
				for(int jdof=0;jdof<m_len_row;jdof++){
					y[iblk*m_len_col+idof] += alpha * vcrs[icrs*blksize+idof*m_len_col+jdof] * x[jblk0*m_len_row+jdof];
				}
				}
			}
			for(int idof=0;idof<m_len_col;idof++){
			for(int jdof=0;jdof<m_len_row;jdof++){
				y[iblk*m_len_col+idof] += alpha * vdia[iblk*blksize+idof*m_len_col+jdof] * x[iblk*m_len_row+jdof];
			}
			}
		}
	}
}

void CMatrixSquareSparse::SetBoundaryCondition
(const int* bc_flag, int nbc_flag)
{
  assert( this->is_dia );
  assert( this->m_nblk_row == this->m_nblk_col );
  assert( this->m_len_row == this->m_len_col );
	const int blksize = m_len_col*m_len_row;
  assert( nbc_flag == m_nblk_col*m_len_col );
	
	for(int iblk=0;iblk<m_nblk_col;iblk++){ // set diagonal
		for(int ilen=0;ilen<m_len_col;ilen++){
      if( bc_flag[iblk*m_len_col+ilen] == 0 ) continue;
			for(int jlen=0;jlen<m_len_row;jlen++){
        m_valDia[iblk*blksize+ilen*m_len_col+jlen] = 0.0;
        m_valDia[iblk*blksize+jlen*m_len_col+ilen] = 0.0;
      }
			m_valDia[iblk*blksize+ilen*m_len_col+ilen] = 1.0;
    }
  }
  /////
  for(int iblk=0;iblk<m_nblk_col;iblk++){ // set row
    for(int icrs=m_colInd[iblk];icrs<m_colInd[iblk+1];icrs++){
      for(int ilen=0;ilen<m_len_col;ilen++){
        if( bc_flag[iblk*m_len_col+ilen] == 0 ) continue;
        for(int jlen=0;jlen<m_len_row;jlen++){
          m_valCrs[icrs*blksize+ilen*m_len_col+jlen] = 0.0;
        }
			}
    }
  }
  /////
  for(int icrs=0;icrs<m_ncrs;icrs++){ // set column
		const int jblk1 = m_rowPtr[icrs];
    for(int jlen=0;jlen<m_len_row;jlen++){
      if( bc_flag[jblk1*m_len_row+jlen] == 0 ) continue;
      for(int ilen=0;ilen<m_len_col;ilen++){
        m_valCrs[icrs*blksize+ilen*m_len_col+jlen] = 0.0;
      }
		}
	}
}

void CMatrixSquareSparse::SetMasterSlave
(const std::vector<int>& aMSFlag)
{
  assert( this->is_dia );
  assert( this->m_nblk_row == this->m_nblk_col );
  assert( this->m_len_row == this->m_len_col );
  const int blksize = m_len_col*m_len_row;
  const int ndof = m_nblk_col*m_len_col;
  if(aMSFlag.size() != m_nblk_col*m_len_col ){ return; }
  /////
  std::vector<int> row2crs(m_nblk_row,-1);
  for(int idof1=0;idof1<ndof;++idof1){ // add row
    int idof0 = aMSFlag[idof1];
    if( idof0 == -1 ) continue;
    int ino0 = idof0 / m_len_col;
    int ilen0 = idof0 - ino0*m_len_col;
    assert( ilen0 >=0 && ilen0 < m_len_col );
    assert( ino0 < m_nblk_col && ilen0 < m_len_col );
    int ino1 = idof1 / m_len_col;
    int ilen1 = idof1 - ino1*m_len_col;
    assert( ino1 < m_nblk_col && ilen1 < m_len_col );
    assert( ilen0 == ilen1 );
    for(int icrs0=m_colInd[ino0];icrs0<m_colInd[ino0+1];++icrs0){
      int jno0 = m_rowPtr[icrs0];
      assert( jno0 >= 0 && jno0 < m_nblk_row );
      row2crs[jno0] = icrs0;
    }
    for(int icrs1=m_colInd[ino1];icrs1<m_colInd[ino1+1];++icrs1){
      int jno1 = m_rowPtr[icrs1]; assert( jno1 >= 0 && jno1 < m_nblk_row );
      assert( jno1 != ino1 );
      if( jno1 != ino0 ){ // add non-diagonal 1 to non-diagonal 0
        const int icrs0 = row2crs[jno1];
        assert( icrs0 >= 0 && icrs0 < m_ncrs );
        for(int jdim=0;jdim<m_len_row;++jdim){
          m_valCrs[icrs0*blksize+ilen0*m_len_col+jdim] += m_valCrs[icrs1*blksize+ilen1*m_len_col+jdim];
        }
      }
      else{ // add non-diagonal 1 to diagonal 0
        for(int jdim=0;jdim<m_len_row;++jdim){
          m_valDia[ino0*blksize+ilen0*m_len_col+jdim] += m_valCrs[icrs1*blksize+ilen1*m_len_col+jdim];
        }
      }
    }
    { // add diagonal 1 to non-diagonal 0
      const int icrs0 = row2crs[ino1];
      assert( icrs0 >= 0 && icrs0 < m_ncrs );
      for(int jdim=0;jdim<m_len_row;++jdim){
        m_valCrs[icrs0*blksize+ilen0*m_len_col+jdim] += m_valDia[ino1*blksize+ilen1*m_len_col+jdim];
      }
    }
    for(int icrs0=m_colInd[ino0];icrs0<m_colInd[ino0+1];++icrs0){
      int jno0 = m_rowPtr[icrs0];
      assert( jno0 >= 0 && jno0 < m_nblk_row );
      row2crs[jno0] = -1;
    }
  }
  //////
  row2crs.assign(m_nblk_row,-1);
  for(int ino=0;ino<m_nblk_col;ino++){
    for(int icrs=m_colInd[ino];icrs<m_colInd[ino+1];++icrs){
      int jno0 = m_rowPtr[icrs]; assert( jno0 >= 0 && jno0 < m_nblk_row );
      row2crs[jno0] = icrs;
    }
    for(int jlen1=0;jlen1<m_len_row;jlen1++){
      int jdof0 = aMSFlag[ino*m_len_row+jlen1];
      if( jdof0 == -1 ) continue;
      int jno0 = jdof0/m_len_row;
      assert( jdof0 - jno0*m_len_row == jlen1 );
      const int icrs0 = row2crs[jno0];
      assert( icrs0 >= 0 && icrs0 < m_ncrs );
      for(int ilen=0;ilen<m_len_col;ilen++){
        m_valCrs[icrs0*blksize+ilen*m_len_col+jlen1] +=m_valDia[ino*blksize+ilen*m_len_col+jlen1];
      }
    }
    for(int icrs1=m_colInd[ino];icrs1<m_colInd[ino+1];icrs1++){
      const int jno1 = m_rowPtr[icrs1];
      assert( jno1 >= 0 && jno1 < m_nblk_row );
      for(int jlen1=0;jlen1<m_len_row;jlen1++){
        int jdof0 = aMSFlag[jno1*m_len_row+jlen1];
        if( jdof0 == -1 ) continue;
        int jno0 = jdof0/m_len_row;
        assert( jno0 >= 0 && jno0 < m_nblk_row );
        assert( jdof0 - jno0*m_len_row == jlen1 );
        if( ino == jno0 ){
          for(int ilen=0;ilen<m_len_col;ilen++){
            m_valDia[jno0*blksize+ilen*m_len_col+jlen1] += m_valCrs[icrs1*blksize+ilen*m_len_col+jlen1];
          }
        }
        else{
          const int icrs0 = row2crs[jno0];
          assert( icrs0 >= 0 && icrs0 < m_ncrs );
          for(int ilen=0;ilen<m_len_col;ilen++){
            m_valCrs[icrs0*blksize+ilen*m_len_col+jlen1] += m_valCrs[icrs1*blksize+ilen*m_len_col+jlen1];
          }
        }
      }
    }
    for(int icrs=m_colInd[ino];icrs<m_colInd[ino+1];++icrs){
      int jno0 = m_rowPtr[icrs]; assert( jno0 >= 0 && jno0 < m_nblk_row );
      row2crs[jno0] = -1;
    }
  }
  //////
  for(int iblk=0;iblk<m_nblk_col;iblk++){
    for(int ilen=0;ilen<m_len_col;ilen++){
      if( aMSFlag[iblk*m_len_col+ilen] == -1 ) continue;
      for(int jlen=0;jlen<m_len_row;jlen++){
        m_valDia[iblk*blksize+ilen*m_len_col+jlen] = 0.0;
        m_valDia[iblk*blksize+jlen*m_len_col+ilen] = 0.0;
      }
      m_valDia[iblk*blksize+ilen*m_len_col+ilen] = 1.0;
    }
  }
  
  ////
  for(int iblk=0;iblk<m_nblk_col;iblk++){
    for(int icrs=m_colInd[iblk];icrs<m_colInd[iblk+1];icrs++){
      for(int idim=0;idim<m_len_col;idim++){
        int idof0 = aMSFlag[iblk*m_len_col+idim];
        if( idof0 == -1 ) continue;
        int jblk = m_rowPtr[icrs];
        for(int jdim=0;jdim<m_len_row;jdim++){
          int idof1 = jblk*m_len_row+jdim;
          if( idof0 != idof1 ){ m_valCrs[icrs*blksize+idim*m_len_col+jdim] = +0.0; }
          else{                 m_valCrs[icrs*blksize+idim*m_len_col+jdim] = -1.0; }
          m_valCrs[icrs*blksize+idim*m_len_col+jdim] = +0.0;
        }
      }
    }
  }
  /////
  for(int iblk=0;iblk<m_nblk_col;iblk++){
    for(int icrs=m_colInd[iblk];icrs<m_colInd[iblk+1];icrs++){
      const int jblk1 = m_rowPtr[icrs];
      for(int jdim=0;jdim<m_len_row;jdim++){
        int idof0 = aMSFlag[jblk1*m_len_row+jdim];
        if( idof0 == -1 ) continue;
        for(int idim=0;idim<m_len_col;idim++){
          int idof1 = iblk*m_len_col+idim;
          if( idof0 != idof1 ){ m_valCrs[icrs*blksize+idim*m_len_col+jdim] = +0.0; }
          else{                 m_valCrs[icrs*blksize+idim*m_len_col+jdim] = -1.0; }
          m_valCrs[icrs*blksize+idim*m_len_col+jdim] = +0.0;
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

double CMatrixSquareSparse::CheckSymmetry() const
{
  assert( this->is_dia );
  assert( this->m_nblk_row == this->m_nblk_col );
  assert( this->m_len_row == this->m_len_col );
  const int blksize = m_len_col*m_len_row;
  const int nlen = m_len_col;
  ////
  double sum = 0;
  for(int ino=0;ino<m_nblk_col;++ino){
    for(int icrs0=m_colInd[ino];icrs0<m_colInd[ino+1];++icrs0){
      int jno = m_rowPtr[icrs0];
      int icrs1 = m_colInd[jno];
      for(;icrs1<m_colInd[jno+1];++icrs1){
        if( m_rowPtr[icrs1] == ino ){ break; }
      }
      if( icrs1 == m_colInd[jno+1] ){ // no counterpart
        sum += MatNorm(m_valCrs+blksize*icrs0,m_len_col,m_len_row);
      }
      else{
        sum += MatNorm_Assym(m_valCrs+blksize*icrs0,m_len_col,m_len_row,
                             m_valCrs+blksize*icrs1);
      }
    }
    sum += MatNorm_Assym(m_valDia+blksize*ino,nlen);
  }
  return sum;
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
  const int n = (int)x.size();
  assert( y.size() == n );
  for(int i=0;i<n;i++){
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
 const CMatrixSquareSparse& mat,
 std::vector<double>& r_vec,
 std::vector<double>& x_vec)
{
  assert( mat.is_dia );
  assert( mat.m_nblk_col == mat.m_nblk_row );
  assert( mat.m_len_col == mat.m_len_row );
  
	const double conv_ratio_tol = conv_ratio;
	const int mx_iter = iteration;
  
	const int nblk = mat.m_nblk_col;
  const int len = mat.m_len_col;
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
 const CMatrixSquareSparse& mat,
 std::vector<double>& r_vec,
 std::vector<double>& x_vec)
{
  assert( mat.is_dia );
  assert( mat.m_nblk_col == mat.m_nblk_row );
  assert( mat.m_len_col == mat.m_len_row );
  
  const int nblk = mat.m_nblk_col;
  const int len = mat.m_len_col;
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

