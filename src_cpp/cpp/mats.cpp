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

#include "delfem2/mats.h"

// Calc Matrix Vector Product
// {y} = alpha*[A]{x} + beta*{y}
template <>
void CMatrixSparse<double>::MatVec
(double alpha,
 const std::vector<double>& x,
 double beta,
 std::vector<double>& y) const
{
	const int blksize = len_col*len_col;

	if( len_col == 1 && len_row == 1 ){
		const double* vcrs  = valCrs.data();
		const double* vdia = valDia.data();
		const unsigned int* colind = colInd.data();
		const unsigned int* rowptr = rowPtr.data();
		////////////////
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
			double& vy = y[iblk];
			vy *= beta;
			const unsigned int colind0 = colind[iblk];
			const unsigned int colind1 = colind[iblk+1];
			for(unsigned int icrs=colind0;icrs<colind1;icrs++){
				assert( icrs < rowPtr.size() );
				const unsigned int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk_row );
				vy += alpha * vcrs[icrs] * x[jblk0];
			}
			vy += alpha * vdia[iblk] * x[iblk];
		}
	}
	else if( len_col == 2 && len_row == 2 ){
		const double* vcrs  = valCrs.data();
		const double* vdia = valDia.data();
		const unsigned int* colind = colInd.data();
		const unsigned int* rowptr = rowPtr.data();
		////////////////
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
			y[iblk*2+0] *= beta;
			y[iblk*2+1] *= beta;
			const unsigned int icrs0 = colind[iblk];
			const unsigned int icrs1 = colind[iblk+1];
			for(unsigned int icrs=icrs0;icrs<icrs1;icrs++){
				assert( icrs < rowPtr.size() );
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
		const double* vcrs  = valCrs.data();
		const double* vdia = valDia.data();
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
				assert( icrs < rowPtr.size() );
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
    const double* vcrs  = valCrs.data();
    const double* vdia = valDia.data();
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
        assert( icrs < rowPtr.size() );
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
		const double* vcrs  = valCrs.data();
		const double* vdia = valDia.data();
		const unsigned int* colind = colInd.data();
		const unsigned int* rowptr = rowPtr.data();
		////////////////
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
			for(unsigned int idof=0;idof<len_col;idof++){ y[iblk*len_col+idof] *= beta; }
			const unsigned int colind0 = colind[iblk];
			const unsigned int colind1 = colind[iblk+1];
			for(unsigned int icrs=colind0;icrs<colind1;icrs++){
				assert( icrs < rowPtr.size() );
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


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void SetMasterSlave
(CMatrixSparse<double>& mat,
 const int* aMSFlag)
{
  assert( !mat.valDia.empty() );
  assert( mat.nblk_row == mat.nblk_col );
  assert( mat.len_row == mat.len_col );
  const unsigned int len = mat.len_col;
  const unsigned int nblk = mat.nblk_col;
  const unsigned int blksize = len*len;
  const unsigned int ndof = nblk*len;
  /////
  std::vector<int> row2crs(nblk,-1);
  for(unsigned int idof1=0;idof1<ndof;++idof1){ // add row
    int idof0 = aMSFlag[idof1];
    if( idof0 == -1 ) continue;
    int ino0 = idof0 / len;
    int ilen0 = idof0 - ino0*len;
    assert( ilen0 >=0 && ilen0 < (int)len );
    assert( ino0 < (int)nblk && ilen0 < (int)len );
    int ino1 = idof1 / len;
    int ilen1 = idof1 - ino1*len;
    assert( ino1 < (int)nblk && ilen1 < (int)len );
    assert( ilen0 == ilen1 );
    for(unsigned int icrs0=mat.colInd[ino0];icrs0<mat.colInd[ino0+1];++icrs0){
      int jno0 = mat.rowPtr[icrs0];
      assert( jno0 >= 0 && jno0 < (int)nblk );
      row2crs[jno0] = icrs0;
    }
    for(unsigned int icrs1=mat.colInd[ino1];icrs1<mat.colInd[ino1+1];++icrs1){
      int jno1 = mat.rowPtr[icrs1];
      assert( jno1 >= 0 && jno1 < (int)nblk );
      assert( jno1 != ino1 );
      if( jno1 != ino0 ){ // add non-diagonal 1 to non-diagonal 0
        const int icrs0 = row2crs[jno1];
        assert( icrs0 >= 0 && icrs0 < (int)mat.rowPtr.size() );
        for(unsigned int jdim=0;jdim<len;++jdim){
          mat.valCrs[icrs0*blksize+ilen0*len+jdim] += mat.valCrs[icrs1*blksize+ilen1*len+jdim];
        }
      }
      else{ // add non-diagonal 1 to diagonal 0
        for(unsigned int jdim=0;jdim<len;++jdim){
          mat.valDia[ino0*blksize+ilen0*len+jdim] += mat.valCrs[icrs1*blksize+ilen1*len+jdim];
        }
      }
    }
    { // add diagonal 1 to non-diagonal 0
      const int icrs0 = row2crs[ino1];
      assert( icrs0 >= 0 && icrs0 < (int)mat.rowPtr.size() );
      for(unsigned int jdim=0;jdim<len;++jdim){
        mat.valCrs[icrs0*blksize+ilen0*len+jdim] += mat.valDia[ino1*blksize+ilen1*len+jdim];
      }
    }
    for(unsigned int icrs0=mat.colInd[ino0];icrs0<mat.colInd[ino0+1];++icrs0){
      int jno0 = mat.rowPtr[icrs0];
      assert( jno0 >= 0 && jno0 < (int)nblk );
      row2crs[jno0] = -1;
    }
  }
  //////
  row2crs.assign(nblk,-1);
  for(unsigned int ino=0;ino<nblk;ino++){
    for(unsigned int icrs=mat.colInd[ino];icrs<mat.colInd[ino+1];++icrs){
      int jno0 = mat.rowPtr[icrs];
      assert( jno0 >= 0 && jno0 < (int)nblk );
      row2crs[jno0] = icrs;
    }
    for(unsigned int jlen1=0;jlen1<len;jlen1++){
      int jdof0 = aMSFlag[ino*len+jlen1];
      if( jdof0 == -1 ) continue;
      int jno0 = jdof0/len;
      assert( jdof0 - jno0*len == jlen1 );
      const int icrs0 = row2crs[jno0];
      assert( icrs0 >= 0 && icrs0 < (int)mat.rowPtr.size() );
      for(unsigned int ilen=0;ilen<len;ilen++){
        mat.valCrs[icrs0*blksize+ilen*len+jlen1] += mat.valDia[ino*blksize+ilen*len+jlen1];
      }
    }
    for(unsigned int icrs1=mat.colInd[ino];icrs1<mat.colInd[ino+1];icrs1++){
      const unsigned int jno1 = mat.rowPtr[icrs1];
      assert( jno1 < nblk );
      for(unsigned int jlen1=0;jlen1<len;jlen1++){
        int jdof0 = aMSFlag[jno1*len+jlen1];
        if( jdof0 == -1 ) continue;
        int jno0 = jdof0/len;
        assert( jno0 >= 0 && jno0 < (int)nblk );
        assert( jdof0 - jno0*len == jlen1 );
        if( (int)ino == jno0 ){
          for(unsigned int ilen=0;ilen<len;ilen++){
            mat.valDia[jno0*blksize+ilen*len+jlen1] += mat.valCrs[icrs1*blksize+ilen*len+jlen1];
          }
        }
        else{
          const int icrs0 = row2crs[jno0];
          assert( icrs0 >= 0 && icrs0 < (int)mat.rowPtr.size() );
          for(unsigned int ilen=0;ilen<len;ilen++){
            mat.valCrs[icrs0*blksize+ilen*len+jlen1] += mat.valCrs[icrs1*blksize+ilen*len+jlen1];
          }
        }
      }
    }
    for(unsigned int icrs=mat.colInd[ino];icrs<mat.colInd[ino+1];++icrs){
      unsigned int jno0 = mat.rowPtr[icrs];
      assert( jno0 < nblk );
      row2crs[jno0] = -1;
    }
  }
  //////
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int ilen=0;ilen<len;ilen++){
      if( aMSFlag[iblk*len+ilen] == -1 ) continue;
      for(unsigned int jlen=0;jlen<len;jlen++){
        mat.valDia[iblk*blksize+ilen*len+jlen] = 0.0;
        mat.valDia[iblk*blksize+jlen*len+ilen] = 0.0;
      }
      mat.valDia[iblk*blksize+ilen*len+ilen] = 1.0;
    }
  }
  
  ////
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int icrs=mat.colInd[iblk];icrs<mat.colInd[iblk+1];icrs++){
      for(unsigned int idim=0;idim<len;idim++){
        int idof0 = aMSFlag[iblk*len+idim];
        if( idof0 == -1 ) continue;
        int jblk = mat.rowPtr[icrs];
        for(unsigned int jdim=0;jdim<len;jdim++){
          int idof1 = jblk*len+jdim;
          if( idof0 != idof1 ){ mat.valCrs[icrs*blksize+idim*len+jdim] = +0.0; }
          else{                 mat.valCrs[icrs*blksize+idim*len+jdim] = -1.0; }
          mat.valCrs[icrs*blksize+idim*len+jdim] = +0.0;
        }
      }
    }
  }
  /////
  for(unsigned int iblk=0;iblk<nblk;iblk++){
    for(unsigned int icrs=mat.colInd[iblk];icrs<mat.colInd[iblk+1];icrs++){
      const int jblk1 = mat.rowPtr[icrs];
      for(unsigned int jdim=0;jdim<len;jdim++){
        int idof0 = aMSFlag[jblk1*len+jdim];
        if( idof0 == -1 ) continue;
        for(unsigned int idim=0;idim<len;idim++){
          int idof1 = iblk*len+idim;
          if( idof0 != idof1 ){ mat.valCrs[icrs*blksize+idim*len+jdim] = +0.0; }
          else{                 mat.valCrs[icrs*blksize+idim*len+jdim] = -1.0; }
          mat.valCrs[icrs*blksize+idim*len+jdim] = +0.0;
        }
      }
    }
  }
}

void ScaleLeftRight
(CMatrixSparse<double>& mat,
 const double* scale){
  assert( mat.nblk_row == mat.nblk_col );
  assert( mat.len_row == mat.len_col );
  const unsigned int nblk = mat.nblk_col;
  const unsigned int len = mat.len_col;
  const unsigned int blksize = len*len;
  for(unsigned int ino=0;ino<nblk;++ino){
    for(unsigned int icrs0=mat.colInd[ino];icrs0<mat.colInd[ino+1];++icrs0){
      const int jno = mat.rowPtr[icrs0];
      const double s0 = scale[ino]*scale[jno];
      for(unsigned int i=0;i<blksize;++i){ mat.valCrs[icrs0*blksize+i] *= s0; }
    }
  }
  if( !mat.valDia.empty() ){
    for(unsigned int ino=0;ino<nblk;++ino){
      double s0 = scale[ino]*scale[ino];
      for(unsigned int i=0;i<blksize;++i){ mat.valDia[ino*blksize+i] *= s0; }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

double Dot
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

double DotX
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
 const CMatrixSparse<double>& mat,
 std::vector<double>& r_vec,
 std::vector<double>& x_vec)
{
  assert( !mat.valDia.empty() );
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
  
  double sqnorm_res = Dot(r_vec,r_vec);
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
			const double pAp = Dot(p_vec,Ap_vec);
			alpha = sqnorm_res / pAp;
		}
    
		// update x
		// {x} = +alpha*{ p} + {x}
		AXPY(alpha,p_vec,x_vec);
    
		// {r} = -alpha*{Ap} + {r}
		AXPY(-alpha,Ap_vec,r_vec);
    
		double sqnorm_res_new = Dot(r_vec,r_vec);
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
 const CMatrixSparse<double>& mat,
 std::vector<double>& r_vec,
 std::vector<double>& x_vec)
{
  assert( !mat.valDia.empty() );
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
    const double sq_norm_res_ini = Dot(r_vec, r_vec);
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
  double r_r2 = Dot(r_vec,r2_vec);
  
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
      const double denominator = Dot(Ap_vec,r2_vec);
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
      const double denominator = Dot(As_vec,As_vec);
      const double numerator = Dot(As_vec,s_vec);
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
      const double sq_norm_res = Dot(r_vec,r_vec);
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
      const double tmp1 = Dot(r_vec,r2_vec);
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

void XPlusAY
(std::vector<double>& X,
 const int nDoF,
 const std::vector<int>& aBCFlag,
 double alpha,
 const std::vector<double>& Y)
{
  for(int i=0;i<nDoF;++i ){
    if( aBCFlag[i] !=0 ) continue;
    X[i] += alpha*Y[i];
  }
}

void XPlusAYBZ
(std::vector<double>& X,
 const int nDoF,
 const std::vector<int>& aBCFlag,
 double alpha,
 const std::vector<double>& Y,
 double beta,
 const std::vector<double>& Z)
{
  for(int i=0;i<nDoF;++i ){
    if( aBCFlag[i] !=0 ) continue;
    X[i] += alpha*Y[i] + beta*Z[i];
  }
}

void XPlusAYBZCW
(std::vector<double>& X,
 const int nDoF,
 const std::vector<int>& aBCFlag,
 double alpha,
 const std::vector<double>& Y,
 double beta,
 const std::vector<double>& Z,
 double gamma,
 const std::vector<double>& W)
{
  for(int i=0;i<nDoF;++i ){
    if( aBCFlag[i] !=0 ) continue;
    X[i] += alpha*Y[i] + beta*Z[i] + gamma*W[i];
  }
}

// set boundary condition
void setRHS_Zero
(std::vector<double>& vec_b,
 const std::vector<int>& aBCFlag,
 int iflag_nonzero)
{
  const int ndof = (int)vec_b.size();
  for (int i=0;i<ndof;++i){
    if (aBCFlag[i]==iflag_nonzero) continue;
    vec_b[i] = 0;
  }
}

void setRHS_MasterSlave
(double* vec_b,
 int nDoF,
 const int* aMSFlag)
{
  for(int idof=0;idof<nDoF;++idof){
    int jdof = aMSFlag[idof];
    if( jdof == -1 ) continue;
    vec_b[jdof] += vec_b[idof];
    vec_b[idof] = 0;
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

double CheckSymmetry(const CMatrixSparse<double>& mat)
{
  assert( mat.nblk_row == mat.nblk_col );
  assert( mat.len_row == mat.len_col );
  const int blksize = mat.len_col*mat.len_row;
  const int nlen = mat.len_col;
  ////
  double sum = 0;
  for(unsigned int ino=0;ino<mat.nblk_col;++ino){
    for(unsigned int icrs0=mat.colInd[ino];icrs0<mat.colInd[ino+1];++icrs0){
      int jno = mat.rowPtr[icrs0];
      unsigned int icrs1 = mat.colInd[jno];
      for(;icrs1<mat.colInd[jno+1];++icrs1){
        if( mat.rowPtr[icrs1] == ino ){ break; }
      }
      if( icrs1 == mat.colInd[jno+1] ){ // no counterpart
        sum += MatNorm(mat.valCrs.data()+blksize*icrs0, mat.len_col, mat.len_row);
      }
      else{
        sum += MatNorm_Assym(mat.valCrs.data()+blksize*icrs0, mat.len_col, mat.len_row,
                             mat.valCrs.data()+blksize*icrs1);
      }
    }
    sum += MatNorm_Assym(mat.valDia.data()+blksize*ino,nlen);
  }
  return sum;
}
