#include <stdio.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/matrix_sparse.h"
#include "delfem2/ilu_sparse.h"
#include "delfem2/fem.h"
#include "delfem2/mshtopo.h"

namespace py = pybind11;

void MatrixSquareSparse_SetPattern
(CMatrixSquareSparse& mss,
 const py::array_t<int>& psup_ind,
 const py::array_t<int>& psup)
{
  assert( psup_ind.ndim()  == 1 );
  assert( psup.ndim()  == 1 );
  const int np = mss.m_nblk_col;
  assert( psup_ind.shape()[0] == np+1 );
  mss.SetPattern(psup_ind.data(), psup_ind.shape()[0],
                 psup.data(),     psup.shape()[0]);
}

void MatrixSquareSparse_SetFixBC
(CMatrixSquareSparse& mss,
 const py::array_t<int>& flagbc)
{
  mss.SetBoundaryCondition(flagbc.data(),flagbc.shape()[0],flagbc.shape()[1]);
}


void PrecondILU0
(CPreconditionerILU&  mat_ilu,
 const CMatrixSquareSparse& mss)
{
  mat_ilu.Initialize_ILU0(mss);
}


std::vector<double> PySolve_PCG
(py::array_t<double>& vec_b,
 py::array_t<double>& vec_x,
 double conv_ratio, double iteration,
 const CMatrixSquareSparse& mat_A,
 const CPreconditionerILU& ilu_A)
{
  //  std::cout << "solve pcg" << std::endl;
  auto buff_vecb = vec_b.request();
  auto buff_vecx = vec_x.request();
  return Solve_PCG((double*)buff_vecb.ptr,
                   (double*)buff_vecx.ptr,
                   conv_ratio,iteration,
                   mat_A,ilu_A);
}

std::vector<double> PySolve_PBiCGStab
(py::array_t<double>& vec_b,
 py::array_t<double>& vec_x,
 double conv_ratio, double iteration,
 const CMatrixSquareSparse& mat_A,
 const CPreconditionerILU& ilu_A)
{
  //  std::cout << "solve pcg" << std::endl;
  auto buff_vecb = vec_b.request();
  auto buff_vecx = vec_x.request();
  return Solve_PBiCGStab((double*)buff_vecb.ptr,
                         (double*)buff_vecx.ptr,
                         conv_ratio,iteration,
                         mat_A,ilu_A);
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void PyMergeLinSys_Poission
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double alpha, double source,
 const py::array_t<double>& aXY,
 const py::array_t<int>& aElm,
 MESHELEM_TYPE elem_type,
 const py::array_t<double>& aVal)
{
  assert( aXY.shape()[1] == 2 || aXY.shape()[1] == 3 );
  assert( nNodeElem(elem_type) == aElm.shape()[1] );
  auto buff_vecb = vec_b.request();
  if( aXY.shape()[1] == 2 ){
    if( elem_type == MESHELEM_TRI ){
      MergeLinSys_Poission_MeshTri2D(mss, (double*)buff_vecb.ptr,
                                     alpha, source,
                                     aXY.data(), aXY.shape()[0],
                                     aElm.data(), aElm.shape()[0],
                                     aVal.data());
    }
  }
  if( aXY.shape()[1] == 3 ){
    if( elem_type == MESHELEM_TET ){
      MergeLinSys_Poission_MeshTet3D(mss, (double*)buff_vecb.ptr,
                                     alpha, source,
                                     aXY.data(), aXY.shape()[0],
                                     aElm.data(), aElm.shape()[0],
                                     aVal.data());
    }
  }
}


void PyMergeLinSys_Diffuse
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double alpha, double rho, double source,
 double dt_timestep, double gamma_newmark,
 const py::array_t<double>& aXY,
 const py::array_t<int>& aElm,
 MESHELEM_TYPE elem_type,
 const py::array_t<double>& aVal,
 const py::array_t<double>& aVelo)
{
  assert( aXY.shape()[1] == 2 || aXY.shape()[1] == 3 );
  assert( nNodeElem(elem_type) == aElm.shape()[1] );
  auto buff_vecb = vec_b.request();
  if( aXY.shape()[1] == 2 ){
    if( elem_type == MESHELEM_TRI ){
      MergeLinSys_Diffusion_MeshTri2D(mss, (double*)buff_vecb.ptr,
                                      alpha, rho, source,
                                      dt_timestep, gamma_newmark,
                                      aXY.data(), aXY.shape()[0],
                                      aElm.data(), aElm.shape()[0],
                                      aVal.data(), aVelo.data());
    }
  }
  else if( aXY.shape()[1] == 3 ){
    if( elem_type == MESHELEM_TET ){
      MergeLinSys_Diffusion_MeshTet3D(mss, (double*)buff_vecb.ptr,
                                      alpha, rho, source,
                                      dt_timestep, gamma_newmark,
                                      aXY.data(), aXY.shape()[0],
                                      aElm.data(), aElm.shape()[0],
                                      aVal.data(), aVelo.data());
    }
  }
}




void PyMergeLinSys_LinearSolidStatic
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double myu, double lambda, double rho,
 std::vector<double>& gravity,
 const py::array_t<double>& aXY,
 const py::array_t<int>& aElm,
 MESHELEM_TYPE elem_type,
 const py::array_t<double>& aVal)
{
  assert( aXY.shape()[1] == 2 || aXY.shape()[1] == 3 );
  assert( nNodeElem(elem_type) == aElm.shape()[1] );
  auto buff_vecb = vec_b.request();
  if( aXY.shape()[1] == 2 ){
    if( elem_type == MESHELEM_TRI ){
      MergeLinSys_SolidStaticLinear_MeshTri2D(mss, (double*)buff_vecb.ptr,
                                              myu,lambda,rho,
                                              gravity[0], gravity[1],
                                              aXY.data(), aXY.shape()[0],
                                              aElm.data(), aElm.shape()[0],
                                              aVal.data());
    }
  }
  if( aXY.shape()[1] == 3 ){
    if( elem_type == MESHELEM_TET ){
      MergeLinSys_SolidStaticLinear_MeshTet3D(mss, (double*)buff_vecb.ptr,
                                              myu,lambda,rho,
                                              gravity[0],gravity[1],gravity[2],
                                              aXY.data(), aXY.shape()[0],
                                              aElm.data(), aElm.shape()[0],
                                              aVal.data());
    }
  }
}


void PyMergeLinSys_LinearSolidDynamic
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double myu, double lambda, double rho,
 std::vector<double>& gravity,
 double dt_timestep, double gamma_newmark, double beta_newmark,
 const py::array_t<double>& aXY,
 const py::array_t<int>& aElm,
 MESHELEM_TYPE elem_type,
 const py::array_t<double>& aVal,
 const py::array_t<double>& aVelo,
 const py::array_t<double>& aAcc)
{
  auto buff_vecb = vec_b.request();
  assert( aXY.shape()[1] == 2 || aXY.shape()[1] == 3 );
  assert( nNodeElem(elem_type) == aElm.shape()[1] );
  if( aXY.shape()[1] == 2 ){
    if( elem_type == MESHELEM_TRI ){
      MergeLinSys_SolidDynamicLinear_MeshTri2D(mss,(double*)buff_vecb.ptr,
                                               myu,lambda,rho,gravity[0],gravity[1],
                                               dt_timestep,gamma_newmark,beta_newmark,
                                               aXY.data(), aXY.shape()[0],
                                               aElm.data(), aElm.shape()[0],
                                               aVal.data(),aVelo.data(),aAcc.data());
    }
  }
  if( aXY.shape()[1] == 3 ){
    if( elem_type == MESHELEM_TET ){
      MergeLinSys_SolidDynamicLinear_MeshTet3D(mss,(double*)buff_vecb.ptr,
                                               myu,lambda,rho,gravity[0],gravity[1],gravity[2],
                                               dt_timestep,gamma_newmark,beta_newmark,
                                               aXY.data(), aXY.shape()[0],
                                               aElm.data(), aElm.shape()[0],
                                               aVal.data(),aVelo.data(),aAcc.data());
    }
  }
}


void PyMergeLinSys_StorksStatic2D
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double myu, double g_x, double g_y,
 const py::array_t<double>& aXY,
 const py::array_t<int>& aTri,
 const py::array_t<double>& aVal)
{
  auto buff_vecb = vec_b.request();
  MergeLinSys_StokesStatic2D(mss,(double*)buff_vecb.ptr,
                             myu,g_x,g_y,
                             aXY.data(), aXY.shape()[0],
                             aTri.data(), aTri.shape()[0],
                             aVal.data());
}

void PyMergeLinSys_StorksDynamic2D
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double myu, double rho, double g_x, double g_y,
 double dt_timestep, double gamma_newmark,
 const py::array_t<double>& aXY,
 const py::array_t<int>& aTri,
 const py::array_t<double>& aVal,
 const py::array_t<double>& aVelo)
{
  auto buff_vecb = vec_b.request();
  MergeLinSys_StokesDynamic2D(mss,(double*)buff_vecb.ptr,
                              myu,rho,g_x,g_y,
                              dt_timestep, gamma_newmark,
                              aXY.data(), aXY.shape()[0],
                              aTri.data(), aTri.shape()[0],
                              aVal.data(),aVelo.data());
}

void PyMergeLinSys_NavierStorks2D
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double myu, double rho, double g_x, double g_y,
 double dt_timestep, double gamma_newmark,
 const py::array_t<double>& aXY,
 const py::array_t<int>& aTri,
 const py::array_t<double>& aVal,
 const py::array_t<double>& aVelo)
{
  auto buff_vecb = vec_b.request();
  MergeLinSys_NavierStokes2D(mss,(double*)buff_vecb.ptr,
                             myu,rho,g_x,g_y,
                             dt_timestep, gamma_newmark,
                             aXY.data(), aXY.shape()[0],
                             aTri.data(), aTri.shape()[0],
                             aVal.data(), aVelo.data());
}

double PyMergeLinSys_Cloth
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double lambda, double myu, double stiff_bend,
 const py::array_t<double>& aPosIni,
 const py::array_t<int>& aTri,
 const py::array_t<int>& aQuad,
 const py::array_t<double>& aXYZ)
{
  auto buff_vecb = vec_b.request();
  double W = MergeLinSys_Cloth(mss,(double*)buff_vecb.ptr,
                               lambda, myu, stiff_bend,
                               aPosIni.data(), aPosIni.shape()[0], aPosIni.shape()[1],
                               aTri.data(), aTri.shape()[0],
                               aQuad.data(), aQuad.shape()[0],
                               aXYZ.data());
  return W;
}

double PyMergeLinSys_MassPoint
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double mass_point,
 double dt,
 const std::vector<double>& gravity,
 const py::array_t<double>& aXYZ,
 const py::array_t<double>& aUVW)
{
  double* pB = (double*)(vec_b.request().ptr);
  double W = 0.0;
  const int np = aXYZ.shape()[0];
  assert(aUVW.shape()[0] == np);
  for(int ip=0;ip<np;ip++){
    const double c[3] = {aXYZ.at(ip,0),aXYZ.at(ip,1),aXYZ.at(ip,2)};
    W -= mass_point*( c[0]*gravity[0] + c[1]*gravity[1] + c[2]*gravity[2] );
    pB[ip*3+0] -= mass_point*gravity[0];
    pB[ip*3+1] -= mass_point*gravity[1];
    pB[ip*3+2] -= mass_point*gravity[2];
  }
  const int ndof = aXYZ.size();
  const double* pUVW = aUVW.data();
  for(int i=0;i<ndof;i++){
    pB[i] = -pB[i] + mass_point*pUVW[i]/dt;
  }
  for(int ip=0;ip<np;ip++){
    mss.m_valDia[ip*9+0*3+0] += mass_point / (dt*dt);
    mss.m_valDia[ip*9+1*3+1] += mass_point / (dt*dt);
    mss.m_valDia[ip*9+2*3+2] += mass_point / (dt*dt);
  }
  return W;
}


void init_fem(py::module &m){
  py::class_<CMatrixSquareSparse>(m,"MatrixSquareSparse")
  .def(py::init<>())
  .def("initialize", &CMatrixSquareSparse::Initialize)
  .def("setZero",    &CMatrixSquareSparse::SetZero);
  
  py::class_<CPreconditionerILU>(m,"PreconditionerILU")
  .def(py::init<>())
  .def("ilu_decomp", &CPreconditionerILU::DoILUDecomp)
  .def("set_value", &CPreconditionerILU::SetValueILU);
  //  .def(py::init<const CPreconditionerILU&>);
  
  m.def("matrixSquareSparse_setPattern", &MatrixSquareSparse_SetPattern);
  m.def("matrixSquareSparse_setFixBC", &MatrixSquareSparse_SetFixBC);
  m.def("precond_ilu0",  &PrecondILU0);
  m.def("linsys_solve_pcg", &PySolve_PCG);
  m.def("linsys_solve_bicgstab",&PySolve_PBiCGStab);
  
  m.def("mergeLinSys_poission", &PyMergeLinSys_Poission);
  m.def("mergeLinSys_diffuse",&PyMergeLinSys_Diffuse);
  m.def("mergeLinSys_linearSolidStatic",&PyMergeLinSys_LinearSolidStatic);
  m.def("mergeLinSys_linearSolidDynamic",&PyMergeLinSys_LinearSolidDynamic);
  m.def("mergeLinSys_storksStatic2D",&PyMergeLinSys_StorksStatic2D);
  m.def("mergeLinSys_storksDynamic2D",&PyMergeLinSys_StorksDynamic2D);
  m.def("mergeLinSys_navierStorks2D",&PyMergeLinSys_NavierStorks2D);
  m.def("mergeLinSys_cloth", &PyMergeLinSys_Cloth);
  m.def("mergeLinSys_massPoint",&PyMergeLinSys_MassPoint);

}
