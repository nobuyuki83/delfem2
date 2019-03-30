#include <stdio.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/matrix_sparse.h"
#include "delfem2/ilu_sparse.h"
#include "delfem2/fem.h"

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

void PyMergeLinSys_Poission2D
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double alpha, double source,
 const py::array_t<double>& aXY,
 const py::array_t<int>& aTri,
 const py::array_t<double>& aVal)
{
  auto buff_vecb = vec_b.request();
  MergeLinSys_Poission2D(mss, (double*)buff_vecb.ptr,
                         alpha, source,
                         aXY.data(), aXY.shape()[0],
                         aTri.data(), aTri.shape()[0],
                         aVal.data());
}




void PyMergeLinSys_Diffuse2D
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double alpha, double rho, double source,
 double dt_timestep, double gamma_newmark,
 const py::array_t<double>& aXY,
 const py::array_t<int>& aTri,
 const py::array_t<double>& aVal,
 const py::array_t<double>& aVelo)
{
  auto buff_vecb = vec_b.request();
  MergeLinSys_Diffusion2D(mss, (double*)buff_vecb.ptr,
                          alpha, rho, source,
                          dt_timestep, gamma_newmark,
                          aXY.data(), aXY.shape()[0],
                          aTri.data(), aTri.shape()[0],
                          aVal.data(), aVelo.data());
}




void PyMergeLinSys_LinearSolidStatic2D
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double myu, double lambda, double rho, double g_x, double g_y,
 const py::array_t<double>& aXY,
 const py::array_t<int>& aTri,
 const py::array_t<double>& aVal)
{
  auto buff_vecb = vec_b.request();
  MergeLinSys_LinearSolid2D_Static(mss, (double*)buff_vecb.ptr,
                                   myu,lambda,rho,g_x,g_y,
                                   aXY.data(), aXY.shape()[0],
                                   aTri.data(), aTri.shape()[0],
                                   aVal.data());
  
}


void PyMergeLinSys_LinearSolidDynamic2D
(CMatrixSquareSparse& mss,
 py::array_t<double>& vec_b,
 double myu, double lambda, double rho, double g_x, double g_y,
 double dt_timestep, double gamma_newmark, double beta_newmark,
 const py::array_t<double>& aXY,
 const py::array_t<int>& aTri,
 const py::array_t<double>& aVal,
 const py::array_t<double>& aVelo,
 const py::array_t<double>& aAcc)
{
  auto buff_vecb = vec_b.request();
  MergeLinSys_LinearSolid2D_Dynamic(mss,(double*)buff_vecb.ptr,
                                    myu,lambda,rho,g_x,g_y,
                                    dt_timestep,gamma_newmark,beta_newmark,
                                    aXY.data(), aXY.shape()[0],
                                    aTri.data(), aTri.shape()[0],
                                    aVal.data(),aVelo.data(),aAcc.data());
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
  
  m.def("mergeLinSys_poission2D", &PyMergeLinSys_Poission2D);
  m.def("mergeLinSys_diffuse2D",&PyMergeLinSys_Diffuse2D);
  m.def("mergeLinSys_linearSolidStatic2D",&PyMergeLinSys_LinearSolidStatic2D);
  m.def("mergeLinSys_linearSolidDynamic2D",&PyMergeLinSys_LinearSolidDynamic2D);
  m.def("mergeLinSys_storksStatic2D",&PyMergeLinSys_StorksStatic2D);
  m.def("mergeLinSys_storksDynamic2D",&PyMergeLinSys_StorksDynamic2D);
  m.def("mergeLinSys_navierStorks2D",&PyMergeLinSys_NavierStorks2D);

}
