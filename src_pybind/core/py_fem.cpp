#include <stdio.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "delfem2/matrix_sparse.h"
#include "delfem2/ilu_sparse.h"
#include "delfem2/fem.h"
#include "delfem2/fem_ematrix.h"
#include "delfem2/mshtopo.h"
#include "delfem2/sdf.h"

#include "delfem2/objfunc_v23.h"
#include "delfem2/dyntri_v2.h"

namespace py = pybind11;

void MatrixSquareSparse_SetPattern
(CMatrixSparse& mss,
 const py::array_t<int>& psup_ind,
 const py::array_t<int>& psup)
{
  assert( mss.nblk_col == mss.nblk_row );
  assert( mss.len_col == mss.len_row );
  assert( psup_ind.ndim()  == 1 );
  assert( psup.ndim()  == 1 );
  const int np = mss.nblk_col;
  assert( psup_ind.shape()[0] == np+1 );
  mss.SetPattern(psup_ind.data(), psup_ind.shape()[0],
                 psup.data(),     psup.shape()[0]);
}

void MatrixSquareSparse_SetFixBC
(CMatrixSparse& mss,
 const py::array_t<int>& flagbc)
{
  assert( mss.nblk_col == mss.nblk_row );
  assert( mss.len_col == mss.len_row );
  assert( flagbc.ndim() == 2 );
  assert( flagbc.shape()[0] == mss.nblk_col );
  assert( flagbc.shape()[1] == mss.len_col );
  mss.SetBoundaryCondition(flagbc.data(),flagbc.shape()[0],flagbc.shape()[1]);
}


void MatrixSquareSparse_ScaleLeftRight
(CMatrixSparse& mss,
 const py::array_t<double>& scale)
{
  assert( mss.nblk_col == mss.nblk_row );
  assert( mss.len_col == mss.len_row );
  assert( scale.ndim() == 1 );
  assert( scale.shape()[0] == mss.nblk_col );
  mss.ScaleLeftRight(scale.data());
}

void LinearSystem_SetMasterSlave
(CMatrixSparse& mss,
 py::array_t<double>& np_b,
 const py::array_t<int>& np_ms)
{
  assert( mss.nblk_col == mss.nblk_row );
  assert( mss.len_col == mss.len_row );
  assert( np_b.ndim() == 2 );
  assert( np_b.shape()[0] == mss.nblk_col );
  assert( np_b.shape()[1] == mss.len_col );
  assert( np_ms.ndim() == 2 );
  assert( np_ms.shape()[0] == np_b.shape()[0] );
  assert( np_ms.shape()[1] == np_b.shape()[1] );
  mss.SetMasterSlave(np_ms.data());
  auto buff_b = np_b.request();
  setRHS_MasterSlave((double*)buff_b.ptr,
                     np_b.shape()[0]*np_b.shape()[1], np_ms.data());
}


void PrecondILU0
(CPreconditionerILU&  mat_ilu,
 const CMatrixSparse& mss)
{
  mat_ilu.Initialize_ILU0(mss);
}


std::vector<double> PySolve_PCG
(py::array_t<double>& vec_b,
 py::array_t<double>& vec_x,
 double conv_ratio, double iteration,
 const CMatrixSparse& mat_A,
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
 const CMatrixSparse& mat_A,
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
(CMatrixSparse& mss,
 py::array_t<double>& vec_b,
 double alpha, double source,
 const py::array_t<double>& aXY,
 const py::array_t<unsigned int>& aElm,
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
(CMatrixSparse& mss,
 py::array_t<double>& vec_b,
 double alpha, double rho, double source,
 double dt_timestep, double gamma_newmark,
 const py::array_t<double>& aXY,
 const py::array_t<unsigned int>& aElm,
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
(CMatrixSparse& mss,
 py::array_t<double>& vec_b,
 double myu, double lambda, double rho,
 std::vector<double>& gravity,
 const py::array_t<double>& aXY,
 const py::array_t<unsigned int>& aElm,
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
(CMatrixSparse& mss,
 py::array_t<double>& vec_b,
 double myu, double lambda, double rho,
 std::vector<double>& gravity,
 double dt_timestep, double gamma_newmark, double beta_newmark,
 const py::array_t<double>& aXY,
 const py::array_t<unsigned int>& aElm,
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
(CMatrixSparse& mss,
 py::array_t<double>& vec_b,
 double myu, double g_x, double g_y,
 const py::array_t<double>& aXY,
 const py::array_t<unsigned int>& aTri,
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
(CMatrixSparse& mss,
 py::array_t<double>& vec_b,
 double myu, double rho, double g_x, double g_y,
 double dt_timestep, double gamma_newmark,
 const py::array_t<double>& aXY,
 const py::array_t<unsigned int>& aTri,
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
(CMatrixSparse& mss,
 py::array_t<double>& vec_b,
 double myu, double rho, double g_x, double g_y,
 double dt_timestep, double gamma_newmark,
 const py::array_t<double>& aXY,
 const py::array_t<unsigned int>& aTri,
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
(CMatrixSparse& mss,
 py::array_t<double>& vec_b,
 double lambda, double myu, double stiff_bend,
 const py::array_t<double>& aPosIni,
 const py::array_t<unsigned int>& aTri,
 const py::array_t<unsigned int>& aQuad,
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


double PyMergeLinSys_Contact
(CMatrixSparse& mss,
 py::array_t<double>& vec_b,
 ////
 double stiff_contact,
 double contact_clearance,
 const std::vector<const CSDF3*>& apSDF,
 const py::array_t<double>& aXYZ)
{
  if( apSDF.size() == 0 ) return 0;
  class CMyInput : public CInput_Contact
  {
  public:
    CMyInput(const std::vector<const CSDF3*>& apSDF){ this->apSDF = apSDF; }
    virtual double penetrationNormal(double& nx, double& ny, double& nz,
                                     double px, double py, double pz) const
    {
      double n[3];
      double max_pd = apSDF[0]->Projection(px,py,pz, n);
      /*
      for(unsigned int ipct=1;ipct<apSDF.size();ipct++){
        double dist0,n0[3];
        dist0 = apSDF[ipct]->Projection(px,py,pz, n0);
        if( dist0 < max_pd ) continue;
        max_pd = dist0;
        n[0] = n0[0];
        n[1] = n0[1];
        n[2] = n0[2];
      }
       */
      nx = -n[0];
      ny = -n[1];
      nz = -n[2];
      return max_pd;
    }
  public:
    std::vector<const CSDF3*> apSDF;
  } input(apSDF);
  auto buff_vecb = vec_b.request();
  double W = MergeLinSys_Contact(mss, (double*)buff_vecb.ptr,
                               stiff_contact,contact_clearance,
                               input,
                               aXYZ.data(), aXYZ.shape()[0]);
  return W;
}

double PyMergeLinSys_MassPoint
(CMatrixSparse& mss,
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
    mss.valDia[ip*9+0*3+0] += mass_point / (dt*dt);
    mss.valDia[ip*9+1*3+1] += mass_point / (dt*dt);
    mss.valDia[ip*9+2*3+2] += mass_point / (dt*dt);
  }
  return W;
}

std::tuple<py::array_t<int>,py::array_t<int>> PyAddMasterSlavePattern
(const py::array_t<int>& ms_flag,
 const py::array_t<int>& np_psup_ind0,
 const py::array_t<int>& np_psup0)
{
  assert(ms_flag.shape()[0] == np_psup_ind0.shape()[0]-1);
  assert(ms_flag.ndim() == 2 );
  std::cout <<  ms_flag.shape()[0] << " " <<  ms_flag.shape()[1] << std::endl;
  std::vector<int> psup_ind, psup;
  JArray_AddMasterSlavePattern(psup_ind, psup,
                               ms_flag.data(), ms_flag.shape()[1],
                               np_psup_ind0.data(), np_psup_ind0.shape()[0], np_psup0.data());
  py::array_t<int> np_psup_ind((int)psup_ind.size(),psup_ind.data());
  py::array_t<int> np_psup((int)psup.size(),psup.data());
  return std::tie(np_psup_ind,np_psup);
}

void PyMasterSlave_DistributeValue
(py::array_t<double>& val,
 const py::array_t<int>& ms_flag)
{
  double* pVal = (double*)(val.request().ptr);
  const int nDoF = ms_flag.size();
  for(int idof=0;idof<nDoF;++idof){
    int jdof = ms_flag.data()[idof];
    if( jdof == -1 ) continue;
    assert( jdof >= 0 && jdof < nDoF );
    pVal[ idof] = pVal[ jdof];
  }
}


void PyPBD_ConstProj_Rigid2D
(py::array_t<double>& npXYt,
 double stiffness,
 const py::array_t<int>& npClstrInd,
 const py::array_t<int>& npClstr,
 const py::array_t<double>& npXY)
{
  PBD_ConstProj_Rigid2D((double*)(npXYt.request().ptr),
                  stiffness,
                  npClstrInd.data(), npClstrInd.size(),
                  npClstr.data(),    npClstr.size(),
                  npXY.data(),       npXY.shape()[0]);
}


void PyConstProj_Rigid3D
(py::array_t<double>& npXYZt,
 double stiffness,
 const py::array_t<int>& npClstrInd,
 const py::array_t<int>& npClstr,
 const py::array_t<double>& npXYZ)
{
  PBD_ConstProj_Rigid3D((double*)(npXYZt.request().ptr),
                        stiffness,
                        npClstrInd.data(), npClstrInd.size(),
                        npClstr.data(),    npClstr.size(),
                        npXYZ.data(),      npXYZ.shape()[0]);
}

void PyPBD_ConstProj_Cloth
(py::array_t<double>& npXYZt,
 const CMeshDynTri2D& mesh)
{
  double* aXYZt = (double*)(npXYZt.request().ptr);
  const std::vector<ETri>& aETri = mesh.aETri;
  const std::vector<CVector2>& aVec2 = mesh.aVec2;
  for(unsigned int it=0;it<aETri.size();++it){
    const int i0 = aETri[it].v[0];
    const int i1 = aETri[it].v[1];
    const int i2 = aETri[it].v[2];
    const double P[3][2] = {
      {aVec2[i0].x,aVec2[i0].y},
      {aVec2[i1].x,aVec2[i1].y},
      {aVec2[i2].x,aVec2[i2].y} };
    double p[3][3] = {
      {aXYZt[i0*3+0], aXYZt[i0*3+1], aXYZt[i0*3+2] },
      {aXYZt[i1*3+0], aXYZt[i1*3+1], aXYZt[i1*3+2] },
      {aXYZt[i2*3+0], aXYZt[i2*3+1], aXYZt[i2*3+2] },
    };
    double C[3], dCdp[3][9];  PBD_CdC_TriStrain2D3D(C, dCdp, P, p);
    double m[3] = {1,1,1};
    const int aIP[3] = {i0,i1,i2};
    PBD_Update_Const3(aXYZt, 3, 3, m, C, &dCdp[0][0], aIP);
  }
}

void PyPointFixBC
(py::array_t<double>& aTmp,
 const py::array_t<int>& aBC,
 const py::array_t<double>& npXY1)
{
  assert( aTmp.ndim() == 2 );
  assert( npXY1.ndim() == 2 );
  assert( aTmp.shape()[1] == npXY1.shape()[1] );
  const int np = aTmp.shape()[0];
  double* ptr = (double*)(aTmp.request().ptr);
  if( npXY1.shape()[1] == 2 ){
    for(int ip=0;ip<np;++ip){
      if( aBC.at(ip) == 0 ){ continue; }
      ptr[ip*2+0] = npXY1.at(ip,0);
      ptr[ip*2+1] = npXY1.at(ip,1);
    }
  }
  if( npXY1.shape()[1] == 3 ){
    for(int ip=0;ip<np;++ip){
      if( aBC.at(ip) == 0 ){ continue; }
      ptr[ip*3+0] = npXY1.at(ip,0);
      ptr[ip*3+1] = npXY1.at(ip,1);
      ptr[ip*3+2] = npXY1.at(ip,2);
    }
  }
}

void init_fem(py::module &m){
  py::class_<CMatrixSparse>(m,"CppMatrixSparse")
  .def(py::init<>())
  .def("initialize", &CMatrixSparse::Initialize)
  .def("set_zero",   &CMatrixSparse::SetZero)
  .def("add_dia",    &CMatrixSparse::AddDia);
  
  py::class_<CPreconditionerILU>(m,"PreconditionerILU")
  .def(py::init<>())
  .def("ilu_decomp", &CPreconditionerILU::DoILUDecomp)
  .def("set_value",  &CPreconditionerILU::SetValueILU);
  
  m.def("matrixSquareSparse_setPattern",     &MatrixSquareSparse_SetPattern);
  m.def("matrixSquareSparse_setFixBC",       &MatrixSquareSparse_SetFixBC);
  m.def("matrixSquareSparse_ScaleLeftRight", &MatrixSquareSparse_ScaleLeftRight);
  
  m.def("addMasterSlavePattern",         &PyAddMasterSlavePattern);
  m.def("precond_ilu0",                  &PrecondILU0);
  m.def("linsys_solve_pcg",              &PySolve_PCG);
  m.def("linsys_solve_bicgstab",         &PySolve_PBiCGStab);
  m.def("linearSystem_setMasterSlave",   &LinearSystem_SetMasterSlave);
  m.def("masterSlave_distributeValue",   &PyMasterSlave_DistributeValue);
  
  m.def("fem_merge_poission",          &PyMergeLinSys_Poission);
  m.def("fem_merge_diffuse",           &PyMergeLinSys_Diffuse);
  m.def("fem_merge_linearSolidStatic", &PyMergeLinSys_LinearSolidStatic);
  m.def("fem_merge_linearSolidDynamic",&PyMergeLinSys_LinearSolidDynamic);
  m.def("fem_merge_storksStatic2D",    &PyMergeLinSys_StorksStatic2D);
  m.def("fem_merge_storksDynamic2D",   &PyMergeLinSys_StorksDynamic2D);
  m.def("fem_merge_navierStorks2D",    &PyMergeLinSys_NavierStorks2D);
  m.def("fem_merge_cloth",             &PyMergeLinSys_Cloth);
  m.def("fem_merge_massPoint",         &PyMergeLinSys_MassPoint);
  m.def("fem_merge_contact",           &PyMergeLinSys_Contact);
  
  m.def("pbd_proj_rigid2d",            &PyPBD_ConstProj_Rigid2D);
  m.def("pbd_proj_rigid3d",            &PyConstProj_Rigid3D);
  m.def("pbd_pointFixBC",              &PyPointFixBC);
  m.def("pbd_proj_cloth",              &PyPBD_ConstProj_Cloth);
}
