/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef FEM_H
#define FEM_H

#include <cassert>
#include <vector>

#include "delfem2/matrix_sparse.h"
#include "delfem2/ilu_sparse.h"
#include "delfem2/fem_ematrix.h"

/*
void FetchData(double* val_to,
               int nno, int ndim,
               const int* aIP,
               const double* val_from, int nstride, int noffset);
 */

void XPlusAY(std::vector<double>& X,
             const int nDoF,
             const std::vector<int>& aBCFlag,
             double alpha,
             const std::vector<double>& Y);

void XPlusAYBZ(std::vector<double>& X,
               const int nDoF,
               const std::vector<int>& aBCFlag,
               double alpha,
               const std::vector<double>& Y,
               double beta,
               const std::vector<double>& Z);

void XPlusAYBZCW(std::vector<double>& X,
                 const int nDoF,
                 const std::vector<int>& aBCFlag,
                 double alpha,
                 const std::vector<double>& Y,
                 double beta,
                 const std::vector<double>& Z,
                 double gamma,
                 const std::vector<double>& W);

// set boundary condition
void setRHS_Zero(std::vector<double>& vec_b,
                 const std::vector<int>& aBCFlag,
                 int iflag_nonzero);

void setRHS_MasterSlave(double* vec_b,
                        int nDoF,
                        const int* aMSFlag);

void SolveLinSys_PCG(const CMatrixSparse& mat_A,
                     std::vector<double>& vec_b,
                     std::vector<double>& vec_x,
                     CPreconditionerILU& ilu_A,
                     double& conv_ratio,
                     int& iteration);

bool SolveLinSys_BiCGStab(CMatrixSparse& mat_A,
                          std::vector<double>& vec_b,
                          std::vector<double>& vec_x,
                          CPreconditionerILU& ilu_A,
                          double& conv_ratio,
                          int& iteration);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MergeLinSys_Poission_MeshTri2D(CMatrixSparse& mat_A,
                            double* vec_b,
                            const double alpha,
                            const double source,
                            const double* aXY1, int np,
                            const unsigned int* aTri1, int nTri,
                            const double* aVal);

void MergeLinSys_Poission_MeshTet3D(CMatrixSparse& mat_A,
                            double* vec_b,
                            const double alpha,
                            const double source,
                            const double* aXYZ, int nXYZ,
                            const unsigned int* aTet, int nTet,
                            const double* aVal);

void MergeLinSys_Diffusion_MeshTri2D(CMatrixSparse& mat_A,
                             double* vec_b,
                             const double alpha,
                             const double rho,
                             const double source,
                             const double dt_timestep,
                             const double gamma_newmark,
                             const double* aXY1, int nXY,
                             const unsigned int* aTri1, int nTri,
                             const double* aVal,
                             const double* aVelo);

void MergeLinSys_Diffusion_MeshTet3D(CMatrixSparse& mat_A,
                             double* vec_b,
                             const double alpha,
                             const double rho,
                             const double source,
                             const double dt_timestep,
                             const double gamma_newmark,
                             const double* aXYZ, int nXYZ,
                             const unsigned int* aTet, int nTet,
                             const double* aVal,
                             const double* aVelo);

void MergeLinSys_SolidStaticLinear_MeshTri2D(CMatrixSparse& mat_A,
                                      double* vec_b,
                                      const double myu,
                                      const double lambda,
                                      const double rho,
                                      const double g_x,
                                      const double g_y,
                                      const double* aXY1, int nXY,
                                      const unsigned int* aTri1, int nTri,
                                      const double* aVal);

void MergeLinSys_SolidDynamicLinear_MeshTri2D(CMatrixSparse& mat_A,
                                      double* vec_b,
                                      const double myu,
                                      const double lambda,
                                      const double rho,
                                      const double g_x,
                                      const double g_y,
                                      const double dt_timestep,
                                      const double gamma_newmark,
                                      const double beta_newmark,
                                      const double* aXY1, int nXY1,
                                      const unsigned int* aTri1, int nTri,
                                      const double* aVal,
                                      const double* aVelo,
                                      const double* aAcc);

void MergeLinSys_StokesStatic2D(CMatrixSparse& mat_A,
                                double* vec_b,
                                const double myu,
                                const double g_x,
                                const double g_y,
                                const double* aXY1, int nXY1,
                                const unsigned int* aTri1, int nTri1,
                                const double* aVal);

void MergeLinSys_StokesDynamic2D(CMatrixSparse& mat_A,
                                  double* vec_b,
                                  const double myu,
                                  const double rho,
                                  const double g_x,
                                  const double g_y,
                                  const double dt_timestep,
                                  const double gamma_newmark,
                                  const double* aXY1, int nXY1,
                                  const unsigned int* aTri1, int nTri1,
                                  const double* aVal,
                                  const double* aVelo);

void MergeLinSys_NavierStokes2D(CMatrixSparse& mat_A,
                                double* vec_b,
                                const double myu,
                                const double rho,
                                const double g_x,
                                const double g_y,
                                const double dt_timestep,
                                const double gamma_newmark,
                                const double* aXY1, int nXY,
                                const unsigned int* aTri1, int nTri,
                                const double* aVal,
                                const double* aVelo);

double MergeLinSys_Cloth(CMatrixSparse& mat_A, // (out) second derivative of energy
                         double* vec_b, // (out) first derivative of energy
                         ////
                         double lambda, // (in) Lame's 1st parameter
                         double myu,  // (in) Lame's 2nd parameter
                         double stiff_bend, // (in) bending stiffness
                         const double* aPosIni, int np, int ndim,
                         const unsigned int* aTri, int nTri, // (in) triangle index
                         const unsigned int* aQuad, int nQuad, // (in) index of 4 vertices required for bending
                         const double* aXYZ);

double MergeLinSys_Contact(CMatrixSparse& ddW,
                         double* dW, // (out) first derivative of energy
                         ////
                         double stiff_contact,
                         double contact_clearance,
                         const CInput_Contact& input,
                         const double* aXYZ,  int nXYZ);

/*
void Solve_LinearSolid_TetP1()
{

unsigned int ndof = aXYZ.size();
MatrixXd Mat(ndof,ndof);
//  Eigen::SparseMatrix<double> Mat(ndof,ndof);
Mat.setZero();

VectorXd Lhs(ndof);
Lhs.setZero();

for(unsigned int itet=0;itet<aTet.size()/4;itet++){
const unsigned int aPo[4] = { aTet[itet*4+0], aTet[itet*4+1], aTet[itet*4+2], aTet[itet*4+3] };
unsigned int i0=aPo[0], i1=aPo[1], i2=aPo[2], i3=aPo[3];
double p0[3] = {aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]};
double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
double p2[3] = {aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2]};
double p3[3] = {aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2]};
////
const double vol = TetVolume3D(p0,p1,p2,p3);
double dldx[4][3];
double zero_order_term[4];
TetDlDx(dldx, zero_order_term,   p0,p1,p2,p3);
////
double disp[4][3] = { {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0} };
double emat[4][4][3][3], eres[4][3];
matRes_LinearSolid_TetP1(emat,eres,
vol,lambda,myu,
0,1,0, 1,
dldx,disp);

for(unsigned int ino=0;ino<4;ino++){
for(unsigned int jno=0;jno<4;jno++){
unsigned int ipo = aPo[ino];
unsigned int jpo = aPo[jno];
for(unsigned int idim=0;idim<3;idim++){
for(unsigned int jdim=0;jdim<3;jdim++){
Mat(ipo*3+idim,jpo*3+jdim) += emat[ino][jno][idim][jdim];
//            Mat.coeffRef(ipo*3+idim,jpo*3+jdim) += emat[ino][jno][idim][jdim];
}
}
}
}
for(unsigned int ino=0;ino<4;ino++){
unsigned int ipo = aPo[ino];
for(unsigned int idim=0;idim<3;idim++){
Lhs(ipo*3+idim) += eres[ino][idim];
}
}
}

std::vector<unsigned int> aFix;
{
for(unsigned int ipo=0;ipo<aXYZ.size()/3;ipo++){
double p[3] = {aXYZ[ipo*3+0], aXYZ[ipo*3+1], aXYZ[ipo*3+2]};
if( p[0] < 0 ){
aFix.push_back(ipo);
}
}
}
for(unsigned int ifix=0;ifix<aFix.size();ifix++){
unsigned int ifix0 = aFix[ifix];
for(int i=0;i<3;i++){
Lhs(ifix0*3+i) = 0;
}
for(unsigned int jpo=0;jpo<aXYZ.size()/3;jpo++){
for(int i=0;i<3;i++){
for(int j=0;j<3;j++){
Mat(ifix0*3+i,jpo*3+j) = 0;
//          Mat.coeffRef(ifix0*3+i,jpo*3+j) = 0;
}
}
}
for(unsigned int jpo=0;jpo<aXYZ.size()/3;jpo++){
for(int i=0;i<3;i++){
for(int j=0;j<3;j++){
Mat(jpo*3+i,ifix0*3+j) = 0;
//          Mat.coeffRef(jpo*3+i,ifix0*3+j) = 0;
}
}
}
{
for(int i=0;i<3;i++){
Mat(ifix0*3+i,ifix0*3+i) = 1;
//        Mat.coeffRef(ifix0*3+i,ifix0*3+i) = 1;
}
}
}


//  sMat = Mat;
//  Mat.makeCompressed();

VectorXd Rhs;
Rhs = Mat.partialPivLu().solve(Lhs);

//  std::cout << ndof << std::endl;
//  ConjugateGradient<SparseMatrix<double>> cg;
//  cg.setMaxIterations(10000);
//  cg.compute(Mat);
//  Rhs = cg.solve(Lhs);

for(unsigned int ip=0;ip<aXYZ.size()/3;ip++){
for(unsigned int i=0;i<3;i++){
aDisp[ip*3+i] += Rhs(ip*3+i);
}
}
}
*/

void MergeLinSys_SolidStaticLinear_MeshTet3D(CMatrixSparse& mat_A,
                                         double* vec_b,
                                         const double myu,
                                         const double lambda,
                                         const double rho,
                                         const double g_x,
                                         const double g_y,
                                         const double g_z,
                                         const double* aXYZ, int nXYZ,
                                         const unsigned int* aTet, int nTet,
                                         const double* aVal);

void MergeLinSys_LinearSolid3D_Static_Q1(CMatrixSparse& mat_A,
                                         std::vector<double>& vec_b,
                                         const double myu,
                                         const double lambda,
                                         const double rho,
                                         const double g_x,
                                         const double g_y,
                                         const double g_z,
                                         const std::vector<double>& aXYZ,
                                         const std::vector<int>& aHex,
                                         const std::vector<double>& aVal);

void MergeLinSys_SolidDynamicLinear_MeshTet3D(CMatrixSparse& mat_A,
                                       double* vec_b,
                                       const double myu,
                                       const double lambda,
                                       const double rho,
                                       const double g_x,
                                       const double g_y,
                                       const double g_z,
                                       const double dt_timestep,
                                       const double gamma_newmark,
                                       const double beta_newmark,
                                       const double* aXYZ, int nXYZ,
                                       const unsigned int* aTet, int nTet,
                                       const double* aVal,
                                       const double* aVelo,
                                       const double* aAcc);

void MergeLinSys_Stokes3D_Static(CMatrixSparse& mat_A,
                                 std::vector<double>& vec_b,
                                 const double myu,
                                 const double rho,
                                 const double g_x,
                                 const double g_y,
                                 const double g_z,
                                 const std::vector<double>& aXYZ,
                                 const std::vector<unsigned int>& aTet,
                                 const std::vector<double>& aVal,
                                 const std::vector<double>& aVelo);

void MergeLinSys_Stokes3D_Dynamic(CMatrixSparse& mat_A,
                                  std::vector<double>& vec_b,
                                  const double myu,
                                  const double rho,
                                  const double g_x,
                                  const double g_y,
                                  const double g_z,
                                  const double dt_timestep,
                                  const double gamma_newmark,
                                  const std::vector<double>& aXYZ,
                                  const std::vector<unsigned int>& aTet,
                                  const std::vector<double>& aVal,
                                  const std::vector<double>& aVelo);

void MergeLinSys_NavierStokes3D_Dynamic(CMatrixSparse& mat_A,
                                        std::vector<double>& vec_b,
                                        const double myu,
                                        const double rho,
                                        const double g_x,
                                        const double g_y,
                                        const double g_z,
                                        const double dt_timestep,
                                        const double gamma_newmark,
                                        const std::vector<double>& aXYZ,
                                        const std::vector<unsigned int>& aTet,
                                        const std::vector<double>& aVal,
                                        const std::vector<double>& aVelo);

/* ------------------------------------------------------------------------- */



/*
void Solve_LinearSolid_TetP2()
{

unsigned int ncorner = aXYZ.size()/3;
unsigned int nedge = aEdge.size()/2;
unsigned int npoint = ncorner + nedge;
unsigned int ndof = npoint*3;
std::cout << "ndof: " << ndof << std::endl;
MatrixXd Mat(ndof,ndof);
//  Eigen::SparseMatrix<double> Mat(ndof,ndof);
//  Mat.reserve(VectorXi::Constant(ndof,100));
Mat.setZero();

VectorXd Lhs(ndof);
Lhs.setZero();

double lambda = 0;
double myu = 1;

for(unsigned int itet=0;itet<aTet.size()/4;itet++){
const unsigned int aPo[10] = {
aTet[itet*4+0],
aTet[itet*4+1],
aTet[itet*4+2],
aTet[itet*4+3],
aTetEdge[itet*6+0] + ncorner,
aTetEdge[itet*6+1] + ncorner,
aTetEdge[itet*6+2] + ncorner,
aTetEdge[itet*6+3] + ncorner,
aTetEdge[itet*6+4] + ncorner,
aTetEdge[itet*6+5] + ncorner };
unsigned int i0=aPo[0], i1=aPo[1], i2=aPo[2], i3=aPo[3];
double p0[3] = {aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]};
double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
double p2[3] = {aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2]};
double p3[3] = {aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2]};
////
const double vol = TetVolume3D(p0,p1,p2,p3);
double dldx[4][3];
double zero_order_term[4];
TetDlDx(dldx, zero_order_term,   p0,p1,p2,p3);
////
double disp[10][3];
for(unsigned int ino=0;ino<10;ino++){
disp[ino][0] = 0;
disp[ino][1] = 0;
disp[ino][2] = 0;
}
double emat[10][10][3][3], eres[10][3];
matRes_LinearSolid_TetP2(emat,eres,
vol,lambda,myu,
0,1,0, 1,
dldx,disp);

for(unsigned int ino=0;ino<10;ino++){
for(unsigned int jno=0;jno<10;jno++){
unsigned int ipo = aPo[ino];
unsigned int jpo = aPo[jno];
for(unsigned int idim=0;idim<3;idim++){
for(unsigned int jdim=0;jdim<3;jdim++){
Mat(ipo*3+idim,jpo*3+jdim) += emat[ino][jno][idim][jdim];
//        Mat.coeffRef(ipo*3+idim, jpo*3+jdim) += emat[ino][jno][idim][jdim];
}
}
}
}
for(unsigned int ino=0;ino<10;ino++){
unsigned int ipo = aPo[ino];
for(unsigned int idim=0;idim<3;idim++){
Lhs(ipo*3+idim) += eres[ino][idim];
}
}
//    std::cout << itet << " " << aTet.size()/4 << std::endl;
}

std::vector<unsigned int> aFix;
{
for(unsigned int ipo=0;ipo<aXYZ.size()/3;ipo++){
double p[3] = {aXYZ[ipo*3+0], aXYZ[ipo*3+1], aXYZ[ipo*3+2]};
if( p[0] < 0 ){
aFix.push_back(ipo);
}
}
for(unsigned int iedge=0;iedge<aEdge.size()/2;iedge++){
unsigned int i0 = aEdge[iedge*2+0];
unsigned int i1 = aEdge[iedge*2+1];
double p0[3] = {aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]};
double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
if( p0[0] < 0 && p1[0] < 0 ){
aFix.push_back(iedge+ncorner);
}
}
}
for(unsigned int iifix=0;iifix<aFix.size();iifix++){
unsigned int ifix0 = aFix[iifix];
for(int i=0;i<3;i++){
Lhs(ifix0*3+i) = 0;
}
for(unsigned int jpo=0;jpo<npoint;jpo++){
for(int i=0;i<3;i++){
for(int j=0;j<3;j++){
Mat(ifix0*3+i,jpo*3+j) = 0;
//          Mat.coeffRef(ifix0*3+i,jpo*3+j) = 0;
}
}
}
for(unsigned int jpo=0;jpo<npoint;jpo++){
for(int i=0;i<3;i++){
for(int j=0;j<3;j++){
Mat(jpo*3+i,ifix0*3+j) = 0;
//          Mat.coeffRef(jpo*3+i,ifix0*3+j) = 0;
}
}
}
{
for(int i=0;i<3;i++){
Mat(ifix0*3+i,ifix0*3+i) = 1;
//        Mat.coeffRef(ifix0*3+i,ifix0*3+i) = 1;
}
}
}
//  Mat.makeCompressed();

VectorXd Rhs;
Rhs = Mat.partialPivLu().solve(Lhs);

//  std::cout << ndof << std::endl;
//  ConjugateGradient<MatrixXd> cg;
//  cg.setMaxIterations(3000);
//  cg.compute(Mat);
//  Rhs = cg.solve(Lhs);
for(unsigned int ip=0;ip<aXYZ.size()/3;ip++){
for(unsigned int i=0;i<3;i++){
aDisp[ip*3+i] += Rhs(ip*3+i);
}
std::cout << aDisp[ip*3+0] << " " << aDisp[ip*3+1] << " " << aDisp[ip*3+2] << std::endl;
}
}
*/

/*
// this was used
void Solve_LinearSolid_TetP2_MyMat()
{
unsigned int ncorner = aXYZ.size()/3;
unsigned int nedge = aEdge.size()/2;
unsigned int npoint = ncorner + nedge;

////
MatVec::CMatDia_BlkCrs matFEM;
MatVec::CVector_Blk vecRes;
MatVec::CVector_Blk vecUpd;
MatVec::CBCFlag bc_flag(npoint,3);
{
std::vector<unsigned int> aEl;
aEl.resize(aTet.size()/4*10);
for(unsigned int itet=0;itet<aTet.size()/4;itet++){
const unsigned int aPo[10] = {
aTet[itet*4+0],
aTet[itet*4+1],
aTet[itet*4+2],
aTet[itet*4+3],
aTetEdge[itet*6+0] + ncorner,
aTetEdge[itet*6+1] + ncorner,
aTetEdge[itet*6+2] + ncorner,
aTetEdge[itet*6+3] + ncorner,
aTetEdge[itet*6+4] + ncorner,
aTetEdge[itet*6+5] + ncorner };
for(unsigned int ino=0;ino<10;ino++){
aEl[itet*10+ino] = aPo[ino];
}
}
matFEM.Initialize(npoint, 3);
CIndexedArray crs_edge;
crs_edge.SetEdge(npoint, aTet.size()/4, 10, aEl);
matFEM.AddPattern(crs_edge);
}
vecRes.Initialize(npoint,3);
vecUpd.Initialize(npoint,3);

matFEM.SetZero();
vecRes.SetVectorZero();

unsigned int ndof = npoint*3;
std::cout << "ndof: " << ndof << std::endl;

for(unsigned int itet=0;itet<aTet.size()/4;itet++){
const unsigned int aPo[10] = {
aTet[itet*4+0],
aTet[itet*4+1],
aTet[itet*4+2],
aTet[itet*4+3],
aTetEdge[itet*6+0] + ncorner,
aTetEdge[itet*6+1] + ncorner,
aTetEdge[itet*6+2] + ncorner,
aTetEdge[itet*6+3] + ncorner,
aTetEdge[itet*6+4] + ncorner,
aTetEdge[itet*6+5] + ncorner };
unsigned int i0=aPo[0], i1=aPo[1], i2=aPo[2], i3=aPo[3];
double p0[3] = {aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]};
double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
double p2[3] = {aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2]};
double p3[3] = {aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2]};
////
const double vol = TetVolume3D(p0,p1,p2,p3);
double dldx[4][3];
double zero_order_term[4];
TetDlDx(dldx, zero_order_term,   p0,p1,p2,p3);
////
double disp[10][3];
for(unsigned int i=0;i<3;i++){
disp[0][i] = aDisp[ aPo[0]*3+i];
disp[1][i] = aDisp[ aPo[1]*3+i];
disp[2][i] = aDisp[ aPo[2]*3+i];
disp[3][i] = aDisp[ aPo[3]*3+i];
disp[4][i] = aDispEdge[ aTetEdge[itet*6+0]*3+i ];
disp[5][i] = aDispEdge[ aTetEdge[itet*6+1]*3+i ];
disp[6][i] = aDispEdge[ aTetEdge[itet*6+2]*3+i ];
disp[7][i] = aDispEdge[ aTetEdge[itet*6+3]*3+i ];
disp[8][i] = aDispEdge[ aTetEdge[itet*6+4]*3+i ];
disp[9][i] = aDispEdge[ aTetEdge[itet*6+5]*3+i ];
}
////
double emat[10][10][3][3], eres[10][3];
matRes_LinearSolid_TetP2(emat,eres,
vol,lambda,myu,
g_x,g_y,g_z, 1,
dldx,disp);

matFEM.Mearge(10, aPo, 10, aPo, 9, &emat[0][0][0][0]);
for(unsigned int ino=0;ino<10;ino++){
unsigned int ipo = aPo[ino];
for(unsigned int idim=0;idim<3;idim++){
vecRes.AddValue(ipo,idim,eres[ino][idim]);
}
}
}
{
vecRes.AddValue(ino_exforce,0,exforce_vector.x);
vecRes.AddValue(ino_exforce,1,exforce_vector.y);
vecRes.AddValue(ino_exforce,2,exforce_vector.z);
}

{
for(unsigned int ipo=0;ipo<aXYZ.size()/3;ipo++){
double p[3] = {aXYZ[ipo*3+0], aXYZ[ipo*3+1], aXYZ[ipo*3+2]};
double h = p[0]*norm.x + p[1]*norm.y + p[2]*norm.z;
if( h < -clearance_height ){
bc_flag.SetBC(ipo, 0);
bc_flag.SetBC(ipo, 1);
bc_flag.SetBC(ipo, 2);
}
}
for(unsigned int iedge=0;iedge<aEdge.size()/2;iedge++){
unsigned int i0 = aEdge[iedge*2+0];
unsigned int i1 = aEdge[iedge*2+1];
double p0[3] = {aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]};
double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
double h0 = p0[0]*norm.x + p0[1]*norm.y + p0[2]*norm.z;
double h1 = p1[0]*norm.x + p1[1]*norm.y + p1[2]*norm.z;
if( p0[0] < -clearance_height && p1[0] < -clearance_height ){
bc_flag.SetBC(iedge+ncorner,0);
bc_flag.SetBC(iedge+ncorner,1);
bc_flag.SetBC(iedge+ncorner,2);
}
}
matFEM.SetBoundaryCondition(bc_flag);
bc_flag.SetZeroToBCDof(vecRes);
}

double conv_ratio = 1.0e-5;
unsigned int iteration = 1000;
MatVec::Sol::Solve_CG(conv_ratio, iteration, matFEM, vecRes, vecUpd);
//  MatVec::Sol::Solve_PCG(conv_ratio, iteration, matFEM, vecRes, vecUpd, matPrec,matFEM);
for(unsigned int ipo=0;ipo<aXYZ.size()/3;ipo++){
aDisp[ipo*3+0] += vecUpd.GetValue(ipo, 0);
aDisp[ipo*3+1] += vecUpd.GetValue(ipo, 1);
aDisp[ipo*3+2] += vecUpd.GetValue(ipo, 2);
}
for(unsigned int ipo=0;ipo<aEdge.size()/2;ipo++){
aDispEdge[ipo*3+0] = vecUpd.GetValue(ipo+ncorner,0);
aDispEdge[ipo*3+1] = vecUpd.GetValue(ipo+ncorner,1);
aDispEdge[ipo*3+2] = vecUpd.GetValue(ipo+ncorner,2);
}

}

*/


#endif /* fem_utility_h */
