#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits>
#include <vector>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/mshtopo.h"
#include "delfem2/msh.h"
#include "delfem2/dyntri_v3.h"

#include "delfem2/ilu_sparse.h"
#include "delfem2/fem.h"

#include "delfem2/color_gl.h"
#include "delfem2/funcs_gl.h"
#include "delfem2/funcs_glut.h"

static double TriArea2D(const double v1[], const double v2[], const double v3[]){
  return 0.5*( (v2[0]-v1[0])*(v3[1]-v1[1]) - (v3[0]-v1[0])*(v2[1]-v1[1]) );
}

double AreaCGCurve(const std::vector<double>& aCV, double cg[2])
{
  const unsigned int nCV = (int)aCV.size()/2;
  double area = 0;
  cg[0] = 0;
  cg[1] = 0;
  for(unsigned int idiv=0;idiv<nCV;idiv++){
    unsigned int ipo0 = idiv;
    unsigned int ipo1 = idiv+1;
    if( idiv == nCV-1 ){ ipo1 = 0; }
    const double p0[2] = { aCV[ipo0*2+0], aCV[ipo0*2+1] };
    const double p1[2] = { aCV[ipo1*2+0], aCV[ipo1*2+1] };
    const double p2[2] = { 0, 0 };
    double a0 = TriArea2D(p0, p1, p2);
    double cg0[2] = { (p0[0]+p1[0]+p2[0])/3.0, (p0[1]+p1[1]+p2[1])/3.0 };
    cg[0] += cg0[0]*a0;
    cg[1] += cg0[1]*a0;
    area += a0;
  }
  cg[0] /= area;
  cg[1] /= area;
  return area;
}

void MakeRandomCV(unsigned int nCV, std::vector<double>& aCV)
{
  aCV.clear();
  for(unsigned int icv=0;icv<nCV;icv++){
    /*
     {
     aCV.push_back((double)rand()/(RAND_MAX+1.0));
     aCV.push_back((double)rand()/(RAND_MAX+1.0));
     }
     */
    {
      double tht = icv*3.1415*2.0/nCV;
      double r = (double)rand()/(RAND_MAX+1.0);
      double px = r*sin(tht);
      double py = r*cos(tht);
      aCV.push_back(px);
      aCV.push_back(py);
    }
  }
  
  {
    double cnt[2];
    double area = AreaCGCurve(aCV,cnt);
    if( area < 0 ){
      std::vector<double> aCV0;
      for(unsigned int idiv=0;idiv<nCV;idiv++){
        unsigned int idiv0 = nCV-1-idiv;
        aCV0.push_back( aCV[idiv0*2+0] );
        aCV0.push_back( aCV[idiv0*2+1] );
      }
      aCV = aCV0;
      area = -area;
    }
    if( area < 1.0e-5 ){
      aCV.clear();
      return;
    }
  }
}


void MakeCurveSpline(const std::vector<double>& aCV, std::vector<double>& aVecCurve, int ndiv=5)
{
  aVecCurve.resize(0);
  const unsigned int nCV = aCV.size()/2;
  for(unsigned int icv=0;icv<nCV;icv++){
    int icv0=(icv+0)%nCV;
    int icv1=(icv+1)%nCV;
    int icv2=(icv+2)%nCV;
    const double p0[2] = { aCV[icv0*2+0], aCV[icv0*2+1] };
    const double p1[2] = { aCV[icv1*2+0], aCV[icv1*2+1] };
    const double p2[2] = { aCV[icv2*2+0], aCV[icv2*2+1] };
    for(unsigned int idiv=0;idiv<ndiv;idiv++){
      const double t = 1.0-(double)idiv/ndiv;
      const double w[3] = {0.5*t*t, -t*t + t + 0.5, 0.5*(1-t)*(1-t) };
      const double px = w[0]*p0[0] + w[1]*p1[0] + w[2]*p2[0];
      const double py = w[0]*p0[1] + w[1]*p1[1] + w[2]*p2[1];
      aVecCurve.push_back(px);
      aVecCurve.push_back(py);
    }
  }
}

void myGlVertex2D(const std::vector<double>& vec, unsigned int i)
{
  ::glVertex3d(vec[i*2],vec[i*2+1],+0.5);
}

void drawCurve
(const std::vector<double>& vec,
 const std::vector<double>& aVecCurve0)
{
  if( aVecCurve0.size() < 2 ){ return; }
  ::glBegin(GL_LINES);
  const unsigned int nvec = (int)vec.size()/2;
  for(unsigned int ivec=0;ivec<nvec;ivec++){
    unsigned int jvec = ivec+1; if( jvec >= nvec ){ jvec -= nvec; }
    myGlVertex2D(vec,ivec);
    myGlVertex2D(vec,jvec);
  }
  ::glEnd();
  
  ::glBegin(GL_POINTS);
  for(unsigned int ivec=0;ivec<nvec;ivec++){
    myGlVertex2D(vec,ivec);
  }
  ::glEnd();
}



////////////////////////////////////////////////////////////////////////////////////

double len = 1.1;
std::vector<double> aVecCurve0; // current test

std::vector<int> aTri1;
std::vector<double> aXY1;
std::vector<int> aInd_IdVtxLoop, aIdVtxLoop; // vtx on loop

std::vector<double> aVal;
std::vector<double> aVelo;
std::vector<double> aAcc;
std::vector<int> aBCFlag; // boundary condition flag
std::vector<int> aMSFlag; // master slave flag

CMatrixSquareSparse mat_A;
std::vector<double> vec_b;
CPreconditionerILU  ilu_A;

int press_button = -1;
double mov_begin_x, mov_begin_y;
double mag = 1.5;

// 0: poisson
// 1: diffusion
int iproblem = 0;
bool is_animation = false;
double dt_timestep = 0.01;
double gamma_newmark = 0.6;
double beta_newmark = 0.36;

//////////////////////////////////////////////////////////////////////////////////////

void MakeMesh(){
  std::vector<double> aCV0; MakeRandomCV(8, aCV0); // current cv
  MakeCurveSpline(aCV0, aVecCurve0,20); // current curve
  std::vector< std::vector<double> > aVecAry;
  {
    aVecAry.resize(1);
    aVecAry[0].push_back(-len); aVecAry[0].push_back(-len);
    aVecAry[0].push_back(+len); aVecAry[0].push_back(-len);
    aVecAry[0].push_back(+len); aVecAry[0].push_back(+len);
    aVecAry[0].push_back(-len); aVecAry[0].push_back(+len);
  }
  aVecAry.push_back(aVecCurve0);
  //  bool res = GenerateTesselation2(aTri1, aXY1, aInd_IdVtxLoop,aIdVtxLoop, 0.045, true,aVecAry);
  bool res = GenerateTesselation2(aTri1, aXY1, aInd_IdVtxLoop, aIdVtxLoop, 0.05, true, aVecAry);
  if( !res ){
    std::cout << "meshing failed" << std::endl;
  }
  std::cout<<"  ntri;"<<aTri1.size()/3<<"  nXY:"<<aXY1.size()/2<<std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
// iproblem: 0, 1
void InitializeProblem_Scalar()
{
  const int np = (int)aXY1.size()/2;
  aBCFlag.assign(np, 0);
  for(int ip=0;ip<np;++ip){
    const double px = aXY1[ip*2+0];
    const double py = aXY1[ip*2+1];
    if( fabs(fabs(px)-len) < 0.0001 || fabs(fabs(py)-len) < 0.0001 ){
      aBCFlag[ip] = 1;
    }
  }
  for(int iip=aInd_IdVtxLoop[1];iip<aInd_IdVtxLoop[2];++iip){
    int ip0 = aIdVtxLoop[iip];
    aBCFlag[ip0] = 1;
  }
  //////
  /*
   CJaggedArray crs;
   crs.SetEdgeOfElem(aTri1, (int)aTri1.size()/3, 3, (int)aXY1.size()/2, false);
   crs.Sort();
   */
  std::vector<int> psup_ind, psup;
  JaggedArray_MeshOneRingNeighborhood(psup_ind, psup,
                                      aTri1.data(), aTri1.size()/3, 3, (int)aXY1.size()/2);
  JaggedArray_Sort(psup_ind, psup);
  ////
  mat_A.Initialize(np, 1, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(),psup.size());
  ilu_A.Initialize_ILU0(mat_A);
}

//////////////////////////////////////////////////////////////////////////////////////////
// iproblem: 0
void SolveProblem_Poisson()
{
  const int np = (int)aXY1.size()/2;
  const int nDoF = np;
  ///////////////////////////
  const double alpha = 1.0;
  const double source = 1.0;
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  MergeLinSys_Poission_MeshTri2D(mat_A,vec_b.data(),
                         alpha,source,
                         aXY1.data(),aXY1.size()/2,
                         aTri1.data(),aTri1.size()/3,
                         aVal.data());
  mat_A.SetBoundaryCondition(aBCFlag.data(),aBCFlag.size(),1);
  setRHS_Zero(vec_b, aBCFlag,0);
  ///////////////////////////
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  ///////////////////////////
  XPlusAY(aVal,nDoF,aBCFlag,
          1.0,vec_x);
}

//////////////////////////////////////////////////////////////////////////////////////////
// iproblem: 1
void SolveProblem_Diffusion()
{
  const int np = (int)aXY1.size()/2;
  const int nDoF = np;
  //////////////////////////
  const double alpha = 1.0;
  const double rho = 1.0;
  const double source = 1.0;
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  MergeLinSys_Diffusion_MeshTri2D(mat_A, vec_b.data(),
                                  alpha, rho, source,
                                  dt_timestep, gamma_newmark,
                                  aXY1.data(), aXY1.size()/2,
                                  aTri1.data(), aTri1.size()/3,
                                  aVal.data(),aVelo.data());
  mat_A.SetBoundaryCondition(aBCFlag.data(),aBCFlag.size(),1);
  setRHS_Zero(vec_b, aBCFlag,0);
  ///////////////////
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  ///////////////////
  XPlusAYBZ(aVal,nDoF,aBCFlag,
            dt_timestep*gamma_newmark,vec_x,
            dt_timestep,aVelo);
  XPlusAY(aVelo,nDoF,aBCFlag,
          1.0,vec_x);
}

/////////////////////////////////////////////////////////////////////
// iproblem: 2, 3
void InitializeProblem_Solid()
{
  const int np = (int)aXY1.size()/2;
  const int nDoF = np*2;
  /////
  aBCFlag.assign(nDoF, 0);
  for(int ip=0;ip<np;++ip){
    const double px = aXY1[ip*2+0];
    const double py = aXY1[ip*2+1];
    if( fabs(py-len) > 0.0001 ){ continue; }
    aBCFlag[ip*2+0] = 1;
    aBCFlag[ip*2+1] = 1;
  }
  aMSFlag.assign(nDoF, -1);
  { // master slave
    int iseed = -1;
    for(int ip=0;ip<np;++ip){
//      const double px = aXY1[ip*2+0];
      const double py = aXY1[ip*2+1];
      if( fabs(py+len) > 0.0001 ){ continue; }
      if( iseed == -1 ){
        iseed = ip;
      }
      else{
        aMSFlag[ip*2+0] = iseed*2+0;
        aMSFlag[ip*2+1] = iseed*2+1;
      }
    }
  }
  //////
  std::vector<int> psup_ind0, psup0;
  JaggedArray_MeshOneRingNeighborhood(psup_ind0, psup0,
                                      aTri1.data(), aTri1.size()/3, 3, (int)aXY1.size()/2);
  std::vector<int> psup_ind, psup;
  addMasterSlavePattern(psup_ind, psup,
                        aMSFlag.data(),2,
                        psup_ind0.data(), psup_ind0.size(), psup0.data());
  JaggedArray_Sort(psup_ind, psup);
  /*
   CJaggedArray crs;
   crs.SetEdgeOfElem(aTri1, (int)aTri1.size()/3, 3, (int)aXY1.size()/2, false);
   crs.addMasterSlavePattern(aMSFlag,2);
   crs.Sort();
   */
  ////
  mat_A.Initialize(np, 2, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(), psup.size());
  ilu_A.Initialize_ILU0(mat_A);
}

/////////////////////////////////////////////////////////////////////
// iproblem: 2
void SolveProblem_LinearSolid_Static()
{
  const int np = (int)aXY1.size()/2;
  const int nDoF = np*2;
  //////////////////////////
  double myu = 10.0;
  double lambda = 10.0;
  double rho = 1.0;
  double g_x = 0.0;
  double g_y = -3.0;
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  MergeLinSys_SolidStaticLinear_MeshTri2D(mat_A,vec_b.data(),
                                          myu,lambda,rho,g_x,g_y,
                                          aXY1.data(), aXY1.size()/2,
                                          aTri1.data(), aTri1.size()/3,
                                          aVal.data());
  mat_A.SetBoundaryCondition(aBCFlag.data(),aBCFlag.size()/2,2);
  setRHS_Zero(vec_b, aBCFlag,0);
  mat_A.SetMasterSlave(aMSFlag.data());
  setRHS_MasterSlave(vec_b.data(),vec_b.size(),aMSFlag.data());
  //////////////////////////
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  //////////////////////////
  XPlusAY(aVal,nDoF,aBCFlag,
          1.0,vec_x);
  for(int idof=0;idof<nDoF;++idof){
    int jdof = aMSFlag[idof];
    if( jdof == -1 ) continue;
    aVal[idof] = aVal[jdof];
  }
}

/////////////////////////////////////////////////////////////////////
// iproblem: 3
void SolveProblem_LinearSolid_Dynamic()
{
  const int np = (int)aXY1.size()/2;
  const int nDoF = np*2;
  //////////////////////////////
  double myu = 10.0;
  double lambda = 10.0;
  double rho = 1.0;
  double g_x = 0.0;
  double g_y = -3.0;
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  MergeLinSys_SolidDynamicLinear_MeshTri2D(mat_A,vec_b.data(),
                                    myu,lambda,rho,g_x,g_y,
                                    dt_timestep,gamma_newmark,beta_newmark,
                                    aXY1.data(), aXY1.size()/2,
                                    aTri1.data(), aTri1.size()/3,
                                    aVal.data(),aVelo.data(),aAcc.data());
  mat_A.SetBoundaryCondition(aBCFlag.data(),aBCFlag.size()/2,2);
  setRHS_Zero(vec_b, aBCFlag,0);
  mat_A.SetMasterSlave(aMSFlag.data());
  setRHS_MasterSlave(vec_b.data(),vec_b.size(),aMSFlag.data());
  //////////////////////////////
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  //////////////////////////////
  XPlusAYBZCW(aVal, nDoF, aBCFlag,
              dt_timestep,aVelo,
              0.5*dt_timestep*dt_timestep,aAcc,
              dt_timestep*dt_timestep*beta_newmark,vec_x);
  XPlusAYBZ(aVelo,nDoF, aBCFlag,
            dt_timestep*gamma_newmark,vec_x,
            dt_timestep,aAcc);
  XPlusAY(aAcc, nDoF, aBCFlag,
          1.0, vec_x);
  for(int idof=0;idof<nDoF;++idof){
    int jdof = aMSFlag[idof];
    if( jdof == -1 ) continue;
    aVal[ idof] = aVal[ jdof];
    aVelo[idof] = aVelo[jdof];
    aAcc[ idof] = aAcc[ jdof];
  }
}

/////////////////////////////////////////////////////////////////////
// iproblem: 4, 5, 6
void InitializeProblem_Fluid()
{
  // set boundary condition
  const int np = (int)aXY1.size()/2;
  const int nDoF = np*3;
  ///////
  aBCFlag.assign(nDoF, 0);
  for(int ip=0;ip<np;++ip){
    const double px = aXY1[ip*2+0];
    const double py = aXY1[ip*2+1];
    if( fabs(fabs(px)-len) < 0.0001 || fabs(fabs(py)-len) < 0.0001 ){
      aBCFlag[ip*3+0] = 1;
      aBCFlag[ip*3+1] = 1;
    }
  }
  for(int iip=aInd_IdVtxLoop[1];iip<aInd_IdVtxLoop[2];++iip){
    int ip0 = aIdVtxLoop[iip];
    aBCFlag[ip0*3+0] = 1;
    aBCFlag[ip0*3+1] = 1;
  }
  aBCFlag[0*3+2] = 1;
  //////
  aMSFlag.clear();
  //////
  for(int ip=0;ip<np;++ip){
//    const double px = aXY1[ip*2+0];
    const double py = aXY1[ip*2+1];
    if( fabs(py-len) < 0.0001 && aBCFlag[ip*3+0] == 1 ){
      aVal[ip*3+0] = 10;
    }
  }
  //////
  std::vector<int> psup_ind, psup;
  JaggedArray_MeshOneRingNeighborhood(psup_ind, psup,
                                      aTri1.data(), aTri1.size()/3, 3, (int)aXY1.size()/2);
  JaggedArray_Sort(psup_ind, psup);
  /*
   CJaggedArray crs;
   crs.SetEdgeOfElem(aTri1, (int)aTri1.size()/3, 3, (int)aXY1.size()/2, false);
   crs.Sort();
   */
  ////
  mat_A.Initialize(np, 3, true);
  mat_A.SetPattern(psup_ind.data(),psup_ind.size(), psup.data(),psup.size());
  ilu_A.Initialize_ILU0(mat_A);
  //  ilu_A.Initialize_ILUk(mat_A, 5);
}

/////////////////////////////////////////////////////////////////////
// iproblem: 7,8,9
void InitializeProblem_Fluid2()
{
  const int np = (int)aXY1.size()/2;
  const int nDoF = np*3;
  // set boundary condition
  aBCFlag.assign(nDoF, 0);
  for(int ip=0;ip<np;++ip){
    const double px = aXY1[ip*2+0];
    const double py = aXY1[ip*2+1];
    if( fabs(fabs(py)-len) < 0.0001 ){
      aBCFlag[ip*3+0] = 1;
      aBCFlag[ip*3+1] = 1;
      aVal[ip*3+0] = 1;
    }
    if( fabs(px-len) < 0.0001 ){
      aBCFlag[ip*3+0] = 1;
      aBCFlag[ip*3+1] = 1;
      aVal[ip*3+0] = 1;
    }
  }
  for(int iip=aInd_IdVtxLoop[1];iip<aInd_IdVtxLoop[2];++iip){
    int ip0 = aIdVtxLoop[iip];
    aBCFlag[ip0*3+0] = 1;
    aBCFlag[ip0*3+1] = 1;
  }
  //  aBCFlag[0*3+2] = 1;
  //////
  aMSFlag.assign(nDoF, -1);
  { // master slave
    int iseed = -1;
    for(int ip=0;ip<np;++ip){
      const double px = aXY1[ip*2+0];
//      const double py = aXY1[ip*2+1];
      if( fabs(px+len) > 0.0001 || aBCFlag[ip*3+0] == 1 ){ continue; }
      if( iseed == -1 ){
        iseed = ip;
      }
      else{
        aMSFlag[ip*3+0] = iseed*3+0;
        aMSFlag[ip*3+1] = iseed*3+1;
        aMSFlag[ip*3+2] = iseed*3+2;
      }
    }
  }
  ///////
  std::vector<int> psup_ind0, psup0;
  JaggedArray_MeshOneRingNeighborhood(psup_ind0, psup0,
                                      aTri1.data(), aTri1.size()/3, 3, (int)aXY1.size()/2);
  std::vector<int> psup_ind, psup;
  addMasterSlavePattern(psup_ind, psup,
                        aMSFlag.data(),3,
                        psup_ind0.data(), psup_ind0.size(), psup0.data());
  JaggedArray_Sort(psup_ind, psup);
  //////
  /*
   CJaggedArray crs;
   crs.SetEdgeOfElem(aTri1, (int)aTri1.size()/3, 3, (int)aXY1.size()/2, false);
   crs.addMasterSlavePattern(aMSFlag,3);
   crs.Sort();
   */
  ////
  mat_A.Initialize(np, 3, true);
  mat_A.SetPattern(psup_ind.data(),psup_ind.size(), psup.data(),psup.size());
  ilu_A.Initialize_ILU0(mat_A);
  //  ilu_A.Initialize_ILUk(mat_A, 5);
}

/////////////////////////////////////////////////////////////////////
// iproblem: 4
void SolveProblem_Stokes_Static()
{
  const int np = (int)aXY1.size()/2;
  const int nDoF = np*3;
  //////////////////////////
  double myu = 1.0;
  double g_x = 0.0;
  double g_y = -0.0;
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  MergeLinSys_StokesStatic2D( mat_A,vec_b.data(),
                              myu,g_x,g_y,
                              aXY1.data(), aXY1.size()/2,
                              aTri1.data(), aTri1.size()/3,
                              aVal.data());
  mat_A.SetBoundaryCondition(aBCFlag.data(),aBCFlag.size()/3,3);
  setRHS_Zero(vec_b, aBCFlag,0);
  if( aMSFlag.size() == vec_b.size() ){
    mat_A.SetMasterSlave(aMSFlag.data());
    setRHS_MasterSlave(vec_b.data(),vec_b.size(),aMSFlag.data());
  }
  //////////////////////////
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  //////////////////////////
  XPlusAY(aVal, nDoF, aBCFlag, 1.0, vec_x);
  if( aMSFlag.size() == nDoF ){
    for(int idof=0;idof<nDoF;++idof){
      int jdof = aMSFlag[idof];
      if( jdof == -1 ) continue;
      assert( jdof >= 0 && jdof < nDoF );
      aVal[ idof] = aVal[ jdof];
    }
  }
}

/////////////////////////////////////////////////////////////////////
// iprob:5
void SolveProblem_Stokes_Dynamic()
{
  const int np = (int)aXY1.size()/2;
  const int nDoF = np*3;
  //////////////////////////
  double myu = 1.0;
  double rho = 10;
  double g_x = 0.0;
  double g_y = -0.0;
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  MergeLinSys_StokesDynamic2D(mat_A,vec_b.data(),
                              myu,rho,g_x,g_y,
                              dt_timestep,gamma_newmark,
                              aXY1.data(), aXY1.size()/2,
                              aTri1.data(), aTri1.size()/3,
                              aVal.data(),aVelo.data());
  mat_A.SetBoundaryCondition(aBCFlag.data(),aBCFlag.size()/3,3);
  setRHS_Zero(vec_b, aBCFlag,0);
  if( aMSFlag.size() == vec_b.size() ){
    mat_A.SetMasterSlave(aMSFlag.data());
    setRHS_MasterSlave(vec_b.data(),vec_b.size(),aMSFlag.data());
  }
  //////////////////////////
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  //////////////////////////
  XPlusAYBZ(aVal,nDoF, aBCFlag,
            dt_timestep*gamma_newmark,vec_x,
            dt_timestep,aVelo);
  XPlusAY(aVelo, nDoF, aBCFlag,
          1.0, vec_x);
  if( aMSFlag.size() == nDoF ){
    for(int idof=0;idof<nDoF;++idof){
      int jdof = aMSFlag[idof];
      if( jdof == -1 ) continue;
      assert( jdof >= 0 && jdof < nDoF );
      aVal[ idof] = aVal[ jdof];
      aVelo[idof] = aVelo[jdof];
    }
  }
}

//////////////////////////////////////////////////////////
// iprob: 6
void SolveProblem_NavierStokes_Dynamic()
{
  const int np = (int)aXY1.size()/2;
  const int nDoF = np*3;
  //////////////////////////////
  double myu = 0.01;
  double rho = 1;
  double g_x = 0.0;
  double g_y = -0.0;
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  MergeLinSys_NavierStokes2D(mat_A,vec_b.data(),
                             myu,rho,g_x,g_y,
                             dt_timestep,gamma_newmark,
                             aXY1.data(), aXY1.size()/2,
                             aTri1.data(), aTri1.size()/3,
                             aVal.data(),aVelo.data());
  mat_A.SetBoundaryCondition(aBCFlag.data(),aBCFlag.size()/3,3);
  setRHS_Zero(vec_b, aBCFlag,0);
  if( aMSFlag.size() == vec_b.size() ){
    mat_A.SetMasterSlave(aMSFlag.data());
    setRHS_MasterSlave(vec_b.data(),vec_b.size(),aMSFlag.data());
  }
  //////////////////////////////
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  SolveLinSys_BiCGStab(mat_A,vec_b,vec_x,ilu_A,
                       conv_ratio, iteration);
  //////////////////////////////
  XPlusAYBZ(aVal,nDoF, aBCFlag,
            dt_timestep*gamma_newmark,vec_x,
            dt_timestep,aVelo);
  XPlusAY(aVelo, nDoF, aBCFlag,
          1.0, vec_x);
  if( aMSFlag.size() == nDoF ){
    for(int idof=0;idof<nDoF;++idof){
      int jdof = aMSFlag[idof];
      if( jdof == -1 ) continue;
      assert( jdof >= 0 && jdof < nDoF );
      aVal[ idof] = aVal[ jdof];
      aVelo[idof] = aVelo[jdof];
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////

void myGlutResize(int w, int h)
{
  glViewport(0, 0, w, h);
  glutPostRedisplay();
}

void myGlutDisplay(void)
{
  ::glClearColor(0.2, .7, 0.7, 1.0);
  //	::glClearColor(0.0, .0, 0.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1, 4.0 );
  
  setGL_Camera2D();
  
  if( iproblem == 0 || iproblem == 1 ){
    {
      std::vector< std::pair<double,CColor> > colorMap;
      makeHeatMap_BlueGrayRed(colorMap, 0, +0.1);
      DrawMeshTri2D_ScalarP1(aXY1.data(),aXY1.size()/2,
                             aTri1.data(),aTri1.size()/3,
                             aVal.data(),1,0,colorMap);
    }
    DrawMeshTri2D_Edge(aTri1,aXY1);
    ::glPointSize(2);
    ::glColor3d(0,0,0);
    DrawPoints2D_Points(aXY1);
  }
  else if( iproblem == 2 || iproblem == 3 ){
    DrawMeshTri2D_FaceDisp2D(aXY1.data(), aXY1.size()/2,
                             aTri1.data(), aTri1.size()/3,
                             aVal.data(), 2);
  }
  else if( iproblem == 4 || iproblem == 5 || iproblem == 6
          || iproblem == 7 || iproblem == 8 || iproblem == 9 )
  {
    std::vector< std::pair<double,CColor> > colorMap;
    makeHeatMap_BlueGrayRed(colorMap, -30, +30);
    DrawMeshTri2D_ScalarP1(aXY1.data(),aXY1.size()/2,
                           aTri1.data(),aTri1.size()/3,
                           aVal.data(),3,2,colorMap);
    ::glColor3d(0,0,0);    
    DrawPoints2D_Vectors(aXY1.data(),aXY1.size()/2, aVal.data(),3,0, 0.1);
    ::glPointSize(2);
    ::glColor3d(0,0,0);
    DrawPoints2D_Points(aXY1);
  }
  
  ShowFPS();
  
  glutSwapBuffers();
}

void myGlutIdle(){
  if( is_animation ){
    if( iproblem == 1){
      SolveProblem_Diffusion();
    }
    if( iproblem == 3){
      SolveProblem_LinearSolid_Dynamic();
    }
    if( iproblem == 5 || iproblem == 8 ){
      SolveProblem_Stokes_Dynamic();
    }
    if( iproblem == 6 || iproblem == 9 ){
      SolveProblem_NavierStokes_Dynamic();
    }
  }
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
  GLint viewport[4];
  ::glGetIntegerv(GL_VIEWPORT,viewport);
  const int win_w = viewport[2];
  const int win_h = viewport[3];
  const double mov_end_x = (2.0*x-win_w)/win_w;
  const double mov_end_y = (win_h-2.0*y)/win_h;
  mov_begin_x = mov_end_x;
  mov_begin_y = mov_end_y;
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  GLint viewport[4];
  ::glGetIntegerv(GL_VIEWPORT,viewport);
  const int win_w = viewport[2];
  const int win_h = viewport[3];
  mov_begin_x = (2.0*x-win_w)/win_w;
  mov_begin_y = (win_h-2.0*y)/win_h;
  press_button = button;
  if( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ){
  }
  else if( button == GLUT_LEFT_BUTTON && state == GLUT_UP){
    press_button = -1;
  }
}

void Solve(bool is_next)
{
  if( is_next ){ iproblem = (iproblem+1)%10; }
  const int np = aXY1.size()/2;
  if( iproblem == 0 ){
    is_animation = false;
    aVal.assign(np, 0.0);
    InitializeProblem_Scalar();
    SolveProblem_Poisson();
  }
  else if( iproblem == 1 ){
    is_animation = true;
    aVal.assign(np, 0.0);
    aVelo.assign(np, 0.0);
    InitializeProblem_Scalar();
    SolveProblem_Diffusion();
  }
  else if( iproblem == 2 ){
    is_animation = false;
    aVal.assign(np*2, 0.0);
    InitializeProblem_Solid();
    SolveProblem_LinearSolid_Static();
  }
  else if( iproblem == 3 ){
    is_animation = true;
    aVal.assign(np*2, 0.0);
    aVelo.assign(np*2, 0.0);
    aAcc.assign(np*2, 0.0);
    InitializeProblem_Solid();
    SolveProblem_LinearSolid_Dynamic();
  }
  else if( iproblem == 4 ){
    is_animation = false;
    aVal.assign(np*3, 0.0);
    InitializeProblem_Fluid();
    SolveProblem_Stokes_Static();
  }
  else if( iproblem == 5 ){
    is_animation = true;
    aVal.assign(np*3, 0.0);
    aVelo.assign(np*3, 0.0);
    InitializeProblem_Fluid();
    SolveProblem_Stokes_Dynamic();
  }
  else if( iproblem == 6 ){
    is_animation = true;
    aVal.assign(np*3, 0.0);
    aVelo.assign(np*3, 0.0);
    InitializeProblem_Fluid();
    SolveProblem_NavierStokes_Dynamic();
  }
  else if( iproblem == 7 ){
    is_animation = false;
    aVal.assign(np*3, 0.0);
    InitializeProblem_Fluid2();
    SolveProblem_Stokes_Static();
  }
  else if( iproblem == 8 ){
    is_animation = true;
    aVal.assign(np*3, 0.0);
    aVelo.assign(np*3, 0.0);
    InitializeProblem_Fluid2();
    SolveProblem_Stokes_Dynamic();
  }
  else if( iproblem == 9 ){
    is_animation = true;
    aVal.assign(np*3, 0.0);
    aVelo.assign(np*3, 0.0);
    InitializeProblem_Fluid2();
    SolveProblem_NavierStokes_Dynamic();
  }
}

void myGlutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
    case 'q':
    case 'Q':
    case '\033':  /* '\033' ÇÕ ESC ÇÃ ASCII ÉRÅ[Éh */
      exit(0);
      break;
    case 'a':
      is_animation = !is_animation; 
      break;
    case 'm':
      MakeMesh();
      Solve(false);
      break;
    case ' ':
      Solve(true);
      break;
    default:
      break;
  }
}

void myGlutSpecial(int key, int x, int y){
  switch(key){
    case GLUT_KEY_PAGE_UP:
      mag *= 1.0/0.9;
      break;
    case GLUT_KEY_PAGE_DOWN:
      mag *= 0.9;
      break;
  }
  ::myGlutResize(-1,-1);
  ::glutPostRedisplay();
}

int main(int argc,char* argv[])
{
  // Initailze GLUT
  ::glutInitWindowPosition(200,200);
  ::glutInitWindowSize(1000, 500);
  ::glutInit(&argc, argv);
  ::glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  ::glutCreateWindow("Cad View");
  
  // Set callback function
  ::glutMotionFunc(myGlutMotion);
  ::glutMouseFunc(myGlutMouse);
  ::glutDisplayFunc(myGlutDisplay);
  ::glutReshapeFunc(myGlutResize);
  ::glutKeyboardFunc(myGlutKeyboard);
  ::glutSpecialFunc(myGlutSpecial);
  ::glutIdleFunc(myGlutIdle);
  
  MakeMesh();
  Solve(false);
  //  InitializeProblem_Scalar();
  //  SolveProblem_Poisson();
  //InitializeProblem_LinearSolid_Static();
  //SolveProblem_LinearSolid_Static();
  //  InitializeProblem_LinearSolid_Dynamic();
  //  SolveProblem_LinearSolid_Dynamic();
  iproblem = 0;
  
  // Enter main loop
  ::glutMainLoop();
  return 0;
}
