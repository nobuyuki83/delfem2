/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#include "delfem2/mat3.h"
#include "delfem2/mat4.h"
#include "delfem2/quat.h"
#include "delfem2/mshtopo.h"
#include "delfem2/primitive.h"
#include "delfem2/vecxitrsol.h"

// ----------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// ------------------------------------

void SetPositionAtFixedBoundary(
    std::vector<double>& aRhs,
    unsigned int iframe,
    const std::vector<double>& aXYZ0,
    const std::vector<int>& aBCFlag)
{
  double A[16];
  {
    dfm2::Mat4_Identity(A);
    const double trans0[3] = {0, -0.8, 0};
    dfm2::Translate_Mat4Affine(A,
                            trans0);
    const double axis0[3] = {0, +2.0*sin(0.03*iframe), 1.0*sin(0.07*iframe)};
    dfm2::Rotate_Mat4AffineRodriguez(A,
                                   axis0);
    const double trans1[3] = {0.2*sin(0.03*iframe), +0.5+0.1*cos(0.05*iframe), 0};
    dfm2::Translate_Mat4Affine(A,
                            trans1);
  }
  const unsigned int np = aRhs.size()/3;
  for(unsigned int ip=0;ip<np;++ip){
    if( aBCFlag[ip*3+0] == 0 ){ continue; }
    if( aBCFlag[ip*3+0] == 1 ){
      aRhs[ip*3+0] = aXYZ0[ip*3+0];
      aRhs[ip*3+1] = aXYZ0[ip*3+1];
      aRhs[ip*3+2] = aXYZ0[ip*3+2];
    }
    if( aBCFlag[ip*3+0] == 2 ) {
      dfm2::Vec3_Mat4Vec3_AffineProjection(aRhs.data()+ip*3, A, aXYZ0.data()+ip*3);
    }
  }
}

class CDeformer_ARAPLinearDisponly {
public:
  CDeformer_ARAPLinearDisponly(
       const std::vector<unsigned int>& psup_ind0,
       const std::vector<unsigned int>& psup0,
       double weight_bc0,
       const std::vector<int>& aBCFlag0) :
  psup_ind(psup_ind0),
  psup(psup0),
  weight_bc(weight_bc0),
  aBCFlag(aBCFlag0)
  {
    const unsigned int np = aBCFlag.size()/3;
    const unsigned int ne = psup.size();
    // -----
    aMatEdge.resize(ne*9*2);
    assert(psup_ind.size()==np+1);
    for(unsigned int ip=0;ip<np;++ip){
      for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
        //          unsigned int jp = psup[ipsup];
        dfm2::Mat3_Identity(aMatEdge.data()+ipsup*18,   +1.0);
        dfm2::Mat3_Identity(aMatEdge.data()+ipsup*18+9, -1.0);
      }
    }
    // ---
    vec_tmp.resize(ne*3);
  }
  void JacobiTVecTmp(double*y ,
                     double alpha, double beta) const {
    const unsigned int np = aBCFlag.size()/3;
    for(int i=0;i<np*3;++i){ y[i] *= beta; }
    for(unsigned int ip=0;ip<np;++ip){
      for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
        unsigned int jp0 = psup[ipsup];
        dfm2::MatTVec3_ScaleAdd(y+ip*3,
                                aMatEdge.data()+ipsup*18,
                                vec_tmp.data()+ipsup*3,
                                alpha, 1.0);
        dfm2::MatTVec3_ScaleAdd(y+jp0*3,
                                aMatEdge.data()+ipsup*18+9,
                                vec_tmp.data()+ipsup*3,
                                alpha, 1.0);
      }
    }
  }
  void MatVec(double* y,
              double alpha, const double* vec,  double beta) const {
    const unsigned int np = aBCFlag.size()/3;
    std::fill(vec_tmp.begin(),vec_tmp.end(), 0.0);
    for(unsigned int ip=0;ip<np;++ip){
      for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
        unsigned int jp0 = psup[ipsup];
        dfm2::MatVec3_ScaleAdd(vec_tmp.data()+ipsup*3,
                               aMatEdge.data()+ipsup*18,
                               vec+ip*3,
                               1.0, 1.0);
        dfm2::MatVec3_ScaleAdd(vec_tmp.data()+ipsup*3,
                               aMatEdge.data()+ipsup*18+9,
                               vec+jp0*3,
                               1.0, 1.0);
      }
    }
    this->JacobiTVecTmp(y,
                        alpha, beta);
    // add diagonal for fixed boundary condition
    for(int i=0;i<aBCFlag.size();++i){
      if( aBCFlag[i] == 0 ){ continue; }
      y[i] += weight_bc*vec[i];
    }
  }
  void MakeLinearSystem(double* aRhs,
                        const double* aXYZ0,
                        const double* aXYZ1,
                        const double* aGoal) const
  {
    const unsigned int np = aBCFlag.size()/3;
    const unsigned int ne = psup.size();
    vec_tmp.assign(ne*3,0);
    for(unsigned int ip=0;ip<np;++ip){
      for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
        const unsigned int jp0 = psup[ipsup];
        const double d0[3] = { aXYZ0[jp0*3+0]-aXYZ0[ip*3+0], aXYZ0[jp0*3+1]-aXYZ0[ip*3+1], aXYZ0[jp0*3+2]-aXYZ0[ip*3+2] };
        const double d1[3] = { aXYZ1[jp0*3+0]-aXYZ1[ip*3+0], aXYZ1[jp0*3+1]-aXYZ1[ip*3+1], aXYZ1[jp0*3+2]-aXYZ1[ip*3+2] };
        vec_tmp[ipsup*3+0] += +(d0[0] - d1[0]);
        vec_tmp[ipsup*3+1] += +(d0[1] - d1[1]);
        vec_tmp[ipsup*3+2] += +(d0[2] - d1[2]);
      }
    }
    this->JacobiTVecTmp(aRhs,
                        -1.0, 0.0);
    // making RHS vector for fixed boundary condition
    for(int i=0;i<np*3;++i){
      if( aBCFlag[i] == 0 ){ continue; }
      aRhs[i] += (aGoal[i]-aXYZ1[i])*weight_bc;
    }
    
  }
public:
  const std::vector<unsigned int> psup_ind;
  const std::vector<unsigned int> psup;
  const double weight_bc;
  const std::vector<int> aBCFlag;
  // -------------
  std::vector<double> aMatEdge;
  mutable std::vector<double> vec_tmp;
};

class CDeformer_ARAP {
public:
  CDeformer_ARAP(const std::vector<unsigned int>& psup_ind0,
                 const std::vector<unsigned int>& psup0,
                 double weight_bc0,
                 const std::vector<int>& aBCFlag0) :
  psup_ind(psup_ind0),
  psup(psup0),
  weight_bc(weight_bc0),
  aBCFlag(aBCFlag0)
  {
    const unsigned int np = aBCFlag.size()/3;
    const unsigned int ne = psup.size();
    // -----
    aMatEdge.resize(ne*27);
    assert(psup_ind.size()==np+1);
    for(unsigned int ip=0;ip<np;++ip){
      for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
        //          unsigned int jp = psup[ipsup];
        dfm2::Mat3_Identity(aMatEdge.data()+ipsup*27+0, +1.0);
        dfm2::Mat3_Identity(aMatEdge.data()+ipsup*27+9, -1.0);
      }
    }
    // ---
    vec_tmp.resize(ne*3);
  }
  void JacobiTVecTmp(double*y ,
                     double alpha, double beta) const {
    const unsigned int np = aBCFlag.size()/3;
    for(int i=0;i<np*6;++i){ y[i] *= beta; }
    for(unsigned int ip=0;ip<np;++ip){
      for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
        unsigned int jp0 = psup[ipsup];
        dfm2::MatTVec3_ScaleAdd(y+ip*3,
                                aMatEdge.data()+ipsup*27+0,
                                vec_tmp.data()+ipsup*3,
                                alpha, 1.0);
        dfm2::MatTVec3_ScaleAdd(y+jp0*3,
                                aMatEdge.data()+ipsup*27+9,
                                vec_tmp.data()+ipsup*3,
                                alpha, 1.0);
        dfm2::MatTVec3_ScaleAdd(y+np*3+ip*3,
                                aMatEdge.data()+ipsup*27+18,
                                vec_tmp.data()+ipsup*3,
                                alpha, 1.0);
      }
    }
  }
  void MatVec(double* y,
              double alpha, const double* vec,  double beta) const {
    const unsigned int np = aBCFlag.size()/3;
    std::fill(vec_tmp.begin(),vec_tmp.end(), 0.0);
    for(unsigned int ip=0;ip<np;++ip){
      for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
        unsigned int jp0 = psup[ipsup];
        dfm2::MatVec3_ScaleAdd(vec_tmp.data()+ipsup*3,
                               aMatEdge.data()+ipsup*27+0,
                               vec+ip*3,
                               1.0, 1.0);
        dfm2::MatVec3_ScaleAdd(vec_tmp.data()+ipsup*3,
                               aMatEdge.data()+ipsup*27+9,
                               vec+jp0*3,
                               1.0, 1.0);
        dfm2::MatVec3_ScaleAdd(vec_tmp.data()+ipsup*3,
                               aMatEdge.data()+ipsup*27+18,
                               vec+(np+ip)*3,
                               1.0, 1.0);
      }
    }
    this->JacobiTVecTmp(y,
                        alpha, beta);
    // add diagonal for fixed boundary condition
    for(int i=0;i<aBCFlag.size();++i){
      if( aBCFlag[i] == 0 ){ continue; }
      y[i] += weight_bc*vec[i];
//      y[np*3+i] += weight_bc*vec[np*3+i];
    }
  }
  void MakeLinearSystem(double* aRhs,
                        const double* aXYZ0,
                        const double* aXYZ1,
                        const double* aGoal,
                        const double* aQuat)
  {
    const unsigned int np = aBCFlag.size()/3;
    const unsigned int ne = psup.size();
    vec_tmp.assign(ne*3,0);
    for(unsigned int ip=0;ip<np;++ip){
      for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
        const unsigned int jp0 = psup[ipsup];
        const double* q0 = aQuat+ip*4;
        const double d0[3] = { aXYZ0[jp0*3+0]-aXYZ0[ip*3+0], aXYZ0[jp0*3+1 ]-aXYZ0[ip*3+1], aXYZ0[jp0*3+2]-aXYZ0[ip*3+2] };
        const double d1[3] = { aXYZ1[jp0*3+0]-aXYZ1[ip*3+0], aXYZ1[jp0*3+1]-aXYZ1[ip*3+1], aXYZ1[jp0*3+2]-aXYZ1[ip*3+2] };
        double Rd0[3]; dfm2::QuatVec(Rd0, q0,d0);
        vec_tmp[ipsup*3+0] += +(Rd0[0] - d1[0]);
        vec_tmp[ipsup*3+1] += +(Rd0[1] - d1[1]);
        vec_tmp[ipsup*3+2] += +(Rd0[2] - d1[2]);
        dfm2::Mat3_Spin_ScaleAdd(
            aMatEdge.data()+ipsup*27+18,
            Rd0,
            -1.0, 0.0);
      }
    }
    this->JacobiTVecTmp(aRhs,
                        -1.0, 0.0);
    // making RHS vector for fixed boundary condition
    for(int i=0;i<np*3;++i){
      if( aBCFlag[i] == 0 ){ continue; }
      aRhs[i] += (aGoal[i]-aXYZ1[i])*weight_bc;
      //aRhs[i+np*3] = 0.0;
    }
  }
  void MakePreconditionerJacobi(){
    const unsigned int np = aBCFlag.size()/3;
    aDiaInv.assign(np*2*9, 0.0);
    for(unsigned int ip=0;ip<np;++ip){
      for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup) {
        const unsigned int jp0 = psup[ipsup];
        dfm2::MatTMat3_ScaleAdd(
            aDiaInv.data()+ip*9,
            aMatEdge.data()+ipsup*27+0,
            aMatEdge.data()+ipsup*27+0,
            1.0,1.0);
        dfm2::MatTMat3_ScaleAdd(
            aDiaInv.data()+jp0*9,
            aMatEdge.data()+ipsup*27+9,
            aMatEdge.data()+ipsup*27+9,
            1.0,1.0);
        dfm2::MatTMat3_ScaleAdd(
            aDiaInv.data()+(np+ip)*9,
            aMatEdge.data()+ipsup*27+18,
            aMatEdge.data()+ipsup*27+18,
            1.0,1.0);
      }
    }
    for(unsigned int ip=0;ip<np;++ip){
      for(int idim=0;idim<3;++idim) {
        if (aBCFlag[ip*3+idim] == 0) { continue; }
        aDiaInv[ip*9+idim*3+idim] += weight_bc;
      }
    }
    for(unsigned int ip=0;ip<np*2;++ip){
      dfm2::Inverse_Mat3(aDiaInv.data()+ip*9);
    }
  }
  void Solve(double* v) const {
    const unsigned int np = aBCFlag.size()/3;
    for(int ip=0;ip<np*2;++ip){
      double tmp[3];
      dfm2::MatVec3(tmp, aDiaInv.data()+ip*9, v+ip*3);
      v[ip*3+0] = tmp[0];
      v[ip*3+1] = tmp[1];
      v[ip*3+2] = tmp[2];
    }
  }
public:
  const std::vector<unsigned int> psup_ind;
  const std::vector<unsigned int> psup;
  const double weight_bc;
  const std::vector<int> aBCFlag;
  // -------------
  std::vector<double> aMatEdge;
  std::vector<double> aDiaInv; // for jacobi preconditining
  mutable std::vector<double> vec_tmp;
};

// ------------------------------------------

void SetProblem(std::vector<double>& aXYZ0,
                std::vector<unsigned int>& aTri,
                std::vector<int>& aBCFlag)

{
  dfm2::MeshTri3D_CylinderClosed(aXYZ0, aTri,
                                 0.2, 1.6,
                                 16, 16);
  {
    const unsigned int np = aXYZ0.size() / 3;
    aBCFlag.assign(np * 3, 0);
    for(unsigned int ip=0;ip<np;++ip) {
      double y0 = aXYZ0[ip*3+1];
      if( y0 < -0.65 ){
        aBCFlag[ip*3+0] = 1;
        aBCFlag[ip*3+1] = 1;
        aBCFlag[ip*3+2] = 1;
      }
      if( y0 > +0.65 ){
        aBCFlag[ip*3+0] = 2;
        aBCFlag[ip*3+1] = 2;
        aBCFlag[ip*3+2] = 2;
      }
    }
  }
}

// -------------------------------------------

void myGlutDisplay_Mesh
(const std::vector<double>& aXYZ0,
const std::vector<double>& aXYZ1,
const std::vector<unsigned int>& aTri)
{
  ::glLineWidth(1);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,0,0);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1,aTri);
  ::glColor3d(0.8,0.8,0.8);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ0.data(),aXYZ0.size()/3,
                                   aTri.data(),aTri.size()/3);
  ::glColor3d(0,0,0);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ1.data(),aXYZ1.size()/3,
                                   aTri.data(),aTri.size()/3);
  
}

void Draw_BCFlag(const std::vector<double>& aXYZ1,
                 const std::vector<int>& aBCFlag)
{ // draw bc as a point
  ::glPointSize(10);
  ::glBegin(GL_POINTS);
  for(unsigned int ip=0;ip<aXYZ1.size()/3;++ip){
    if( aBCFlag[ip*3+0] == 0 ){ continue; }
    if( aBCFlag[ip*3+0] == 1 ){ ::glColor3d(0,0,1); }
    if( aBCFlag[ip*3+0] == 2 ){ ::glColor3d(0,1,0); }
    ::glVertex3dv(aXYZ1.data()+ip*3);
  }
  ::glEnd();
}



// --------------------------------------------------

int main(int argc,char* argv[])
{
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  delfem2::opengl::setSomeLighting();
  
  std::vector<unsigned int> aTri;
  std::vector<double> aXYZ0, aXYZ1;
  std::vector<int> aBCFlag;
  std::vector<unsigned int> psup_ind, psup;
  {
    SetProblem(aXYZ0, aTri, aBCFlag);
    dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                               aTri.data(), aTri.size()/3, 3,
                               (int)aXYZ0.size()/3);
    dfm2::JArray_Sort(psup_ind, psup);
    aXYZ1 = aXYZ0;
  }

  
  while (true){
    const double weight_bc = 100.0;
    int iframe = 0;
    {
      CDeformer_ARAPLinearDisponly def0(psup_ind, psup, weight_bc, aBCFlag);
      const unsigned int np = aXYZ0.size()/3;
      glfwSetWindowTitle(viewer.window, "Linear Disponly");
      for(;iframe<50;++iframe)
      {
        for(int i=0;i<np*3;++i){ // adding noise for debuggng purpose
          aXYZ1[i] += 0.02*(double)rand()/(RAND_MAX+1.0)-0.01;
        }
        std::vector<double> aGoal(np*3);
        SetPositionAtFixedBoundary(aGoal,
                                       iframe,aXYZ0,aBCFlag);
        std::vector<double> aRhs(np*3,0.0);
        def0.MakeLinearSystem(aRhs.data(),
                              aXYZ0.data(), aXYZ1.data(), aGoal.data());
        std::vector<double> aUpd(np*3,0.0);
        std::vector<double> aRes = dfm2::Solve_CG(aRhs.data(), aUpd.data(),
                                                  np*3, 1.0e-4, 300, def0);
        std::cout << "iframe: " << iframe << "   nitr:" << aRes.size() << std::endl;
        for(int i=0;i<np*3;++i){ aXYZ1[i] += aUpd[i]; }
        // ------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(aXYZ0,aXYZ1,aTri);
        Draw_BCFlag(aXYZ1,aBCFlag);
        viewer.DrawEnd_oldGL();
        if( glfwWindowShouldClose(viewer.window) ){ goto CLOSE; }
      }
    } // end linear disponly
    // -------------------------------------------------------
    { // begin lienar disprot without preconditioner
      glfwSetWindowTitle(viewer.window, "Linear Disprot without Prec");
      unsigned int np = aXYZ0.size()/3;
      CDeformer_ARAP def1(psup_ind, psup, weight_bc, aBCFlag);
      std::vector<double> aQuat(np*4); // array of quaternion
      for(;iframe<100;++iframe){
        for(int ip=0;ip<np;++ip){ dfm2::Quat_Identity(aQuat.data()+ip*4); }
        for(int i=0;i<np*3;++i){ // adding noise for debuggng purpose
          aXYZ1[i] += 0.02*(double)rand()/(RAND_MAX+1.0)-0.01;
        }
        std::vector<double> aGoal(np*3);
        SetPositionAtFixedBoundary(aGoal,
                                       iframe,aXYZ0,aBCFlag);
        std::vector<double> aRhs(np*6,0.0);
        def1.MakeLinearSystem(aRhs.data(),
                              aXYZ0.data(), aXYZ1.data(), aGoal.data(),
                              aQuat.data());
        std::vector<double> aUpd(np*6,0.0);
        std::vector<double> aRes = dfm2::Solve_CG(aRhs.data(), aUpd.data(),
                                                  np*6, 1.0e-4, 400, def1);
        std::cout << "iframe:" << iframe << "   itr:" << aRes.size() << std::endl;
        for(int ip=0;ip<np;++ip){
          dfm2::Add3(aXYZ1.data()+ip*3, aUpd.data()+ip*3);
          double q0[4]; dfm2::Quat_CartesianAngle(q0, aUpd.data()+np*3+ip*3);
          double q1[4]; dfm2::QuatQuat(q1, q0, aQuat.data()+ip*4);
          dfm2::Copy_Quat(aQuat.data()+ip*4, q1);
        }
        // ------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(aXYZ0,aXYZ1, aTri);
        Draw_BCFlag(aXYZ1,aBCFlag);
        dfm2::opengl::Draw_QuaternionsCoordinateAxes(aXYZ1,aQuat,0.04);
        viewer.DrawEnd_oldGL();
        if( glfwWindowShouldClose(viewer.window) ){ goto CLOSE; }
      } // end of frame loop
    } // end linear disprot without preconditioner
    // -------------------------------
    { // begin lienar disprot with preconditioner
      glfwSetWindowTitle(viewer.window, "Linear Disprot with Prec");
      const unsigned int np = aXYZ0.size()/3;
      CDeformer_ARAP def1(psup_ind, psup, weight_bc, aBCFlag);
      std::vector<double> aQuat(np*4);
      for(;iframe<200;++iframe){
        for(int ip=0;ip<np;++ip){ dfm2::Quat_Identity(aQuat.data()+ip*4); }
        for(int i=0;i<np*3;++i){ // adding noise for debuggng purpose
          aXYZ1[i] += 0.02*(double)rand()/(RAND_MAX+1.0)-0.01;
        }
        std::vector<double> aGoal(np*3);
        SetPositionAtFixedBoundary(aGoal,
                                       iframe,aXYZ0,aBCFlag);
        std::vector<double> aRhs(np*6,0.0);
        def1.MakeLinearSystem(aRhs.data(),
                              aXYZ0.data(), aXYZ1.data(), aGoal.data(),
                              aQuat.data());
        def1.MakePreconditionerJacobi();
        std::vector<double> aUpd(np*6,0.0);
        std::vector<double> aRes = dfm2::Solve_PCG(aRhs.data(), aUpd.data(),
                                                  np*6, 1.0e-4, 400, def1, def1);
        std::cout << "iframe:" << iframe << "   itr:" << aRes.size() << std::endl;
        for(int ip=0;ip<np;++ip){
          dfm2::Add3(aXYZ1.data()+ip*3, aUpd.data()+ip*3);
          double q0[4]; dfm2::Quat_CartesianAngle(q0, aUpd.data()+np*3+ip*3);
          double q1[4]; dfm2::QuatQuat(q1, q0,aQuat.data()+ip*4);
          dfm2::Copy_Quat(aQuat.data()+ip*4, q1);
        }
        // ------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(aXYZ0,aXYZ1, aTri);
        Draw_BCFlag(aXYZ1,aBCFlag);
        dfm2::opengl::Draw_QuaternionsCoordinateAxes(aXYZ1,aQuat,0.04);
        viewer.DrawEnd_oldGL();
        if( glfwWindowShouldClose(viewer.window) ){ goto CLOSE; }
      } // end of frame loop
    } // end linear disprot with preconditioner
    // -------------------------------
    { // begin nonlienar disprot with preconditioner
      glfwSetWindowTitle(viewer.window, "NonLinear Disprot with Prec");
      const unsigned int np = aXYZ0.size()/3;
      CDeformer_ARAP def1(psup_ind, psup, weight_bc, aBCFlag);
      std::vector<double> aQuat(np*4);
      for(int ip=0;ip<np;++ip){ dfm2::Quat_Identity(aQuat.data()+ip*4); }
      for(;iframe<400;++iframe){
        /*
        for(int i=0;i<np*3;++i){ // adding noise for debuggng purpose
          aXYZ1[i] += 0.02*(double)rand()/(RAND_MAX+1.0)-0.01;
        }
         */
        std::vector<double> aGoal(np*3);
        SetPositionAtFixedBoundary(aGoal,
                                       iframe,aXYZ0,aBCFlag);
        std::vector<double> aRhs(np * 6, 0.0);
        def1.MakeLinearSystem(aRhs.data(),
                              aXYZ0.data(), aXYZ1.data(), aGoal.data(),
                              aQuat.data());
        def1.MakePreconditionerJacobi();
        std::vector<double> aUpd(np * 6, 0.0);
        std::vector<double> aRes = dfm2::Solve_PCG(aRhs.data(), aUpd.data(),
                                                   np * 6, 1.0e-5, 400, def1, def1);
        std::cout << "iframe:" << iframe << "   itr:" << aRes.size() << std::endl;
        for (int ip = 0; ip < np; ++ip) {
          dfm2::Add3(aXYZ1.data() + ip * 3, aUpd.data() + ip * 3);
          double q0[4];
          dfm2::Quat_CartesianAngle(q0, aUpd.data() + np * 3 + ip * 3);
          double q1[4];
          dfm2::QuatQuat(q1, q0, aQuat.data() + ip * 4);
          dfm2::Copy_Quat(aQuat.data() + ip * 4, q1);
        }
        // ------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(aXYZ0,aXYZ1, aTri);
        Draw_BCFlag(aXYZ1,aBCFlag);
        dfm2::opengl::Draw_QuaternionsCoordinateAxes(aXYZ1,aQuat,0.04);
        viewer.DrawEnd_oldGL();
        if( glfwWindowShouldClose(viewer.window) ){ goto CLOSE; }
      } // end of frame loop
    } // end linear disprot with preconditioner
  }
  
CLOSE:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


