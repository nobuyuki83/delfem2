#ifndef CLOTH_INTERNAL_H
#define CLOTH_INTERNAL_H

#include "delfem2/lsilu_mats.h"
#include "delfem2/lsmats.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/lsvecx.h"
#include "delfem2/femcloth.h"

namespace dfm2 = delfem2;

// -------------------------------------

/**
 * @brief compute total energy and its first and second derivatives
 */
void AddWdW_Cloth(
    double& W, // (out) energy
    std::vector<double>& dW, // (out) first derivative of energy
    //
    const std::vector<double>& aXYZ, // (in) deformed vertex positions
    const std::vector<double>& aXYZ0, // (in) initial vertex positions
    const std::vector<unsigned int>& aTri, // (in) triangle index
    const std::vector<unsigned int>& aQuad, // (in) index of 4 vertices required for bending
    double lambda, // (in) Lame's 1st parameter
    double myu,  // (in) Lame's 2nd parameter
    double stiff_bend // (in) bending stiffness
 )
{
  // marge element in-plane strain energy
  for(size_t itri=0;itri<aTri.size()/3;itri++){
    const unsigned int aIP[3] = { aTri[itri*3+0], aTri[itri*3+1], aTri[itri*3+2] };
    double C[3][3]; double c[3][3];
    for(int ino=0;ino<3;ino++){
      const int ip = aIP[ino];
      for(int i=0;i<3;i++){ C[ino][i] = aXYZ0[ip*3+i]; }
      for(int i=0;i<3;i++){ c[ino][i] = aXYZ [ip*3+i]; }
    }
    double e, de[3][3], dde[3][3][3][3];
    dfm2::WdWddW_CST( e,de,dde, C,c, lambda,myu );
    W += e;  // marge energy
    // marge de
    for(int ino=0;ino<3;ino++){
      const int ip = aIP[ino];
      for(int i =0;i<3;i++){ dW[ip*3+i] += de[ino][i]; }
    }
  }
  // marge element bending energy
  for(size_t iq=0;iq<aQuad.size()/4;iq++){
    const unsigned int aIP[4] = { aQuad[iq*4+0], aQuad[iq*4+1], aQuad[iq*4+2], aQuad[iq*4+3] };
    double C[4][3]; double c[4][3];
    for(int ino=0;ino<4;ino++){
      const int ip = aIP[ino];
      for(int i=0;i<3;i++){ C[ino][i] = aXYZ0[ip*3+i]; }
      for(int i=0;i<3;i++){ c[ino][i] = aXYZ [ip*3+i]; }
    }
    double e, de[4][3], dde[4][4][3][3];
    dfm2::WdWddW_Bend( e,de,dde, C,c, stiff_bend );
    W += e;  // marge energy
    // marge de
    for(int ino=0;ino<4;ino++){
      const int ip = aIP[ino];
      for(int i =0;i<3;i++){ dW[ip*3+i] += de[ino][i]; }
    }
  }
}


/**
 * @brief compute total energy and its first and second derivatives
 * @param dW first derivative of energy
 * @param aXYZ deformed vertex positions
 * @param gravity gravitational accerelation in xyz directions stored in a native array
 * @param mass_point mass of a point
 */
double AddWdW_Gravity
(std::vector<double>& dW,
 const std::vector<double>& aXYZ,
 const double gravity[3],
 double mass_point)
{
  double W = 0;
  for(unsigned int ip=0;ip<aXYZ.size()/3;ip++){
    const double c[3] = {aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]};
    W -= mass_point*( c[0]*gravity[0] + c[1]*gravity[1] + c[2]*gravity[2] );
    dW[ip*3+0] -= mass_point*gravity[0];
    dW[ip*3+1] -= mass_point*gravity[1];
    dW[ip*3+2] -= mass_point*gravity[2];
  }
  return W;
}

void StepTime_InternalDynamics(
    std::vector<double>& aXYZ, // (in,out) deformed vertex positions
    std::vector<double>& aUVW, // (in,out) deformed vertex velocity
    delfem2::CMatrixSparse<double>& mat_A,
    //
    const std::vector<double>& aXYZ0,// (in) initial vertex positions
    const std::vector<int>& aBCFlag, // (in) boundary condition flag (0:free 1:fixed)
    const std::vector<unsigned int>& aTri, // (in) triangle index
    const std::vector<unsigned int>& aQuad, // (in) index of 4 vertices required for bending
    const double dt, // (in) size of time step
    double lambda, // (in) Lame's 1st parameter
    double myu, // (in) Lame's 2nd parameter
    double stiff_bend, // (in) bending stiffness
    const double gravity[3], // (in) gravitatinal accereration
    double mass_point, // (in) mass for a point
    double stiff_contact,
    double contact_clearance,
    const dfm2::CInput_Contact& input_contact)
{
  const unsigned int np = static_cast<unsigned int>(aXYZ.size()/3); // number of point
  const unsigned int nDof = np*3; // degree of freedom
  // compute total energy and its first and second derivatives
  double W = 0;
  std::vector<double> vec_b(nDof,0);
	mat_A.setZero();
  std::vector<int> tmp_buffer(np,-1);
  W += delfem2::MergeLinSys_Cloth(
      mat_A,vec_b.data(),
      lambda,myu,stiff_bend,
      aXYZ0.data(), static_cast<unsigned int>(aXYZ0.size()/3), 3,
      aTri.data(),  static_cast<unsigned int>(aTri.size()/3),
      aQuad.data(), static_cast<unsigned int>(aQuad.size()/4),
      aXYZ.data());
  W += delfem2::MergeLinSys_Contact(
      mat_A, vec_b.data(),
      stiff_contact,contact_clearance,
      input_contact,
      aXYZ.data(), static_cast<unsigned int>(aXYZ.size()/3));
  W += AddWdW_Gravity(
      vec_b,
      aXYZ,
      gravity,
      mass_point);
  //std::cout << "energy : " << W << std::endl;
  // compute coefficient matrix and left-hand-side vector
  // Back-ward Eular time integration
  for(unsigned int i=0;i<nDof;i++){
    vec_b[i] = -vec_b[i] + mass_point*aUVW[i]/dt;
  }
  for(unsigned int ip=0;ip<np;ip++){
    mat_A.val_dia_[ip*9+0*3+0] += mass_point / (dt*dt);
    mat_A.val_dia_[ip*9+1*3+1] += mass_point / (dt*dt);
    mat_A.val_dia_[ip*9+2*3+2] += mass_point / (dt*dt);
  }
  mat_A.SetFixedBC(aBCFlag.data());
  for(unsigned int ip=0;ip<np;ip++){
    if( aBCFlag[ip] == 0 ) continue;
    vec_b[ip*3+0] = 0;
    vec_b[ip*3+1] = 0;
    vec_b[ip*3+2] = 0;
  }
  // solve linear system
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  {
    std::size_t n = vec_b.size();
    vec_x.resize(n);
    std::vector<double> tmp0(n), tmp1(n);
    auto vb = delfem2::CVecXd(vec_b);
    auto vx = delfem2::CVecXd(vec_x);
    auto vs = delfem2::CVecXd(tmp0);
    auto vt = delfem2::CVecXd(tmp1);
    Solve_CG(
        vx,vb,vs,vt,
        conv_ratio, iteration, mat_A);
  }
  std::cout << "  conv_ratio:" << conv_ratio << "  iteration:" << iteration << std::endl;
  // update position
  for(unsigned int i=0;i<nDof;i++){ aXYZ[i] += vec_x[i]; }
  // update velocity
  for(unsigned int i=0;i<nDof;i++){ aUVW[i] = vec_x[i]/dt; }
}

void StepTime_InternalDynamicsILU(
    std::vector<double>& aXYZ, // (in,out) deformed vertex positions
    std::vector<double>& aUVW, // (in,out) deformed vertex velocity
    delfem2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    //
    const std::vector<double>& aXYZ0,// (in) initial vertex positions
    const std::vector<int>& aBCFlag, // (in) boundary condition flag (0:free 1:fixed)
    const std::vector<unsigned int>& aTri, // (in) triangle index
    const std::vector<unsigned int>& aQuad, // (in) index of 4 vertices required for bending
    const double dt, // (in) size of time step
    double lambda, // (in) Lame's 1st parameter
    double myu, // (in) Lame's 2nd parameter
    double stiff_bend, // (in) bending stiffness
    const double gravity[3], // (in) gravitatinal accereration
    double mass_point, // (in) mass for a point
    double stiff_contact,
    double contact_clearance,
    const dfm2::CInput_Contact& input_contact)
{
  const unsigned int np = static_cast<unsigned int>(aXYZ.size()/3); 
  const unsigned int nDof = np*3;
  // compute total energy and its first and second derivatives
  double W = 0;
  std::vector<double> vec_b(nDof,0);
	mat_A.setZero();
  std::vector<int> tmp_buffer(np,-1);
  W += delfem2::MergeLinSys_Cloth(
      mat_A,vec_b.data(),
      lambda,myu,stiff_bend,
      aXYZ0.data(), static_cast<unsigned int>(aXYZ0.size()/3), 3,
      aTri.data(), static_cast<unsigned int>(aTri.size()/3),
      aQuad.data(), static_cast<unsigned int>(aQuad.size()/4),
      aXYZ.data());
  W += delfem2::MergeLinSys_Contact(
      mat_A,vec_b.data(),
      stiff_contact,contact_clearance,
      input_contact,
      aXYZ.data(), static_cast<unsigned int>(aXYZ.size()/3));
  W += AddWdW_Gravity(
      vec_b,
      aXYZ,
      gravity,mass_point);
//  std::cout << "energy : " << W << std::endl;
  // compute coefficient matrix and left-hand-side vector
  // Back-ward Eular time integration
  for(unsigned int i=0;i<nDof;i++){
    vec_b[i] = -vec_b[i] + mass_point*aUVW[i]/dt;
  }
  for(unsigned int ip=0;ip<np;ip++){
    mat_A.val_dia_[ip*9+0*3+0] += mass_point / (dt*dt);
    mat_A.val_dia_[ip*9+1*3+1] += mass_point / (dt*dt);
    mat_A.val_dia_[ip*9+2*3+2] += mass_point / (dt*dt);
  }
  mat_A.SetFixedBC(aBCFlag.data());
  for(unsigned int i=0;i<nDof;i++){
    if( aBCFlag[i] == 0 ) continue;
    vec_b[i] = 0;
  }
  ilu_A.CopyValue(mat_A);
  ilu_A.Decompose();
  // solve linear system
  double conv_ratio = 1.0e-4;
  int iteration = 100;
  std::vector<double> vec_x(vec_b.size());
  {
    const std::size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    auto vr = dfm2::CVecXd(vec_b);
    auto vu = dfm2::CVecXd(vec_x);
    auto vt = dfm2::CVecXd(tmp0);
    auto vs = dfm2::CVecXd(tmp1);
    Solve_PCG(
        vr, vu, vt, vs,
        conv_ratio, iteration, mat_A, ilu_A);
//    Solve_CG(
//        vr, vu, vt, vs,
//        conv_ratio, iteration, mat_A);
  }
//  std::cout << "  conv_ratio:" << conv_ratio << "  iteration:" << iteration << std::endl;
  // update position
  for(unsigned int i=0;i<nDof;i++){
    if( aBCFlag[i] != 0 ) continue;
    aXYZ[i] += vec_x[i];
  }
  // update velocity
  for(unsigned int i=0;i<nDof;i++){
    if( aBCFlag[i] != 0 ) continue;
    aUVW[i] = vec_x[i]/dt;
  }
}




/**
 * @brief Update intermediate velocity for cloth
 * @param aUVW  (in,out) deformed vertex velocity
 * @param aXYZ (in,out) deformed vertex position
 * @param aXYZ0 (in) initial vertex positions
 * @param aBCFlag (in) boundary condition flag (0:free else:fixed)
 * @param aTri (in) triangle index
 * @param aQuad (in) index of 4 vertices required for bending
 * @param dt (in) size of time step
 * @param lambda (in) Lame's 1st parameter
 * @param myu (in) Lame's 2nd parameter
 * @param stiff_bend (in) bending stiffness
 * @param gravity (in) gravitatinal accereration
 * @param mass_point  (in) mass for a point
 */
void UpdateIntermidiateVelocity
(std::vector<double>& aUVW,
 ////
 const std::vector<double>& aXYZ, // (in,out)
 const std::vector<double>& aXYZ0,// (in) initial vertex positions
 const std::vector<int>& aBCFlag, // (in) boundary condition flag (0:free else:fixed)
 const std::vector<unsigned int>& aTri, // (in) triangle index
 const std::vector<unsigned int>& aQuad, // (in) index of 4 vertices required for bending
 const double dt, // (in) size of time step
 double lambda, // (in) Lame's 1st parameter
 double myu, // (in) Lame's 2nd parameter
 double stiff_bend, // (in) bending stiffness
 const double gravity[3], // (in) gravitatinal accereration
 double mass_point // (in) mass for a point
 )
{
  const unsigned int np = static_cast<unsigned int>(aXYZ.size()/3); // number of point
  const unsigned int nDof = np*3; // degree of freedom
  
  // compute total energy and its first and second derivatives
  double W = 0;
  std::vector<double> dW(nDof,0);
  AddWdW_Cloth(W,dW,
               aXYZ,aXYZ0,
               aTri,aQuad,
               lambda,myu,stiff_bend);
  W += AddWdW_Gravity(dW,
                      aXYZ,
                      gravity,mass_point);
  for(unsigned int ip=0;ip<np;ip++){
    if( aBCFlag[ip] == 0 ) continue;
    aUVW[ip*3+0] = 0;
    aUVW[ip*3+1] = 0;
    aUVW[ip*3+2] = 0;
  }
  for(unsigned int ip=0;ip<np;ip++){
    aUVW[ip*3+0] += -(0.5*dt/mass_point)*dW[ip*3+0];
    aUVW[ip*3+1] += -(0.5*dt/mass_point)*dW[ip*3+1];
    aUVW[ip*3+2] += -(0.5*dt/mass_point)*dW[ip*3+2];
  }
}


// Setting problem here
void SetClothShape_Square
(std::vector<double>& aXYZ0, // (out) undeformed vertex positions
 std::vector<int>& aBCFlag, // (out) boundary condition flag (0:free 1:fixed)
 std::vector<unsigned int>& aTri, // (out) index of triangles
 std::vector<unsigned int>& aQuad, // (out) index of 4 vertices required for bending
 // ---------------------
 int ndiv, // (in) number of division of the square cloth edge
 double cloth_size) // (in) size of square cloth
{
  // make vertex potision array
  const double elem_length = cloth_size/ndiv; // size of an element
  const int nxyz =(ndiv+1)*(ndiv+1); // number of points
  aXYZ0.reserve( nxyz*3 );
  for(int ix=0;ix<ndiv+1;ix++){
    for(int iy=0;iy<ndiv+1;iy++){
      aXYZ0.push_back( ix*elem_length );
      aXYZ0.push_back( iy*elem_length );
      aXYZ0.push_back( 0.0 );
    }
  }
  
  // make triangle index array
  const int ntri = ndiv*ndiv*2;
  aTri.reserve( ntri*3 );
  for(int ix=0;ix<ndiv;ix++){
    for(int iy=0;iy<ndiv;iy++){
      aTri.push_back(  ix   *(ndiv+1)+ iy    );
      aTri.push_back( (ix+1)*(ndiv+1)+ iy    );
      aTri.push_back(  ix   *(ndiv+1)+(iy+1) );
      ////
      aTri.push_back( (ix+1)*(ndiv+1)+(iy+1) );
      aTri.push_back(  ix   *(ndiv+1)+(iy+1) );
      aTri.push_back( (ix+1)*(ndiv+1)+ iy    );
    }
  }
  
  // make quad index array
  const int nquad = ndiv*ndiv + ndiv*(ndiv-1)*2;
  aQuad.reserve( nquad*4 );
  for(int ix=0;ix<ndiv;ix++){
    for(int iy=0;iy<ndiv;iy++){
      aQuad.push_back( (ix+0)*(ndiv+1)+(iy+0) );
      aQuad.push_back( (ix+1)*(ndiv+1)+(iy+1) );
      aQuad.push_back( (ix+1)*(ndiv+1)+(iy+0) );
      aQuad.push_back( (ix+0)*(ndiv+1)+(iy+1) );
    }
  }
  for(int ix=0;ix<ndiv;ix++){
    for(int iy=0;iy<ndiv-1;iy++){
      aQuad.push_back( (ix+1)*(ndiv+1)+(iy+0) );
      aQuad.push_back( (ix+0)*(ndiv+1)+(iy+2) );
      aQuad.push_back( (ix+1)*(ndiv+1)+(iy+1) );
      aQuad.push_back( (ix+0)*(ndiv+1)+(iy+1) );
    }
  }
  for(int ix=0;ix<ndiv-1;ix++){
    for(int iy=0;iy<ndiv;iy++){
      aQuad.push_back( (ix+0)*(ndiv+1)+(iy+1) );
      aQuad.push_back( (ix+2)*(ndiv+1)+(iy+0) );
      aQuad.push_back( (ix+1)*(ndiv+1)+(iy+0) );
      aQuad.push_back( (ix+1)*(ndiv+1)+(iy+1) );
    }
  }
  
  aBCFlag = std::vector<int>(nxyz*3,0);
  for(int iy=0;iy<ndiv+1;iy++){
    aBCFlag[iy*3+0] = 1;
    aBCFlag[iy*3+1] = 1;
    aBCFlag[iy*3+2] = 1;
  }
}


#endif
