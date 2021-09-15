/*
 * Copyright (c) 2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/femsolidlinear.h"

namespace delfem2 {
namespace femsolidlinear{

template <int nno>
void AddEmat_LinearSolid2(
    double emat[nno][nno][2][2],
    double lambda,
    double myu,
    const double dldx[nno][2],
    double w)
{
  for(unsigned int in=0;in<nno;++in){
    for(unsigned int jn=0;jn<nno;++jn){
      emat[in][jn][0][0] += w*(lambda+myu)*dldx[in][0]*dldx[jn][0];
      emat[in][jn][0][1] += w*(lambda*dldx[in][0]*dldx[jn][1]+myu*dldx[jn][0]*dldx[in][1]);
      emat[in][jn][1][0] += w*(lambda*dldx[in][1]*dldx[jn][0]+myu*dldx[jn][1]*dldx[in][0]);
      emat[in][jn][1][1] += w*(lambda+myu)*dldx[in][1]*dldx[jn][1];
      const double dtmp1 = w*myu*(dldx[in][1]*dldx[jn][1]+dldx[in][0]*dldx[jn][0]);
      emat[in][jn][0][0] += dtmp1;
      emat[in][jn][1][1] += dtmp1;
    }
  }
}

template <int nno>
void SetEmat_LinearSolid2(
    double emat[nno][nno][2][2],
    double lambda,
    double myu,
    const double dldx[nno][2],
    double w)
{
  for(unsigned int in=0;in<nno;++in){
    for(unsigned int jn=0;jn<nno;++jn){
      emat[in][jn][0][0] = w*(lambda+myu)*dldx[in][0]*dldx[jn][0];
      emat[in][jn][0][1] = w*(lambda*dldx[in][0]*dldx[jn][1]+myu*dldx[jn][0]*dldx[in][1]);
      emat[in][jn][1][0] = w*(lambda*dldx[in][1]*dldx[jn][0]+myu*dldx[jn][1]*dldx[in][0]);
      emat[in][jn][1][1] = w*(lambda+myu)*dldx[in][1]*dldx[jn][1];
      const double dtmp1 = w*myu*(dldx[in][1]*dldx[jn][1]+dldx[in][0]*dldx[jn][0]);
      emat[in][jn][0][0] += dtmp1;
      emat[in][jn][1][1] += dtmp1;
    }
  }
}

template <int nno>
void AddEmat_LinearSolid3(
    double emat[nno][nno][3][3],
    double lambda,
    double myu,
    const double dNdx[nno][3],
    double w)
{
  for (unsigned int ino = 0; ino<nno; ino++){
    for (unsigned int jno = 0; jno<nno; jno++){
      emat[ino][jno][0][0] += w*(lambda*dNdx[ino][0]*dNdx[jno][0]+myu*dNdx[jno][0]*dNdx[ino][0]);
      emat[ino][jno][0][1] += w*(lambda*dNdx[ino][0]*dNdx[jno][1]+myu*dNdx[jno][0]*dNdx[ino][1]);
      emat[ino][jno][0][2] += w*(lambda*dNdx[ino][0]*dNdx[jno][2]+myu*dNdx[jno][0]*dNdx[ino][2]);
      emat[ino][jno][1][0] += w*(lambda*dNdx[ino][1]*dNdx[jno][0]+myu*dNdx[jno][1]*dNdx[ino][0]);
      emat[ino][jno][1][1] += w*(lambda*dNdx[ino][1]*dNdx[jno][1]+myu*dNdx[jno][1]*dNdx[ino][1]);
      emat[ino][jno][1][2] += w*(lambda*dNdx[ino][1]*dNdx[jno][2]+myu*dNdx[jno][1]*dNdx[ino][2]);
      emat[ino][jno][2][0] += w*(lambda*dNdx[ino][2]*dNdx[jno][0]+myu*dNdx[jno][2]*dNdx[ino][0]);
      emat[ino][jno][2][1] += w*(lambda*dNdx[ino][2]*dNdx[jno][1]+myu*dNdx[jno][2]*dNdx[ino][1]);
      emat[ino][jno][2][2] += w*(lambda*dNdx[ino][2]*dNdx[jno][2]+myu*dNdx[jno][2]*dNdx[ino][2]);
      const double dtmp1 = dNdx[ino][0]*dNdx[jno][0]+dNdx[ino][1]*dNdx[jno][1]+dNdx[ino][2]*dNdx[jno][2];
      emat[ino][jno][0][0] += w*myu*dtmp1;
      emat[ino][jno][1][1] += w*myu*dtmp1;
      emat[ino][jno][2][2] += w*myu*dtmp1;
    }
  }
}

}
}



// ----------------------------
// below: tri

DFM2_INLINE void delfem2::EMat_SolidStaticLinear_Tri2D(
    double eres[3][2],
    double emat[3][3][2][2],
    const double myu,
    const double lambda,
    const double rho,
    const double g_x,
    const double g_y,
    const double disp[3][2],
    const double coords[3][2])
{
  constexpr int nno = 3;
  constexpr int ndim = 2;

  const double area = femutil::TriArea2D(coords[0],coords[1],coords[2]);
  double dldx[nno][ndim];
  {
    double zero_order_term[nno];
    TriDlDx(dldx, zero_order_term, coords[0], coords[1], coords[2]);
  }
  femsolidlinear::SetEmat_LinearSolid2<3>(
      emat,
      lambda,myu,dldx,area);

  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = area*rho*g_x*0.33333333333333333;
    eres[ino][1] = area*rho*g_y*0.33333333333333333;
  }

  // ----------------------------------------------
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1];
      eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1];
    }
  }
}

DFM2_INLINE void delfem2::EMat_SolidDynamicLinear_Tri2D(
    double eres[3][2],
    double emat[3][3][2][2],
    const double myu,
    const double lambda,
    const double rho,
    const double g_x,
    const double g_y,
    const double dt_timestep,
    [[maybe_unused]] const double gamma_newmark,
    const double beta_newmark,
    const double disp[3][2],
    const double velo[3][2],
    const double acc[3][2],
    const double coords[3][2],
    bool is_initial)
{
  constexpr int nno = 3;
  constexpr int ndim = 2;

  const double area = femutil::TriArea2D(coords[0],coords[1],coords[2]);
  double dldx[nno][ndim];   // spatial derivative of linear shape function
  {
    double zero_order_term[nno];  // const term of shape function
    TriDlDx(dldx, zero_order_term,   coords[0], coords[1], coords[2]);
  }

  double eKmat[nno][nno][ndim][ndim];
  femsolidlinear::SetEmat_LinearSolid2<3>(
      eKmat,
      lambda,myu,dldx,area);

  double eMmat[nno][nno][ndim][ndim];
  EmatConsistentMassTri2<2>(
      eMmat,
      area*rho,false);

  // calc external force
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = area*rho*g_x*0.33333333333333333333333333;
    eres[ino][1] = area*rho*g_y*0.33333333333333333333333333;
  }


  { // calc coeff matrix for newmark-beta
    double dtmp1 = beta_newmark*dt_timestep*dt_timestep;
    for(int i=0;i<nno*nno*ndim*ndim;i++){
      (&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i]+dtmp1*(&eKmat[0][0][0][0])[i];
    }
  }

  // calc element redisual vector
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres[ino][0] -= eKmat[ino][jno][0][0]*disp[jno][0]+eKmat[ino][jno][0][1]*disp[jno][1];
      eres[ino][1] -= eKmat[ino][jno][1][0]*disp[jno][0]+eKmat[ino][jno][1][1]*disp[jno][1];
    }
    for(int jno=0;jno<nno;jno++){
      eres[ino][0] -= eMmat[ino][jno][0][0]*acc[jno][0]+eMmat[ino][jno][0][1]*acc[jno][1];
      eres[ino][1] -= eMmat[ino][jno][1][0]*acc[jno][0]+eMmat[ino][jno][1][1]*acc[jno][1];
    }
  }
  if( is_initial ){
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eres[ino][0] -= dt_timestep*(eKmat[ino][jno][0][0]*velo[jno][0]+eKmat[ino][jno][0][1]*velo[jno][1]);
        eres[ino][1] -= dt_timestep*(eKmat[ino][jno][1][0]*velo[jno][0]+eKmat[ino][jno][1][1]*velo[jno][1]);
      }
      for(int jno=0;jno<nno;jno++){
        eres[ino][0] -= 0.5*dt_timestep*dt_timestep*(eKmat[ino][jno][0][0]*acc[jno][0]+eKmat[ino][jno][0][1]*acc[jno][1]);
        eres[ino][1] -= 0.5*dt_timestep*dt_timestep*(eKmat[ino][jno][1][0]*acc[jno][0]+eKmat[ino][jno][1][1]*acc[jno][1]);
      }
    }
  }
}

// above: tri
// -----------------------------------------
// below: quad

void delfem2::EMat_SolidLinear2_QuadOrth_GaussInt(
    double emat[4][4][2][2],
    double lx,
    double ly,
    double myu,
    double lambda,
    unsigned int ngauss)
{
  namespace lcl = delfem2::femutil;
  for(unsigned int i=0;i<16*4;++i){ (&emat[0][0][0][0])[i] = 0.0; }
  unsigned int nw = NIntLineGauss[ngauss];
  for(unsigned int iw=0;iw<nw;++iw){
    for(unsigned int jw=0;jw<nw;++jw){
      const double w = lx*ly*0.25*LineGauss<double>[ngauss][iw][1]*LineGauss<double>[ngauss][jw][1];
      const double x1 = (1-LineGauss<double>[ngauss][iw][0])*0.5;
      const double y1 = (1-LineGauss<double>[ngauss][jw][0])*0.5;
      const double x2 = 1 - x1;
      const double y2 = 1 - y1;
      // u = u1*(x1y1) + u2*(x2y1) + u3*(x2y2) + u4*(x1y2)
      // l1 = x1y1, l2=x2y1, l3=x2y2, l4=x1y2
      const double dldx[4][2] = {
        {-y1/lx, -x1/ly},
        {+y1/lx, -x2/ly},
        {+y2/lx, +x2/ly},
        {-y2/lx, +x1/ly} };
      femsolidlinear::AddEmat_LinearSolid2<4>(emat,lambda,myu,dldx,w);
    }
  }
}

DFM2_INLINE void delfem2::elemMatRes_LinearSolidGravity3_Static_Q1(
    const double myu,
    const double lambda,
    const double rho,
    const double g[3],
    const double coords[8][3],
    const double disp[8][3],
    //
    double emat[8][8][3][3],
    double eres[8][3])
{
  constexpr int iGauss = 2;
  constexpr int nInt = NIntLineGauss[iGauss];
  const double (*gauss)[2] = LineGauss<double>[iGauss];

  std::fill_n(&emat[0][0][0][0],8*8*3*3, 0.0);
  std::fill_n(&eres[0][0], 8*3, 0.0);

  double vol = 0.0;
  for(int ir1=0;ir1<nInt;ir1++){
    for(int ir2=0;ir2<nInt;ir2++){
      for(int ir3=0;ir3<nInt;ir3++){
        const double r1 = gauss[ir1][0];
        const double r2 = gauss[ir2][0];
        const double r3 = gauss[ir3][0];
        double dndx[8][3], an[8], detjac;
        ShapeFunc_Hex8(r1, r2, r3, coords, detjac, dndx, an);
        const double detwei = detjac * gauss[ir1][1] * gauss[ir2][1] * gauss[ir3][1];
        vol += detwei;
        femsolidlinear::AddEmat_LinearSolid3<8>(emat,lambda,myu,dndx,detwei);
        for(int ino=0;ino<8;ino++){
          eres[ino][0] += detwei*rho*g[0]*an[ino];
          eres[ino][1] += detwei*rho*g[1]*an[ino];
          eres[ino][2] += detwei*rho*g[2]*an[ino];
        }
      }
    }
  }
  AddEmatEvecScale3<8>(
      eres,
      emat,disp,-1.);
}


// above: quad
// -------------------------------
// below: tet

DFM2_INLINE void delfem2::stress_LinearSolid_TetP2(
    double stress[3][3],
    const double l0,
    const double l1,
    const double l2,
    const double l3,
    [[maybe_unused]] const double vol,
    const double lambda,
    const double myu,
    [[maybe_unused]] const double g_x,
    [[maybe_unused]] const double g_y,
    [[maybe_unused]] const double g_z,
    const double dldx[4][3],
    const double disp[10][3])
{
  /*
  double N[10] = {
    l0*(2*l0-1),
    l1*(2*l1-1),
    l2*(2*l2-1),
    l3*(2*l3-1),
    4*l0*l1,
    4*l1*l2,
    4*l0*l2,
    4*l0*l3,
    4*l1*l3,
    4*l2*l3,
  };
   */
  double dNdx[10][3];
  for (unsigned int i = 0; i<3; i++){
    dNdx[0][i] = (4*l0-1)*dldx[0][i];
    dNdx[1][i] = (4*l1-1)*dldx[1][i];
    dNdx[2][i] = (4*l2-1)*dldx[2][i];
    dNdx[3][i] = (4*l3-1)*dldx[3][i];
    dNdx[4][i] = 4*dldx[0][i]*l1+4*l0*dldx[1][i];
    dNdx[5][i] = 4*dldx[1][i]*l2+4*l1*dldx[2][i];
    dNdx[6][i] = 4*dldx[0][i]*l2+4*l0*dldx[2][i];
    dNdx[7][i] = 4*dldx[0][i]*l3+4*l0*dldx[3][i];
    dNdx[8][i] = 4*dldx[1][i]*l3+4*l1*dldx[3][i];
    dNdx[9][i] = 4*dldx[2][i]*l3+4*l2*dldx[3][i];
  }
  double dudx[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
  for (unsigned int ino = 0; ino<10; ino++){
    for (unsigned int i = 0; i<3; i++){
      for (unsigned int j = 0; j<3; j++){
        dudx[i][j] += disp[ino][i]*dNdx[ino][j];
      }
    }
  }
  double strain[3][3];
  for (unsigned int i = 0; i<3; i++){
    for (unsigned int j = 0; j<3; j++){
      strain[i][j] = 0.5*(dudx[i][j]+dudx[j][i]);
    }
  }
  {
    for (unsigned int i = 0; i<3; i++){
      for (unsigned int j = 0; j<3; j++){
        stress[i][j] = 2*myu*strain[i][j];
      }
      stress[i][i] += lambda*(strain[0][0]+strain[1][1]+strain[2][2]);
    }
  }
}


DFM2_INLINE void delfem2::matRes_LinearSolid_TetP2(
    double emat[10][10][3][3],
    double eres[10][3],
    const double vol, const double lambda, const double myu,
    const double g_x, const double g_y, const double g_z,
    const double rho,
    const double dldx[4][3],
    const double disp[10][3])
{
  constexpr unsigned int nOrder = 2;
  constexpr unsigned int nInt = NIntTetGauss[nOrder];
  const double (*gauss)[4] = TetGauss<double>[nOrder];
  std::fill_n(&emat[0][0][0][0], 10*10*3*3, 0.0);
  std::fill_n( &eres[0][0], 10*3, 0.0);
  for (unsigned int iint = 0; iint<nInt; iint++){
    double l0 = gauss[iint][0];
    double l1 = gauss[iint][1];
    double l2 = gauss[iint][2];
    double l3 = (1-l0-l1-l2);
    double w = gauss[iint][3];
    double N[10] = {
      l0*(2*l0-1),
      l1*(2*l1-1),
      l2*(2*l2-1),
      l3*(2*l3-1),
      4*l0*l1,
      4*l1*l2,
      4*l0*l2,
      4*l0*l3,
      4*l1*l3,
      4*l2*l3,
    };
    double dNdx[10][3];
    for (unsigned int i = 0; i<3; i++){
      dNdx[0][i] = (4*l0-1)*dldx[0][i];
      dNdx[1][i] = (4*l1-1)*dldx[1][i];
      dNdx[2][i] = (4*l2-1)*dldx[2][i];
      dNdx[3][i] = (4*l3-1)*dldx[3][i];
      dNdx[4][i] = 4*dldx[0][i]*l1+4*l0*dldx[1][i];
      dNdx[5][i] = 4*dldx[1][i]*l2+4*l1*dldx[2][i];
      dNdx[6][i] = 4*dldx[0][i]*l2+4*l0*dldx[2][i];
      dNdx[7][i] = 4*dldx[0][i]*l3+4*l0*dldx[3][i];
      dNdx[8][i] = 4*dldx[1][i]*l3+4*l1*dldx[3][i];
      dNdx[9][i] = 4*dldx[2][i]*l3+4*l2*dldx[3][i];
    }
    femsolidlinear::AddEmat_LinearSolid3<10>(emat,lambda,myu,dNdx,w*vol);
    for (unsigned int ino = 0; ino<10; ino++){
      eres[ino][0] += w*vol*rho*g_x*N[ino];
      eres[ino][1] += w*vol*rho*g_y*N[ino];
      eres[ino][2] += w*vol*rho*g_z*N[ino];
    }
  } // iint
  AddEmatEvecScale3<10>(
      eres,
      emat,disp,-1.);
}


DFM2_INLINE void delfem2::EMat_SolidLinear_Static_Tet(
    double emat[4][4][3][3],
    double eres[4][3],
    const double myu,
    const double lambda,
    const double P[4][3],
    const double disp[4][3],
    bool is_add)
{
  const double vol = femutil::TetVolume3D(P[0], P[1], P[2], P[3]);
  double dldx[4][3];
  {
    double const_term[4];    
    TetDlDx(dldx, const_term, P[0], P[1], P[2], P[3]);
  }
  // ----------------------
  if( is_add ){
    femsolidlinear::AddEmat_LinearSolid3<4>(emat,lambda,myu,dldx,vol);
  }
  else{
   SetEmat_LinearSolid3<4>(emat,lambda,myu,dldx,vol);
  }
  if( !is_add ){
    for(int i=0;i<12;++i){ (&eres[0][0])[i] = 0.0; }
  }
  AddEmatEvecScale3<4>(
      eres,
      emat,disp,-1.);
}

DFM2_INLINE void delfem2::EMat_SolidLinear_NewmarkBeta_MeshTet3D(
    double eres[4][3],
    double emat[4][4][3][3],
    const double myu,
    const double lambda,
    const double rho,
    const double g_x,
    const double g_y,
    const double g_z,
    const double dt,
    [[maybe_unused]] const double gamma_newmark,
    const double beta_newmark,
    const double disp[4][3],
    const double velo[4][3],
    const double acc[4][3],
    const double P[4][3],
    bool is_initial_iter)
{
  constexpr int nno = 4;
  constexpr int ndim = 3;
  
  const double vol = femutil::TetVolume3D(P[0],P[1],P[2],P[3]);
  double dldx[nno][ndim];		// spatial derivative of linear shape function
  {
    double zero_order_term[nno];	// const term of shape function
    TetDlDx(dldx, zero_order_term,   P[0],P[1],P[2],P[3]);
  }
  
  double eKmat[nno][nno][ndim][ndim];
  SetEmat_LinearSolid3<4>(eKmat,lambda,myu,dldx,vol);
  
  double eMmat[nno][nno][ndim][ndim];
  SetEmatConsistentMassTet(
      eMmat,
      rho*vol);
  
  //
  {	// calc coeff matrix for newmark-beta
    double dtmp1 = beta_newmark*dt*dt;
    for(int i=0;i<nno*nno*ndim*ndim;i++){
      (&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i]+dtmp1*(&eKmat[0][0][0][0])[i];
    }
  }

  // calc element redisual vector
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = vol*rho*g_x*0.25;
    eres[ino][1] = vol*rho*g_y*0.25;
    eres[ino][2] = vol*rho*g_z*0.25;
  }
  AddEmatEvecScale3<4>(
      eres,
      eKmat,disp,-1.);
  AddEmatEvecScale3<4>(
      eres,
      eMmat,acc,-1.);
  if( is_initial_iter ){
    AddEmatEvecScale3<4>(
        eres,
        eKmat,velo,-dt);
    AddEmatEvecScale3<4>(
        eres,
        eKmat,acc,-0.5*dt*dt);
  }
}

