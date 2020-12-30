/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_HEADER_ONLY
// Merge use explicitly use the template so for static library we need to include the template itself.
#  include "delfem2/lsmats.h"
#endif
#include "delfem2/femsolidlinear.h"

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
      const double w = lx*ly*0.25*LineGauss[ngauss][iw][1]*LineGauss[ngauss][jw][1];
      const double x1 = (1-LineGauss[ngauss][iw][0])*0.5;
      const double y1 = (1-LineGauss[ngauss][jw][0])*0.5;
      const double x2 = 1 - x1;
      const double y2 = 1 - y1;
      // u = u1*(x1y1) + u2*(x2y1) + u3*(x2y2) + u4*(x1y2)
      // l1 = x1y1, l2=x2y1, l3=x2y2, l4=x1y2
      const double dldx[4][2] = {
        {-y1/lx, -x1/ly},
        {+y1/lx, -x2/ly},
        {+y2/lx, +x2/ly},
        {-y2/lx, +x1/ly} };
      for(unsigned int in=0;in<4;++in){
        for(unsigned int jn=0;jn<4;++jn){
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
  }
}


DFM2_INLINE void delfem2::ddW_SolidLinear_Tet3D(
    double* eKmat,
    double lambda,
    double myu,
    double vol,
    double dldx[4][3],
    bool is_add,
    unsigned int nstride)
{
  if( !is_add ){
    for(unsigned int i=0;i<4*4*nstride*nstride;++i){ eKmat[i] = 0.0; }
  }
  for (int ino = 0; ino<4; ino++){
    for (int jno = 0; jno<4; jno++){
      double* pK = eKmat+(nstride*nstride)*(ino*4+jno);
      pK[0*nstride+0] += vol*(lambda*dldx[ino][0]*dldx[jno][0]+myu*dldx[jno][0]*dldx[ino][0]);
      pK[0*nstride+1] += vol*(lambda*dldx[ino][0]*dldx[jno][1]+myu*dldx[jno][0]*dldx[ino][1]);
      pK[0*nstride+2] += vol*(lambda*dldx[ino][0]*dldx[jno][2]+myu*dldx[jno][0]*dldx[ino][2]);
      pK[1*nstride+0] += vol*(lambda*dldx[ino][1]*dldx[jno][0]+myu*dldx[jno][1]*dldx[ino][0]);
      pK[1*nstride+1] += vol*(lambda*dldx[ino][1]*dldx[jno][1]+myu*dldx[jno][1]*dldx[ino][1]);
      pK[1*nstride+2] += vol*(lambda*dldx[ino][1]*dldx[jno][2]+myu*dldx[jno][1]*dldx[ino][2]);
      pK[2*nstride+0] += vol*(lambda*dldx[ino][2]*dldx[jno][0]+myu*dldx[jno][2]*dldx[ino][0]);
      pK[2*nstride+1] += vol*(lambda*dldx[ino][2]*dldx[jno][1]+myu*dldx[jno][2]*dldx[ino][1]);
      pK[2*nstride+2] += vol*(lambda*dldx[ino][2]*dldx[jno][2]+myu*dldx[jno][2]*dldx[ino][2]);
      const double dtmp1 = dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]+dldx[ino][2]*dldx[jno][2];
      pK[0*nstride+0] += vol*myu*dtmp1;
      pK[1*nstride+1] += vol*myu*dtmp1;
      pK[2*nstride+2] += vol*myu*dtmp1;
    }
  }
}


DFM2_INLINE void delfem2::stress_LinearSolid_TetP2(
    double stress[3][3],
    const double l0,
    const double l1,
    const double l2,
    const double l3,
    const double vol,
    const double lambda,
    const double myu,
    const double g_x,
    const double g_y,
    const double g_z,
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
  for (unsigned int i = 0; i<10*10*3*3; i++){ (&emat[0][0][0][0])[i] = 0; }
  for (unsigned int i = 0; i<10*3; i++){ (&eres[0][0])[i] = 0; }
  unsigned int nOrder = 2;
  unsigned int nInt = NIntTetGauss[nOrder];
  for (unsigned int iint = 0; iint<nInt; iint++){
    double l0 = TetGauss[nOrder][iint][0];
    double l1 = TetGauss[nOrder][iint][1];
    double l2 = TetGauss[nOrder][iint][2];
    double l3 = (1-l0-l1-l2);
    double w = TetGauss[nOrder][iint][3];
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
    /*
    {
    double tN[4] = {0,0,0,0};
    for(unsigned int ino=0;ino<10;ino++){
    tN[0] += N[ino];
    tN[1] += dNdx[ino][0];
    tN[2] += dNdx[ino][1];
    tN[3] += dNdx[ino][2];
    }
    std::cout << tN[0] << "   " << tN[1] << " " << tN[2] << " " << tN[3] << std::endl;
    }
    */
    for (unsigned int ino = 0; ino<10; ino++){
      for (unsigned int jno = 0; jno<10; jno++){
        emat[ino][jno][0][0] += w*vol*(lambda*dNdx[ino][0]*dNdx[jno][0]+myu*dNdx[jno][0]*dNdx[ino][0]);
        emat[ino][jno][0][1] += w*vol*(lambda*dNdx[ino][0]*dNdx[jno][1]+myu*dNdx[jno][0]*dNdx[ino][1]);
        emat[ino][jno][0][2] += w*vol*(lambda*dNdx[ino][0]*dNdx[jno][2]+myu*dNdx[jno][0]*dNdx[ino][2]);
        emat[ino][jno][1][0] += w*vol*(lambda*dNdx[ino][1]*dNdx[jno][0]+myu*dNdx[jno][1]*dNdx[ino][0]);
        emat[ino][jno][1][1] += w*vol*(lambda*dNdx[ino][1]*dNdx[jno][1]+myu*dNdx[jno][1]*dNdx[ino][1]);
        emat[ino][jno][1][2] += w*vol*(lambda*dNdx[ino][1]*dNdx[jno][2]+myu*dNdx[jno][1]*dNdx[ino][2]);
        emat[ino][jno][2][0] += w*vol*(lambda*dNdx[ino][2]*dNdx[jno][0]+myu*dNdx[jno][2]*dNdx[ino][0]);
        emat[ino][jno][2][1] += w*vol*(lambda*dNdx[ino][2]*dNdx[jno][1]+myu*dNdx[jno][2]*dNdx[ino][1]);
        emat[ino][jno][2][2] += w*vol*(lambda*dNdx[ino][2]*dNdx[jno][2]+myu*dNdx[jno][2]*dNdx[ino][2]);
        const double dtmp1 = dNdx[ino][0]*dNdx[jno][0]+dNdx[ino][1]*dNdx[jno][1]+dNdx[ino][2]*dNdx[jno][2];
        emat[ino][jno][0][0] += w*vol*myu*dtmp1;
        emat[ino][jno][1][1] += w*vol*myu*dtmp1;
        emat[ino][jno][2][2] += w*vol*myu*dtmp1;
      }
    }
    for (unsigned int ino = 0; ino<10; ino++){
      eres[ino][0] += w*vol*rho*g_x*N[ino];
      eres[ino][1] += w*vol*rho*g_y*N[ino];
      eres[ino][2] += w*vol*rho*g_z*N[ino];
      for (unsigned int jno = 0; jno<10; jno++){
        eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1]+emat[ino][jno][0][2]*disp[jno][2];
        eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1]+emat[ino][jno][1][2]*disp[jno][2];
        eres[ino][2] -= emat[ino][jno][2][0]*disp[jno][0]+emat[ino][jno][2][1]*disp[jno][1]+emat[ino][jno][2][2]*disp[jno][2];
      }
    }
  }
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
  ddW_SolidLinear_Tet3D(&emat[0][0][0][0],
                        lambda, myu, vol, dldx, is_add, 3);
  if( !is_add ){
    for(int i=0;i<12;++i){ (&eres[0][0])[i] = 0.0; }
  }
  for (int ino = 0; ino<4; ino++){
    for (int jno = 0; jno<4; jno++){
      eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1]+emat[ino][jno][0][2]*disp[jno][2];
      eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1]+emat[ino][jno][1][2]*disp[jno][2];
      eres[ino][2] -= emat[ino][jno][2][0]*disp[jno][0]+emat[ino][jno][2][1]*disp[jno][1]+emat[ino][jno][2][2]*disp[jno][2];
    }
  }
}

DFM2_INLINE void delfem2::MakeMat_LinearSolid3D_Static_Q1(
    const double myu,
    const double lambda,
    const double rho,
    const double g_x,
    const double g_y,
    const double g_z,
    const double coords[8][3],
    const double disp[8][3],
    //
    double emat[8][8][3][3],
    double eres[8][3])
{
  const int nDegInt = 2;
  const int nInt = NIntLineGauss[nDegInt];
  const double (*Gauss)[2] = LineGauss[nDegInt];
  
  for(unsigned int i=0;i<8*8*3*3;i++){ *( &emat[0][0][0][0]+i) = 0.0; }
  for(unsigned int i=0;i<    8*3;i++){ *( &eres[0][0]      +i) = 0.0; }
  
  double vol = 0.0;
  for(int ir1=0;ir1<nInt;ir1++){
  for(int ir2=0;ir2<nInt;ir2++){
  for(int ir3=0;ir3<nInt;ir3++){
    const double r1 = Gauss[ir1][0];
    const double r2 = Gauss[ir2][0];
    const double r3 = Gauss[ir3][0];
    double detjac, detwei, dndx[8][3], an[8];
    ShapeFunc_Hex8(r1,r2,r3,coords,detjac,dndx,an);
    detwei = detjac*Gauss[ir1][1]*Gauss[ir2][1]*Gauss[ir3][1];
    vol += detwei;
    for(int ino=0;ino<8;ino++){
    for(int jno=0;jno<8;jno++){
      double dtmp1 = 0.0;
      for(int idim=0;idim<3;idim++){
        for(int jdim=0;jdim<3;jdim++){
          emat[ino][jno][idim][jdim] += detwei*( lambda*dndx[ino][idim]*dndx[jno][jdim]
                                                +myu*dndx[jno][idim]*dndx[ino][jdim] );
        }
        dtmp1 += dndx[ino][idim]*dndx[jno][idim];
      }
      for(int idim=0;idim<3;idim++){
        emat[ino][jno][idim][idim] += detwei*myu*dtmp1;
      }
    }
    }
    for(int ino=0;ino<8;ino++){
      eres[ino][0] += detwei*rho*g_x*an[ino];
      eres[ino][1] += detwei*rho*g_y*an[ino];
      eres[ino][2] += detwei*rho*g_z*an[ino];
    }
  }
  }
  }
  for (int ino = 0; ino<8; ino++){
  for (int jno = 0; jno<8; jno++){
    eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1]+emat[ino][jno][0][2]*disp[jno][2];
    eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1]+emat[ino][jno][1][2]*disp[jno][2];
    eres[ino][2] -= emat[ino][jno][2][0]*disp[jno][0]+emat[ino][jno][2][1]*disp[jno][1]+emat[ino][jno][2][2]*disp[jno][2];
  }
  }
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
    const double gamma_newmark,
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
  ddW_SolidLinear_Tet3D(&eKmat[0][0][0][0],
                        lambda,myu,
                        vol, dldx, false, 3);
  
  double eMmat[nno][nno][ndim][ndim];
  ddW_MassConsistentVal3D_Tet3D(
      &eMmat[0][0][0][0],
      rho,vol,false,3);
  
  // calc external force
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = vol*rho*g_x*0.25;
    eres[ino][1] = vol*rho*g_y*0.25;
    eres[ino][2] = vol*rho*g_z*0.25;
  }
  
  ////////////////////////////////////////////////////////////////////////////////////
  
  {	// calc coeff matrix for newmark-beta
    double dtmp1 = beta_newmark*dt*dt;
    for(int i=0;i<nno*nno*ndim*ndim;i++){
      (&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i]+dtmp1*(&eKmat[0][0][0][0])[i];
    }
  }
  
  // calc element redisual vector
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres[ino][0] -= eKmat[ino][jno][0][0]*disp[jno][0]+eKmat[ino][jno][0][1]*disp[jno][1]+eKmat[ino][jno][0][2]*disp[jno][2];
      eres[ino][1] -= eKmat[ino][jno][1][0]*disp[jno][0]+eKmat[ino][jno][1][1]*disp[jno][1]+eKmat[ino][jno][1][2]*disp[jno][2];
      eres[ino][2] -= eKmat[ino][jno][2][0]*disp[jno][0]+eKmat[ino][jno][2][1]*disp[jno][1]+eKmat[ino][jno][2][2]*disp[jno][2];
    }
    for(int jno=0;jno<nno;jno++){
      eres[ino][0] -= eMmat[ino][jno][0][0]*acc[jno][0]+eMmat[ino][jno][0][1]*acc[jno][1]+eMmat[ino][jno][0][2]*acc[jno][2];
      eres[ino][1] -= eMmat[ino][jno][1][0]*acc[jno][0]+eMmat[ino][jno][1][1]*acc[jno][1]+eMmat[ino][jno][1][2]*acc[jno][2];
      eres[ino][2] -= eMmat[ino][jno][2][0]*acc[jno][0]+eMmat[ino][jno][2][1]*acc[jno][1]+eMmat[ino][jno][2][2]*acc[jno][2];
    }
  }
  if( is_initial_iter ){
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eres[ino][0] -= dt*(eKmat[ino][jno][0][0]*velo[jno][0]+eKmat[ino][jno][0][1]*velo[jno][1]+eKmat[ino][jno][0][2]*velo[jno][2]);
        eres[ino][1] -= dt*(eKmat[ino][jno][1][0]*velo[jno][0]+eKmat[ino][jno][1][1]*velo[jno][1]+eKmat[ino][jno][1][2]*velo[jno][2]);
        eres[ino][2] -= dt*(eKmat[ino][jno][2][0]*velo[jno][0]+eKmat[ino][jno][2][1]*velo[jno][1]+eKmat[ino][jno][2][2]*velo[jno][2]);
      }
      for(int jno=0;jno<nno;jno++){
        eres[ino][0] -= 0.5*dt*dt*(eKmat[ino][jno][0][0]*acc[jno][0]+eKmat[ino][jno][0][1]*acc[jno][1]+eKmat[ino][jno][0][2]*acc[jno][2]);
        eres[ino][1] -= 0.5*dt*dt*(eKmat[ino][jno][1][0]*acc[jno][0]+eKmat[ino][jno][1][1]*acc[jno][1]+eKmat[ino][jno][1][2]*acc[jno][2]);
        eres[ino][2] -= 0.5*dt*dt*(eKmat[ino][jno][2][0]*acc[jno][0]+eKmat[ino][jno][2][1]*acc[jno][1]+eKmat[ino][jno][2][2]*acc[jno][2]);
      }
    }
  }
}


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
  double dldx[nno][ndim], zero_order_term[nno];
  TriDlDx(dldx, zero_order_term, coords[0],coords[1],coords[2]);
  
  for(int ino=0;ino<nno;ino++){
     for(int jno=0;jno<nno;jno++){
       emat[ino][jno][0][0] = area*(lambda+myu)*dldx[ino][0]*dldx[jno][0];
       emat[ino][jno][0][1] = area*(lambda*dldx[ino][0]*dldx[jno][1]+myu*dldx[jno][0]*dldx[ino][1]);
       emat[ino][jno][1][0] = area*(lambda*dldx[ino][1]*dldx[jno][0]+myu*dldx[jno][1]*dldx[ino][0]);
       emat[ino][jno][1][1] = area*(lambda+myu)*dldx[ino][1]*dldx[jno][1];
       const double dtmp1 = (dldx[ino][1]*dldx[jno][1]+dldx[ino][0]*dldx[jno][0])*area*myu;
       emat[ino][jno][0][0] += dtmp1;
       emat[ino][jno][1][1] += dtmp1;
     }
  }
  
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
    const double gamma_newmark,
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
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eKmat[ino][jno][0][0] = area*(lambda+myu)*dldx[ino][0]*dldx[jno][0];
      eKmat[ino][jno][0][1] = area*(lambda*dldx[ino][0]*dldx[jno][1]+myu*dldx[jno][0]*dldx[ino][1]);
      eKmat[ino][jno][1][0] = area*(lambda*dldx[ino][1]*dldx[jno][0]+myu*dldx[jno][1]*dldx[ino][0]);
      eKmat[ino][jno][1][1] = area*(lambda+myu)*dldx[ino][1]*dldx[jno][1];
      const double dtmp1 = (dldx[ino][1]*dldx[jno][1]+dldx[ino][0]*dldx[jno][0])*area*myu;
      eKmat[ino][jno][0][0] += dtmp1;
      eKmat[ino][jno][1][1] += dtmp1;
    }
  }
  
  double eMmat[nno][nno][ndim][ndim];
  {
    const double dtmp1 = area*rho*0.0833333333333333333333333;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat[ino][jno][0][0] = dtmp1;
        eMmat[ino][jno][0][1] = 0.0;
        eMmat[ino][jno][1][0] = 0.0;
        eMmat[ino][jno][1][1] = dtmp1;
      }
      eMmat[ino][ino][0][0] += dtmp1;
      eMmat[ino][ino][1][1] += dtmp1;
    }
  }
  
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
