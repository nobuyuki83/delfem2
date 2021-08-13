/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/fempoisson.h"

// --------------------------------------------------------
// below tri

DFM2_INLINE void delfem2::EMat_Poisson_Tri2D(
    double eres[3],
    double emat[3][3],
    const double alpha,
    const double source,
    const double coords[3][2],
    const double value[3])
{
  constexpr int nno = 3;
  constexpr int ndim = 2;
  //
  const double area = femutil::TriArea2D(coords[0], coords[1], coords[2]);
  //
  double dldx[nno][ndim], const_term[nno];
  TriDlDx(dldx, const_term, coords[0], coords[1], coords[2]);
  //
  for (int i = 0; i<9; ++i){ (&emat[0][0])[i] = 0.0; }
  for (int ino = 0; ino<nno; ino++){
    for (int jno = 0; jno<nno; jno++){
      emat[ino][jno] = alpha*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
    }
  }
  eres[0] = 0;  eres[1] = 0;  eres[2] = 0;
  for (int ino = 0; ino<nno; ino++){
    eres[ino] = source*area*0.33333333333333333;
  }
  for (int ino = 0; ino<nno; ino++){
    for (int jno = 0; jno<nno; jno++){
      eres[ino] -= emat[ino][jno]*value[jno];
    }
  }
}


DFM2_INLINE void delfem2::EMat_Diffusion_Tri2D(
    double eres[3],
    double emat[3][3],
    const double alpha,
    const double source,
    const double dt_timestep,
    const double gamma_newmark,
    const double rho,
    const double coords[3][2],
    const double value[3],
    const double velo[3])
{
  constexpr int nno = 3;
  constexpr int ndim = 2;

  const double area = femutil::TriArea2D(coords[0],coords[1],coords[2]);
  double dldx[nno][ndim], const_term[nno];
  TriDlDx(dldx,const_term,coords[0],coords[1],coords[2]);

  // --------------------------

  for (int i = 0; i<9; ++i){ (&emat[0][0])[i] = 0.0; }
  double eCmat[nno][nno];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat[ino][jno] = alpha*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
    }
  }
  double eMmat[nno][nno];
  {
    const double dtmp1 = rho*area*0.08333333333333333;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat[ino][jno] = dtmp1;
      }
      eMmat[ino][ino] += dtmp1;
    }
  }

  eres[0] = 0;  eres[1] = 0;  eres[2] = 0;
  for(int ino=0;ino<nno;ino++){
    eres[ino] = source*area*0.333333333333333333;
  }

  {
    const double dtmp1 = gamma_newmark*dt_timestep;
    for(int i=0;i<nno*nno;i++){
      (&emat[0][0])[i] = (&eMmat[0][0])[i]+dtmp1*(&eCmat[0][0])[i];
    }
  }
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres[ino]	-= eCmat[ino][jno]*(value[jno]+dt_timestep*velo[jno]) + eMmat[ino][jno]*velo[jno];
    }
  }

}

// above: tri
// ------------------------------------------------------------
// below: quad

void delfem2::EMat_Poission2_QuadOrth(
    double emat[4][4],
    double lx,
    double ly)
{
  const double lxy = ly / lx;
  const double lyx = lx / ly;
  emat[0][0] = emat[1][1] = emat[2][2] = emat[3][3] = (+1 / 3.0) * (lxy + lyx);
  emat[0][1] = emat[1][0] = emat[2][3] = emat[3][2] = (+1 / 6.0) * (-2 * lxy + lyx);
  emat[0][2] = emat[2][0] = emat[1][3] = emat[3][1] = (-1 / 6.0) * (lxy + lyx);
  emat[0][3] = emat[3][0] = emat[1][2] = emat[2][1] = (+1 / 6.0) * (lxy - 2 * lyx);
}


void delfem2::EMat_Poisson2_QuadOrth_GaussInt(
    double emat[4][4],
    double lx,
    double ly,
    unsigned int ngauss)
{
  namespace lcl = delfem2::femutil;
  for (unsigned int i = 0; i < 16; ++i) { (&emat[0][0])[i] = 0.0; }
  unsigned int nw = NIntLineGauss[ngauss];
  for (unsigned int iw = 0; iw < nw; ++iw) {
    for (unsigned int jw = 0; jw < nw; ++jw) {
      const double w = lx * ly * 0.25 * LineGauss<double>[ngauss][iw][1] * LineGauss<double>[ngauss][jw][1];
      const double x1 = (1 - LineGauss<double>[ngauss][iw][0]) * 0.5;
      const double y1 = (1 - LineGauss<double>[ngauss][jw][0]) * 0.5;
      const double x2 = 1 - x1;
      const double y2 = 1 - y1;
      // u = u1x1y1 + u2x2y1 + u3x2y2 + u4x1y2
      const double dldx[4][2] = {
          {-y1 / lx, -x1 / ly},
          {+y1 / lx, -x2 / ly},
          {+y2 / lx, +x2 / ly},
          {-y2 / lx, +x1 / ly}};
      for (unsigned int in = 0; in < 4; ++in) {
        for (unsigned int jn = 0; jn < 4; ++jn) {
          emat[in][jn] += w * (dldx[in][0] * dldx[jn][0] + dldx[in][1] * dldx[jn][1]);
        }
      }
    }
  }
}

// above: quad
// -------------------------------------------------------------
// below: tet

DFM2_INLINE void delfem2::EMat_Poisson_Tet3D(
    double eres[4],
    double emat[4][4],
    const double alpha,
    const double source,
    const double coords[4][3],
    const double value[4])
{
  constexpr int nno = 4;
  constexpr int ndim = 3;
  //
  const double area = femutil::TetVolume3D(coords[0], coords[1], coords[2], coords[3]);
  //
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx, const_term, coords[0], coords[1], coords[2], coords[3]);

  for (int i = 0; i<16; ++i){ (&emat[0][0])[i] = 0.0; }
  for (int ino = 0; ino<nno; ino++){
    for (int jno = 0; jno<nno; jno++){
      emat[ino][jno] = alpha*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]+dldx[ino][2]*dldx[jno][2]);
    }
  }
  eres[0] = 0;  eres[1] = 0;  eres[2] = 0;  eres[3] = 0;
  for (int ino = 0; ino<nno; ino++){
    eres[ino] = source*area*0.25;
  }
  for (int ino = 0; ino<nno; ino++){
    for (int jno = 0; jno<nno; jno++){
      eres[ino] -= emat[ino][jno]*value[jno];
    }
  }
}

DFM2_INLINE void delfem2::EMat_Diffusion_Newmark_Tet3D(
    double eres[4],
    double emat[4][4],
    const double alpha,
    const double source,
    const double dt_timestep,
    const double gamma_newmark,
    const double rho,
    const double coords[4][3],
    const double value[4],
    const double velo[4])
{
  constexpr int nno = 4;
  constexpr int ndim = 3;

  const double vol = femutil::TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx,const_term,coords[0],coords[1],coords[2],coords[3]);
  
  // ----------------------

  for (int i=0; i<16; ++i){ (&emat[0][0])[i] = 0.0; }
  double eCmat[nno][nno];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat[ino][jno] = alpha*vol*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]+dldx[ino][2]*dldx[jno][2]);
    }
  }
  double eMmat[nno][nno];
  {
    const double dtmp1 = rho*vol*0.05;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat[ino][jno] = dtmp1;
      }
      eMmat[ino][ino] += dtmp1;
    }
  }

  eres[0] = 0;  eres[1] = 0;  eres[2] = 0;  eres[3] = 0;
  for(int ino=0;ino<nno;ino++){
    eres[ino] = source*vol*0.25;
  }
  {
    const double dtmp1 = gamma_newmark*dt_timestep;
    for(int i=0;i<nno*nno;i++){
      (&emat[0][0])[i] = (&eMmat[0][0])[i]+dtmp1*(&eCmat[0][0])[i];
    }
  }
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres[ino]	-= eCmat[ino][jno]*(value[jno]+dt_timestep*velo[jno]) + eMmat[ino][jno]*velo[jno];
    }
  }
}




