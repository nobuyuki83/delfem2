/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_EMAT_H
#define DFM2_EMAT_H

#include "delfem2/dfm2_inline.h"
#include <complex>

namespace delfem2 {

/**
 * @brief derivative of a shape function of a triangle and constant compornent
 */
void TriDlDx(
    double dldx[3][2],
    double const_term[3],
    const double p0[2],
    const double p1[2],
    const double p2[2]);

void ddW_MassConsistentVal3D_Tet3D(
    double* eMmat,
    double rho, double vol,
    bool is_add,
    unsigned int nstride = 3);

void ddW_SolidLinear_Tet3D(
    double* eKmat,
    double lambda, double myu,
    double vol, double dldx[4][3],
    bool is_add,
    unsigned int nstride = 3);

// -------------------------------------------------------------------

// -[\alpha]\nabla^2[value] = [source]
void EMat_Poisson_Tri2D(
    double eres[3],
    double emat[3][3],
    const double alpha, const double source,
    const double coords[3][2],
    const double value[3]);

void EMat_Helmholtz_Tri2D(
    std::complex<double> eres[3],
    std::complex<double> emat[][3],
    const double wave_length,
    const double coords[3][2],
    const std::complex<double> value[3]);

void EMat_SommerfeltRadiationBC_Line2D(
    std::complex<double> eres[2],
    std::complex<double> emat[2][2],
    double wave_length,
    const double P[2][2],
    const std::complex<double> val[2]);

// [\rho][velo] - [\alpha]\nabla^2[value] = [source]
void EMat_Diffusion_Tri2D(
    double eres[3],
    double emat[3][3],
    const double alpha, const double source,
    const double dt, const double gamma, const double rho,
    const double coords[3][2],
    const double value[3], const double velo[3]);

void EMat_SolidStaticLinear_Tri2D(
    double eres[3][2],
    double emat[][3][2][2],
    const double myu, const double lambda,
    const double rho, const double g_x, const double g_y,
    const double disp[3][2],
    const double coords[3][2]);

void EMat_SolidDynamicLinear_Tri2D(
    double eres[3][2],
    double emat[3][3][2][2],
    const double myu, const double lambda,
    const double rho, const double g_x, const double g_y,
    const double dt_timestep,
    const double gamma_newmark,  const double beta_newmark,
    const double disp[3][2],
    const double velo[3][2],
    const double acc[3][2],
    const double coords[3][2],
    bool is_initial);

DFM2_INLINE void stress_LinearSolid_TetP2
 (double stress[3][3],
  const double l0, const double l1, const double l2, const double l3,
  const double vol, const double lambda, const double myu,
  const double g_x, const double g_y, const double g_z,
  const double dldx[4][3],
  const double disp[10][3]);

void EMat_Stokes2D_Static_P1(
    double alpha, double g_x, double g_y,
    const double coords[][2],
    const double velo_press[3][3],
    double emat[3][3][3][3],
    double eres[3][3]);

void MakeMat_Stokes2D_Static_P1P1
 (double alpha, double g_x, double g_y,
  const double coords[][2],
  const double velo[3][2], const double press[3],
  double emat_uu[][3][2][2], double emat_up[][3][2], double emat_pu[][3][2], double emat_pp[][3],
  double eres_u[][2], double eres_p[3]);

void EMat_Stokes2D_Dynamic_P1(
    double alpha, double rho, double g_x, double g_y,
    const double dt_timestep, const double gamma_newmark,
    const double coords[][2],
    const double velo_press[3][3], const double acc_apress[3][3],
    double emat[3][3][3][3],
    double eres[3][3]);

DFM2_INLINE void MakeMat_Stokes2D_Dynamic_Newmark_P1P1
 (double alpha, double rho, double g_x, double g_y,
  const double dt_timestep, const double gamma_newmark,
  const double coords[3][2],
  const double velo[3][2], const double press[3], const double acc[3][2], const double apress[3],
  double emat_uu[3][3][2][2], double emat_up[][3][2], double emat_pu[][3][2], double emat_pp[][3],
  double eres_u[3][2], double eres_p[3]);

DFM2_INLINE void MakeMat_NavierStokes2D_Dynamic_Newmark_P1P1
 (double rho, double myu, double g_x, double g_y,
  double dt, double gamma,
  const double coords[][2],
  const double velo[][2], const double press[], const double acc[][2], const double apress[],
  double emat_uu[][3][2][2], double emat_up[][3][2], double emat_pu[][3][2], double emat_pp[][3],
  double eres_u[3][2], double eres_p[3]);

void EMat_NavierStokes2D_NonStatic_Newmark_P1P1(
    double dt, double gamma,
    double rho, double myu, double g_x, double g_y,
    const double coords[][2],
    const double velo[][2], const double acc[][2], const double press[], const double apress[],
    double eres_u[3][2], double eres_p[3],
    double eCmat_uu[][3][2][2], double eCmat_up[][3][2], double eCmat_pu[][3][2], double eCmat_pp[][3],
    double eMmat_uu[][3][2][2], double eMmat_pu[][3][2]);

void EMat_NavierStokes2D_Dynamic_P1(
    double myu, double rho, double g_x, double g_y,
    const double dt_timestep, const double gamma_newmark,
    const double coords[][2],
    const double velo_press[3][3], const double acc_apress[3][3],
    double emat[3][3][3][3],
    double eres[3][3]);

/**
 * @brief compute energy and its 1st and 2nd derivative for cloth bending
 * @param[out] W strain energy
 * @param[out] dW 1st derivative of energy
 * @param[out] ddW 2nd derivative of energy
 * @param[in] C undeformed quad vertex positions
 * @param[in] c deformed quad vertex positions
 */
void WdWddW_Bend(
    double& W,
    double dW[4][3],
    double ddW[4][4][3][3],
    //
    const double C[4][3],
    const double c[4][3],
    double stiff);

void WdWddW_CST(
    double& W, // (out) energy
    double dW[3][3], // (out) 1st derivative of energy
    double ddW[3][3][3][3], // (out) 2nd derivative of energy
    //
    const double C[3][3], // (in) undeformed triangle vertex positions
    const double c[3][3], // (in) deformed triangle vertex positions
    const double lambda, // (in) Lame's 1st parameter
    const double myu);   // (in) Lame's 2nd parameter

/**
 * @brief compute energy and its 1st and 2nd derivative for contact against object
 */
class CInput_Contact
{
public:
  virtual double penetrationNormal(
      double& nx, double& ny, double& nz,
      double px, double py, double pz) const = 0;
};

void WdWddW_Contact(
    double& W,  // (out) energy
    double dW[3], // (out) 1st derivative of energy
    double ddW[3][3], // (out) 2nd derivative of energy
    //
    const double c[3], // (in) deformed triangle vertex positions
    double stiff_contact,
    double contact_clearance,
    const CInput_Contact& input);

// --------------------------------------------------------------------------

void EMat_Poisson_Tet3D(double eres[4],
    double emat[][4],
    const double alpha, const double source,
    const double coords[4][3],
    const double value[4]);

void EMat_Diffusion_Newmark_Tet3D(double eres[4],
    double emat[4][4],
    const double alpha, const double source,
    const double dt_timestep, const double gamma_newmark, const double rho,
    const double coords[4][3],
    const double value[4], const double velo[4]);

void EMat_SolidLinear_Static_Tet(double emat[4][4][3][3],
    double eres[4][3],
    const double myu, const double lambda,
    const double coords[4][3],
    const double disp[4][3],
    bool is_add);

void MakeMat_LinearSolid3D_Static_Q1(const double myu, const double lambda,
    const double rho, const double g_x, const double g_y, const double g_z,
    const double coords[8][3],
    const double disp[8][3],
    //
    double emat[8][8][3][3],
    double eres[8][3]);

DFM2_INLINE void matRes_LinearSolid_TetP2
 (double emat[10][10][3][3],
  double eres[10][3],
  const double vol, const double lambda, const double myu,
  const double g_x, const double g_y, const double g_z,
  const double rho,
  const double dldx[4][3],
  const double disp[10][3]);

void EMat_SolidLinear_NewmarkBeta_MeshTet3D(
    double eres[4][3],
    double emat[4][4][3][3],
    const double myu, const double lambda,
    const double rho, const double g_x, const double g_y, const double g_z,
    const double dt, const double gamma_newmark,  const double beta_newmark,
    const double disp[4][3], const double velo[4][3], const double acc[4][3],
    const double coords[4][3],
    bool is_initial);

void MakeMat_Stokes3D_Static_P1(
    double alpha, double g_x, double g_y, double g_z,
    const double coords[4][3],
    const double velo_press[4][4],
    double emat[4][4][4][4],
    double eres[4][4]);

void MakeMat_Stokes3D_Static_P1P1
 (double alpha, double g_x, double g_y, double g_z,
  const double coords[4][3],
  const double velo[4][3], const double press[4],
  double emat_uu[4][4][3][3], double emat_up[4][4][3], double emat_pu[4][4][3], double emat_pp[4][4],
  double eres_u[4][3], double eres_p[4]);

void MakeMat_Stokes3D_Dynamic_P1(
    double alpha, double rho, double g_x, double g_y, double g_z,
    const double dt_timestep, const double gamma_newmark,
    const double coords[4][3],
    const double velo_press[4][4], const double acc_apress[4][4],
    double emat[4][4][4][4],
    double eres[4][4]);

void MakeMat_Stokes3D_Dynamic_Newmark_P1P1
 (double alpha, double rho, double g_x, double g_y, double g_z,
  const double dt_timestep, const double gamma_newmark,
  const double coords[4][3],
  const double velo[4][3], const double press[4], const double acc[4][3], const double apress[4],
  double emat_uu[4][4][3][3], double emat_up[4][4][3], double emat_pu[4][4][3], double emat_pp[4][4],
  double eres_u[4][3], double eres_p[4]);

void MakeMat_NavierStokes3D_Dynamic_P1(
    double myu, double rho, double g_x, double g_y, double g_z,
    const double dt_timestep, const double gamma_newmark,
    const double coords[4][3],
    const double velo_press[4][4], const double acc_apress[4][4],
    double emat[4][4][4][4],
    double eres[4][4]);

void MakeMat_NavierStokes3D_Dynamic_Newmark_P1P1
 (double rho, double myu, double g_x, double g_y, double g_z,
  double dt, double gamma,
  const double coords[4][3],
  const double velo[4][3], const double press[4], const double acc[4][3], const double apress[4],
  double emat_uu[4][4][3][3], double emat_up[4][4][3], double emat_pu[4][4][3], double emat_pp[4][4],
  double eres_u[4][3], double eres_p[4]);

DFM2_INLINE void MakeMat_PlateBendingDKT
 (double emat_ww[3][3],
  double emat_wr[3][3][2],
  double emat_rw[3][3][2],
  double emat_rr[3][3][2][2],
  double eres_w[3],
  double eres_r[3][2],
  const double young, const double poisson, const double thickness,
  const double coord[][2], const double w[], const double rot[][2]);

void WdWddW_PlateBendingMITC3(double& W,
    double dW[3][3],
    double ddW[3][3][3][3],
    const double C[3][2], // initial XY position
    const double u[3][3], // z displacement + xy axis rotation
    double thk,
    double lambda,
    double myu);
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/emat.cpp"
#endif
  
#endif /* fem_ematrix_h */
