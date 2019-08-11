/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef fem_ematrix_h
#define fem_ematrix_h

// derivative of a shape function of a triangle and constant compornent
void TriDlDx
(double dldx[][2], double const_term[],
 const double p0[], const double p1[], const double p2[]);

// -[\alpha]\nabla^2[value] = [source]
void MakeMat_Poisson2D_P1
(const double alpha, const double source,
 const double coords[3][2],
 const double value[3],
 double eres[3],
 double emat[][3]);

// [\rho][velo] - [\alpha]\nabla^2[value] = [source]
void MakeMat_Diffusion2D_P1
(const double alpha, const double source,
 const double dt, const double gamma, const double rho,
 const double coords[3][2], 
 const double value[3], const double velo[3],
 double eres[3],
 double emat[3][3]);

void MakeMat_LinearSolid2D_Static_P1
(const double myu, const double lambda,
 const double rho, const double g_x, const double g_y,
 const double disp[3][2],
 const double coords[3][2],
 double eres[3][2],
 double emat[][3][2][2]);

void MakeMat_LinearSolid2D_Dynamic_P1
(const double myu, const double lambda,
 const double rho, const double g_x, const double g_y,
 const double dt_timestep, const double gamma_newmark,  const double beta_newmark,
 const double disp[3][2], const double velo[3][2], const double acc[3][2],
 const double coords[3][2],
 double eres[3][2],
 double emat[3][3][2][2],
 bool is_initial);

void MakeMat_Stokes2D_Static_P1
(double alpha, double g_x, double g_y,
 const double coords[][2],
 const double velo_press[3][3],
 double emat[3][3][3][3],
 double eres[3][3]);

void MakeMat_Stokes2D_Dynamic_P1
(double alpha, double rho, double g_x, double g_y,
 const double dt_timestep, const double gamma_newmark,
 const double coords[][2],
 const double velo_press[3][3], const double acc_apress[3][3],
 double emat[3][3][3][3],
 double eres[3][3]);

void MakeMat_NavierStokes2D_NonStatic_Newmark_P1P1
(double dt, double gamma,
 double rho, double myu, double g_x, double g_y,
 const double coords[][2],
 const double velo[][2], const double acc[][2], const double press[], const double apress[],
 double eres_u[3][2], double eres_p[3],
 double eCmat_uu[][3][2][2], double eCmat_up[][3][2], double eCmat_pu[][3][2], double eCmat_pp[][3],
 double eMmat_uu[][3][2][2], double eMmat_pu[][3][2]);

void MakeMat_NavierStokes2D_Dynamic_P1
(double myu, double rho, double g_x, double g_y,
 const double dt_timestep, const double gamma_newmark,
 const double coords[][2],
 const double velo_press[3][3], const double acc_apress[3][3],
 double emat[3][3][3][3],
 double eres[3][3]);

// compute energy and its 1st and 2nd derivative for cloth bending
void WdWddW_Bend
(double& W,  // (out) strain energy
 double dW[4][3], // (out) 1st derivative of energy
 double ddW[4][4][3][3], // (out) 2nd derivative of energy
 ////
 const double C[4][3], // (in) undeformed triangle vertex positions
 const double c[4][3], // (in) deformed triangle vertex positions
 double stiff);

void WdWddW_CST
(double& W, // (out) energy
 double dW[3][3], // (out) 1st derivative of energy
 double ddW[3][3][3][3], // (out) 2nd derivative of energy
 ////
 const double C[3][3], // (in) undeformed triangle vertex positions
 const double c[3][3], // (in) deformed triangle vertex positions
 const double lambda, // (in) Lame's 1st parameter
 const double myu);   // (in) Lame's 2nd parameter

// compute energy and its 1st and 2nd derivative for contact against object

class CInput_Contact
{
public:
  virtual double penetrationNormal(double& nx, double& ny, double& nz,
                                   double px, double py, double pz) const = 0;
};

void WdWddW_Contact
(double& W,  // (out) energy
 double dW[3], // (out) 1st derivative of energy
 double ddW[3][3], // (out) 2nd derivative of energy
 ////
 const double c[3], // (in) deformed triangle vertex positions
 double stiff_contact,
 double contact_clearance,
 const CInput_Contact& input);



//////////////////////////////

void MakeMat_Poisson3D_P1
(const double alpha, const double source,
 const double coords[4][3],
 const double value[4],
 double eres[4],
 double emat[][4]);

void MakeMat_Diffusion3D_P1
(const double alpha, const double source,
 const double dt_timestep, const double gamma_newmark, const double rho,
 const double coords[4][3],
 const double value[4], const double velo[4],
 double eres[4],
 double emat[4][4]);

void MakeMat_LinearSolid3D_Static_P1
(const double myu, const double lambda,
 const double rho, const double g_x, const double g_y, const double g_z,
 const double coords[4][3],
 const double disp[4][3],
 ////
 double emat[4][4][3][3],
 double eres[4][3]);

void MakeMat_LinearSolid3D_Static_Q1
(const double myu, const double lambda,
 const double rho, const double g_x, const double g_y, const double g_z,
 const double coords[8][3],
 const double disp[8][3],
 ////
 double emat[8][8][3][3],
 double eres[8][3]);

void MakeMat_LinearSolid3D_Dynamic_P1
(const double myu, const double lambda,
 const double rho, const double g_x, const double g_y, const double g_z,
 const double dt, const double gamma_newmark,  const double beta_newmark,
 const double disp[4][3], const double velo[4][3], const double acc[4][3],
 const double coords[4][3],
 double eres[4][3],
 double emat[4][4][3][3],
 bool is_initial);

void MakeMat_Stokes3D_Static_P1
(double alpha, double g_x, double g_y, double g_z,
 const double coords[4][3],
 const double velo_press[4][4],
 double emat[4][4][4][4],
 double eres[4][4]);

void MakeMat_Stokes3D_Dynamic_P1
(double alpha, double rho, double g_x, double g_y, double g_z,
 const double dt_timestep, const double gamma_newmark,
 const double coords[4][3],
 const double velo_press[4][4], const double acc_apress[4][4],
 double emat[4][4][4][4],
 double eres[4][4]);

void MakeMat_NavierStokes3D_Dynamic_P1
(double myu, double rho, double g_x, double g_y, double g_z,
 const double dt_timestep, const double gamma_newmark,
 const double coords[4][3],
 const double velo_press[4][4], const double acc_apress[4][4],
 double emat[4][4][4][4],
 double eres[4][4]);


#endif /* fem_ematrix_h */
