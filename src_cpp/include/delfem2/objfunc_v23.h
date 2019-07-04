//
//  pbd_v23.h
//  delfem2
//
//  Created by Nobuyuki Umetani on 2019-06-08.
//

#ifndef OBJFUNC_v23_h
#define OBJFUNC_v23_h

void Energy_MIPS(double& E,
                 double dE[3][3],
                 double ddE[3][3][3][3],
                 const double c[3][3],
                 const double C[3][3]);


void PBD_Pre3D(std::vector<double>& aXYZt,
               double dt,
               const double gravity[3],
               const std::vector<double>& aXYZ,
               const std::vector<double>& aUVW,
               const std::vector<int>& aBCFlag);

void PBD_Post(std::vector<double>& aXYZ,
              std::vector<double>& aUVW,
              double dt,
              const std::vector<double>& aXYZt,
              const std::vector<int>& aBCFlag);

void PBD_Update_Const3_Point3_Dim3(std::vector<double>& aXYZt,
                                   const double m[3],
                                   const double C[3],
                                   const double dCdp[3][9],
                                   const int aIP[3]);
void PBD_Update_Const3(double* aXYZt,
                       const int np,
                       const int ndim,
                       const double* m,
                       const double* C,
                       const double* dCdp,
                       const int* aIP);

void PBD_ConstProj_Rigid2D(double* aXYt,
                                  double stiffness,
                                  const int* clstr_ind, int nclstr_ind,
                                  const int* clstr,     int nclstr0,
                                  const double* aXY0,   int nXY0);

void PBD_ConstProj_Rigid3D(double* aXYZt,
                                  double stiffness,
                                  const int* clstr_ind, int nclstr_ind,
                                  const int* clstr,     int nclstr0,
                                  const double* aXYZ0,   int nXYZ0);

void PBD_CdC_TriStrain2D3D(double C[3],
                       double dCdp[3][9],
                       const double P[3][2], // (in) undeformed triangle vertex positions
                       const double p[3][3]); // (in) deformed triangle vertex positions
void Check_CdC_TriStrain(const double P[3][2], // (in) undeformed triangle vertex positions
                         const double p[3][3]); // (in) deformed triangle vertex positions)

void PBD_ConstraintProjection_DistanceTri2D3D(double C[3],
                                              double dCdp[3][9],
                                              const double P[3][2], // (in) undeformed triangle vertex positions
                                              const double p[3][3]); // (in) deformed triangle vertex positions
void PBD_ConstraintProjection_EnergyStVK(double& C,
                                         double dCdp[9],
                                         const double P[3][2], // (in) undeformed triangle vertex positions
                                         const double p[3][3], // (in) deformed triangle vertex positions)
                                         const double lambda,
                                         const double myu);
void PBD_ConstraintProjection_DistanceTet(double C[6],
                                          double dCdp[6][12],
                                          const double P[4][3], // (in) undeformed triangle vertex positions
                                          const double p[4][3]); // (in) deformed triangle vertex positions

void Check_ConstraintProjection_DistanceTri2D3D(const double P[3][2], // (in) undeformed triangle vertex positions
                                                const double p[3][3]); // (in) deformed triangle vertex positions)
void Check_ConstraintProjection_EnergyStVK(const double P[3][2], // (in) undeformed triangle vertex positions
                                           const double p[3][3], // (in) deformed triangle vertex positions)
                                           const double lambda,
                                           const double myu);

void PBD_CdC_QuadBend(double C[3],
                      double dCdp[3][12],
                      const double P[4][3],
                      const double p[4][3]);


#endif /* pbd_v23_h */
