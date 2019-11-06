/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


/**
 * @details this file should not depends on anything other than std::vector
 */

#ifndef DFM2_MSH_H
#define DFM2_MSH_H

#include <vector>

// -----------------
// work on points

void updateMinMaxXYZ(double& x_min, double& x_max,
                     double& y_min, double& y_max,
                     double& z_min, double& z_max,
                     double x, double y, double z);
void GetCenterWidth(double& cx, double& cy, double& cz,
                    double& wx, double& wy, double& wz,
                    const int nXYZ, const double* paXYZ);
void GetCenterWidth(double& cx, double& cy, double& cz,
                    double& wx, double& wy, double& wz,
                    const std::vector<double>& aXYZ);
void GetCenterWidth(double cw[6],
                    const std::vector<double>& aXYZ);
void GetCenterWidthGroup(double& cx, double& cy, double& cz,
                         double& wx, double& wy, double& wz,
                         const std::vector<double>& aXYZ,
                         const std::vector<int>& aElem,
                         const int nnoel,
                         int igroup,
                         const std::vector<int>& aIndGroup);
void GetCenterWidthGroup(double& cx, double& cy, double& cz,
                         double& wx, double& wy, double& wz,
                         const std::vector<double>& aXYZ,
                         const std::vector<int>& aElemInd,
                         const std::vector<int>& aElem,
                         int igroup,
                         const std::vector<int>& aIndGroup);
void GetCenterWidth3DGroup(double cw[6],
                           ////
                           const std::vector<double>& aXYZ,
                           const std::vector<int>& aElemInd,
                           const std::vector<int>& aElem,
                           int igroup,
                           const std::vector<int>& aIndGroup);
// local coordinate
void GetCenterWidthLocal(double& lcx, double& lcy, double& lcz,
                         double& lwx, double& lwy, double& lwz,
                         const std::vector<double>& aXYZ,
                         const double lex[3],
                         const double ley[3],
                         const double lez[3]);
void MinMaxXYZ(double mm[6],
               const std::vector<double>& aXYZ);


// ------------------------------------------

namespace delfem2{

/**
 * @brief rotate with the Bryant angle (in the  order of XYZ) around the origin.
 * @details the angles are in the radian.
 */
void Rotate(std::vector<double>& aXYZ,
            double radx, double rady, double radz);
void Translate(std::vector<double>& aXYZ,
               double tx, double ty, double tz);
void Translate(double tx, double ty, double tz,
               const unsigned int nnode_, double* pXYZs_);
void Translate(std::vector<double>& aXY,
               double tx, double ty);
void Scale(std::vector<double>& aXYZ,
           double s);
void Scale(double s,
           const unsigned int nnode_, double* pXYZs_);
double Size(const std::vector<double>& aXYZ);
void Normalize(std::vector<double>& aXYZ,
               double s = 1.0);
  
}

// ------------------------------------------------------------

void CenterOfGravity(double& cgx, double& cgy, double& cgz,
                     const std::vector<double>& aXYZ);
void CenterOfGravity_Tri(double& cgx, double& cgy, double& cgz,
                         int itri,
                         const std::vector<double>& aXYZ,
                         const std::vector<int>& aTri);
double CenterOfGravity_TriMsh3DFlg_Shell(double& cgx, double& cgy, double& cgz,
                                         const std::vector<double>& aXYZ,
                                         const std::vector<int>& aTri,
                                         int iflg,
                                         const std::vector<int>& aFlg);
void CenterOfGravity_Shell(double& cgx, double& cgy, double& cgz,
                           const std::vector<double>& aXYZ,
                           const std::vector<int>& aTri);
void CenterOfGravity_Solid(double& cgx, double& cgy, double& cgz,
                           const std::vector<double>& aXYZ,
                           const std::vector<int>& aTri);
void CenterOfGravity_Tet(double& v_tot,
                         double& cgx, double& cgy, double& cgz,
                         const std::vector<double>& aXYZC,
                         const std::vector<int>& aTetC);

// -------------------------------------------------------------

void RemoveUnreferencedPoints_MeshElem(std::vector<double>& aXYZOut,
                                       std::vector<unsigned int>& aElemOut,
                                       std::vector<int>& aMapInOut,
                                       int ndim,
                                       const std::vector<double>& aXYZIn,
                                       const std::vector<unsigned int>& aElemIn);
void Normal_MeshTri3D(double* aNorm,
                      const double* aXYZ, int nXYZ,
                      const unsigned int* aTri, int nTri);
void Quality_MeshTri2D(double& max_aspect, double& min_area,
                       const double* aXY,
                       const unsigned int* aTri, int nTri);

// ----------------------
// set primitive mesh

void SetTopology_ExtrudeTri2Tet(unsigned int* aTet,
                                int nXY,
                                const unsigned int* aTri, int nTri,
                                int nlayer);
void ExtrudeTri2Tet(int nlayer, double h,
                    std::vector<double>& aXYZ,
                    std::vector<unsigned int>& aTet,
                    const std::vector<double>& aXY,
                    const std::vector<unsigned int>& aTri);

// -----------------------------------------------------
// considering moving these file to other locations

void makeSolidAngle(std::vector<double>& aSolidAngle,
                    const std::vector<double>& aXYZ,
                    const std::vector<int>& aTri);
void makeSolidAngle(std::vector<double>& aSolidAngle,
                    const std::vector<double>& aXYZ,
                    const std::vector<int>& aTri,
                    const std::vector<double>& aNorm,
                    std::vector<int>& elsup_ind,
                    std::vector<int>& elsup);

/**
 * @brief Compute mass of the points (lumped mass) for 3D tet mesh
 * @details aMassMatrixLumped need to be allocated before in the size of nXY
 * this is here because both "fem" and "pbd" use this function
 */
void MassPoint_Tet3D(double* aMassMatrixLumped,
                      double rho,
                      const double* aXYZ, int nXYZ,
                      const unsigned int* aTet, int nTet);

/**
 * @brief Compute mass of the points (lumped mass) for 2D triangle mesh
 * @param aMassMatrixLumped (out) this need to be allocated before in the size of nXY
 * @details this is here because both "fem" and "pbd" use this function
 */
void MassPoint_Tri2D(double* aMassMatrixLumped,
                     //
                     double rho,
                     const double* aXY, int nXY,
                     const unsigned int* aTri, int nTri);


void LaplacianSmoothing(std::vector<double>& aXYZ,
                        const std::vector<int>& aTri,
                        const std::vector<int>& elsup_ind,
                        const std::vector<int> elsup);

// ---------------------------------------------------------

void SubdivisionPoints_QuadCatmullClark(std::vector<double>& aXYZ1,
                                        //
                                        const std::vector<unsigned int>& aQuad1,
                                        const std::vector<int>& aEdgeFace0,
                                        const std::vector<int>& psupIndQuad0,
                                        const std::vector<int>& psupQuad0,
                                        const unsigned int* aQuad0, int nQuad0,
                                        const double* aXYZ0, int nXYZ0);
void SubdivisionPoints_Quad(std::vector<double>& aXYZ1,
                            //
                            const std::vector<int>& aQuad1,
                            const std::vector<int>& aEdgeFace0,
                            const std::vector<int>& psupIndQuad0,
                            const std::vector<int>& psupQuad0,
                            const std::vector<int>& aQuad0,
                            const std::vector<double>& aXYZ0);

void SubdivisionPoints_Hex(std::vector<double>& aXYZ1,
                           //
                           const std::vector<int>& psupHex0,
                           const std::vector<int>& psupIndHex0,
                           const std::vector<unsigned int>& aQuadHex0,
                           const unsigned int* aHex0, int nHex0,
                           const double* aXYZ0, int nXYZ0);


#endif
