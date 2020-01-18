/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


/**
 * @details this file should not depends on anything except for  std::vector
 */

#ifndef DFM2_MSH_H
#define DFM2_MSH_H

#include <vector>

// -----------------
// work on points

namespace delfem2{

/**
 * @brief update minimum and maximum coordinates
 * @details implemented for "float" and "double"
 */
template<typename T>
void updateMinMaxXYZ(T& x_min, T& x_max,
                     T& y_min, T& y_max,
                     T& z_min, T& z_max,
                     T x, T y, T z);

/**
 * @param bb3 (out) bounding box in the order of <minx, miny, minz, maxx, maxy, maxz>
 * @param aXYZ (in) array of 3D coordinates of points
 * @param nXYZ (in) number of points
 * @details implemented for "float" and "double"
 */
template<typename T>
void Min3Max3_Points3(
    T min3[3],
    T max3[3],
    const T* aXYZ,
    const unsigned int nXYZ);


void MinMaxXYZ(
    double mm[6],
    const std::vector<double>& aXYZ);


// center & width
void GetCenterWidth(double& cx, double& cy, double& cz,
                    double& wx, double& wy, double& wz,
                    const int nXYZ, const double* paXYZ);
void CenterWidth_Points3D(double& cx, double& cy, double& cz,
                          double& wx, double& wy, double& wz,
                          const std::vector<double>& aXYZ);
void CenterWidth_Points3D(double cw[6],
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
                           //
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



// ------------------------------------------


/**
 * @brief rotate with the Bryant angle (in the  order of XYZ) around the origin.
 * @details the angles are in the radian.
 */
void Rotate_Points3D(std::vector<double>& aXYZ,
                     double radx, double rady, double radz);
  
template <typename T>
void Translate_Points3(std::vector<T>& aXYZ,
                       T tx, T ty, T tz);
template <typename T>
void Translate_Points3(T* pXYZs_,
                       const unsigned int nnode_,
                       T tx, T ty, T tz);
template <typename T>
void Translate_Points2(std::vector<T>& aXY,
                       T tx, T ty);
  
template <typename T>
void Scale_PointsX(std::vector<T>& aXYZ,
                   T s);
template <typename T>
void Scale_Points3(T* pXYZs_,
                   const unsigned int nnode_,
                   T s);
  
double Size_Points3D_LongestAABBEdge(const std::vector<double>& aXYZ);
void Normalize_Points3D(std::vector<double>& aXYZ,
                        double s = 1.0);
  
/**
 * @details implemented for "float" and "double"
 */
template <typename T>
void CG_Point3 (T cg[3],
                const std::vector<T>& aXYZ);

// points above here
// ----------------------------------------------------------------------------------------------
// mesh from here
  

void CG_Tri(double& cgx, double& cgy, double& cgz,
            int itri,
            const std::vector<double>& aXYZ,
            const std::vector<int>& aTri);


/**
 * @brief center positions of each triangle and the maximum radius of the triangle
 * @details this funciton is implemented for "float" and double.
 * the aXYZ_c0 will be resized to aTri.size()/3
 */
template <typename T>
T CentsMaxRad_MeshTri3(std::vector<T>& aXYZ_c0,
                       const std::vector<T>& aXYZ,
                       const std::vector<unsigned int>& aTri);
  
template <typename T>
void CG_MeshTri3_Shell(T cg[3],
                       const std::vector<T>& aXYZ,
                       const std::vector<unsigned int>& aTri);
  
template <typename T>
T CG_TriMsh3Flg_Shell(T cg[3],
                      const std::vector<T>& aXYZ,
                      const std::vector<unsigned int>& aTri,
                      int iflg,
                      const std::vector<int>& aFlg);
  
template <typename T>
void CG_MeshTri3_Solid(T cg[3],
                       const std::vector<T>& aXYZ,
                       const std::vector<unsigned int>& aTri);
template <typename T>
void CG_MeshTet3(T& v_tot,
                 T cg[3],
                 const std::vector<T>& aXYZ,
                 const std::vector<unsigned int>& aTet);
  
// ---------------------------------
  
void RemoveUnreferencedPoints_MeshElem(std::vector<double>& aXYZ1,
                                       std::vector<unsigned int>& aElem1,
                                       std::vector<int>& aMap01,
                                       unsigned int ndim,
                                       const std::vector<double>& aXYZ0,
                                       const std::vector<unsigned int>& aElem0);

void Normal_MeshTri3D(double *aNorm,
                      const double *aXYZ, unsigned int nXYZ,
                      const unsigned int *aTri, unsigned int nTri);

void Quality_MeshTri2D(
    double &max_aspect, double &min_area,
    const double *aXY,
    const unsigned int *aTri, unsigned int nTri);

// ----------------------
// set primitive mesh

void SetTopology_ExtrudeTri2Tet(
    unsigned int *aTet,
    int nXY,
    const unsigned int *aTri, int nTri,
    int nlayer);

void ExtrudeTri2Tet(
    int nlayer, double h,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTet,
    const std::vector<double> &aXY,
    const std::vector<unsigned int> &aTri);

void makeSolidAngle(
    std::vector<double> &aSolidAngle,
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aNorm,
    std::vector<int> &elsup_ind,
    std::vector<int> &elsup);

/**
 * @brief Compute mass of the points (lumped mass) for 3D tet mesh
 * @details aMassMatrixLumped need to be allocated before in the size of nXY
 * this is here because both "fem" and "pbd" use this function
 */
void MassPoint_Tet3D(
    double *aMassMatrixLumped,
    double rho,
    const double *aXYZ, unsigned int nXYZ,
    const unsigned int *aTet, unsigned int nTet);

/**
 * @brief Compute mass of the points (lumped mass) for 2D triangle mesh
 * @param aMassMatrixLumped (out) this need to be allocated before in the size of nXY
 * @details this is here because both "fem" and "pbd" use this function
 */
void MassPoint_Tri2D(
    double *aMassMatrixLumped,
    //
    double rho,
    const double *aXY, unsigned int nXY,
    const unsigned int *aTri, unsigned int nTri);


void LaplacianSmoothing(std::vector<double>& aXYZ,
                        const std::vector<int>& aTri,
                        const std::vector<int>& elsup_ind,
                        const std::vector<int> elsup);

// ---------------------------------------------------------

void SubdivisionPoints_QuadCatmullClark(std::vector<double>& aXYZ1,
                                        //
                                        const std::vector<unsigned int>& aQuad1,
                                        const std::vector<int>& aEdgeFace0,
                                        const std::vector<unsigned int> &psupIndQuad0,
                                        const std::vector<unsigned int> &psupQuad0,
                                        const unsigned int* aQuad0, unsigned int nQuad0,
                                        const double* aXYZ0, unsigned int nXYZ0);
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
                           const std::vector<unsigned int> &psupIndHex0,
                           const std::vector<unsigned int> &psupHex0,
                           const std::vector<unsigned int>& aQuadHex0,
                           const unsigned int* aHex0, unsigned int nHex0,
                           const double* aXYZ0, unsigned int nXYZ0);

} // delfem2

#endif
