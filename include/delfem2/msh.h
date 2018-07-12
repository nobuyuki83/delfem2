#if !defined(SURFACE_MESH_READER_H)
#define SURFACE_MESH_READER_H

#include <string>
#include <vector>
#include <fstream>

////////////////////////////////////////////////////////////////////////
// work on points

void GetCenterWidth(double& cx, double& cy, double& cz,
                    double& wx, double& wy, double& wz,
                    const int nXYZ, const double* paXYZ);
void GetCenterWidth(double& cx, double& cy, double& cz,
                    double& wx, double& wy, double& wz,
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
// Bryant angle
void Rotate(std::vector<double>& aXYZ,
            double radx, double rady, double radz);
void CenterOfGravity(double& cgx, double& cgy, double& cgz,
                     const std::vector<double>& aXYZ);
void Translate(double tx, double ty, double tz,
               std::vector<double>& aXYZ);
void Translate(std::vector<double>& aXYZ,
               double tx, double ty, double tz);
void Scale(double s,
           std::vector<double>& aXYZ);
double Size(const std::vector<double>& aXYZ);
void Normalize(std::vector<double>& aXYZ,
               double s = 1.0);

//////////////////////////////////

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
void RemoveUnreferencedPoints3D(std::vector<double>& aXYZOut,
                                std::vector<int>& aTriOut,
                                const std::vector<double>& aXYZIn,
                                const std::vector<int>& aTriIn);
void MakeNormal(std::vector<double>& aNorm,
                const std::vector<double>& aXYZ,
                const std::vector<int>& aTri);

////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// set primitive mesh

// y axis is the pole
void setSphereMesh(std::vector<double>& aXYZ, std::vector<int>& aTri,
                   double r,
                   int nla, int nlo);
// y axis is the axis of cylinder
void setOpenCylinderMesh(std::vector<double>& aXYZ, std::vector<int>& aTri,
                         double r, double l,
                         int nr, int nl);
void setClosedCylinderMesh(std::vector<double>& aXYZ, std::vector<int>& aTri,
                           double r, double l,
                           int nr, int nl);
void setQuad_CubeTopo(std::vector<int>& aQuad);
void setCubeMesh(std::vector<double>& aXYZ, std::vector<int>& aTri,
                 int n);
void setDiskMesh(std::vector<double>& aXYZ, std::vector<int>& aTri,
                 double r, int nr, int nth);
void ExtrudeTri2Tet(int nlayer, double h,
                    std::vector<double>& aXYZ,
                    std::vector<int>& aTet,
                    const std::vector<double>& aXY,
                    const std::vector<int>& aTri);


///////////////////
void makeSolidAngle(std::vector<double>& aSolidAngle,
                    const std::vector<double>& aXYZ,
                    const std::vector<int>& aTri);
void LaplacianSmoothing(std::vector<double>& aXYZ,
                        const std::vector<int>& aTri,
                        const std::vector<int>& elsup_ind,
                        const std::vector<int> elsup);
void makeSolidAngle(std::vector<double>& aSolidAngle,
                    const std::vector<double>& aXYZ,
                    const std::vector<int>& aTri,
                    const std::vector<double>& aNorm,
                    std::vector<int>& elsup_ind,
                    std::vector<int>& elsup);
void SubdivisionMovePoints_QuadCatmullClark(std::vector<double>& aXYZ1,
                                            ///
                                            const std::vector<int>& aQuad1,
                                            const std::vector<int>& aEdgeFace0,
                                            const std::vector<int>& psupIndQuad0,
                                            const std::vector<int>& psupQuad0,
                                            const std::vector<int>& aQuad0,
                                            const std::vector<double>& aXYZ0);

//////////////////////////////////////////////////////////////
// raw mesh functions

void MakeNormal
(double*& aNorm_,
 const unsigned int nnode_, const double* pXYZs_,
 const unsigned int ntri_,  const unsigned int* aTri_);

void Translate
(double tx, double ty, double tz,
 const unsigned int nnode_, double* pXYZs_);

void Scale
(double s,
 const unsigned int nnode_, double* pXYZs_);

#endif
