#ifndef DFM2_BEM
#define DFM2_BEM

#include <complex>
#include <vector>
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"

typedef std::complex<double> COMPLEX;
const COMPLEX IMG = COMPLEX(0.0, 1.0);

// -----------------------------------

bool Solve_BiCGSTAB(double& conv_ratio, int& iteration,
                    std::vector<double>& u_vec,
                    const std::vector<double>& A,
                    const std::vector<double>& y_vec);

// {y} = [A]{x}
void matVec(std::vector<double>& y,
            const std::vector<double>& A,
            const std::vector<double>& x);

double squaredNorm(const std::vector<double>& v);

// -----------------------------


inline delfem2::CVec3d MidPoint(
    int itri,
    const std::vector<unsigned int>& aTri,
    const std::vector<double>& aXYZ)
{
  delfem2::CVec3d p;
  int i0 = aTri[itri*3+0];
  int i1 = aTri[itri*3+1];
  int i2 = aTri[itri*3+2];
  p.p[0] = (aXYZ[i0*3+0]+aXYZ[i1*3+0]+aXYZ[i2*3+0])/3.0;
  p.p[1] = (aXYZ[i0*3+1]+aXYZ[i1*3+1]+aXYZ[i2*3+1])/3.0;
  p.p[2] = (aXYZ[i0*3+2]+aXYZ[i1*3+2]+aXYZ[i2*3+2])/3.0;
  return p;
}

// ----------------------------------

void makeLinearSystem_PotentialFlow_Order0th(std::vector<double>& A,
                                             std::vector<double>& f,
                                             //
                                             const delfem2::CVec3d& velo_inf,
                                             int ngauss,
                                             const std::vector<double>& aXYZ,
                                             const std::vector<unsigned int> &aTri);

//! @brief evaluate BEM solution where the value is constant over a triangle
void evaluateField_PotentialFlow_Order0th(double& phi_pos,
                                          delfem2::CVec3d& gradphi_pos,
                                          //
                                          const delfem2::CVec3d& pos,
                                          const delfem2::CVec3d& velo_inf,
                                          int ngauss,
                                          const std::vector<double>& aValTri,
                                          const std::vector<double>& aXYZ,
                                          const std::vector<unsigned int> &aTri);

// --------------------------------------

delfem2::CVec3d evaluateField_PotentialFlow_Order1st(double& phi_pos,
                                              const delfem2::CVec3d& pos,
                                              const delfem2::CVec3d& velo_inf,
                                              int ngauss,
                                              const std::vector<double>& aValSrf,
                                              const std::vector<double>& aXYZ,
                                              const std::vector<int>& aTri);

void makeLinearSystem_PotentialFlow_Order1st(std::vector<double>& A,
                                             std::vector<double>& f,
                                             //
                                             const delfem2::CVec3d& velo_inf,
                                             int ngauss,
                                             const std::vector<double>& aXYZ,
                                             const std::vector<int>& aTri,
                                             const std::vector<double>& aSolidAngle);

// -----------------------------------

void makeLinearSystem_VortexSheet_Order0th(std::vector<double>& A,
                                           std::vector<double>& f,
                                           //
                                           const delfem2::CVec3d& velo,
                                           int ngauss,
                                           const std::vector<double>& aXYZ,
                                           const std::vector<int>& aTri);

delfem2::CVec3d evaluateField_VortexSheet_Order0th(const delfem2::CVec3d& pos,
                                            const std::vector<double>& aValSrf,
                                            //
                                            int ngauss,
                                            const std::vector<double>& aXYZ,
                                            const std::vector<int>& aTri,
                                            int jtri_exclude);

class CVortexParticle
{
public:
  delfem2::CVec3d pos; // position
  delfem2::CVec3d circ; // circulation
  delfem2::CVec3d circ0;
  double rad; // radius
  delfem2::CMat3d m;
  //
  delfem2::CVec3d velo;
  delfem2::CVec3d velo_pre;
  delfem2::CMat3d gradvelo;
  delfem2::CMat3d gradvelo_pre;
public:
  void stepTime_AdamsBashforth(double dt){
//           pos += velo*dt;
    //      circ += gradvelo*circ*dt;
    pos += (1.5*velo-0.5*velo_pre)*dt;
    //    circ += (1.5*gradvelo*circ-0.5*gradvelo_pre*circ)*dt;
  }
};

delfem2::CVec3d veloVortexParticles(const delfem2::CVec3d& p0,
                             const std::vector<CVortexParticle>& aVortexParticle,
                             int ivp_self);

delfem2::CMat3d gradveloVortexParticles(delfem2::CVec3d& velo,
                                          const delfem2::CVec3d& p0,
                                          const std::vector<CVortexParticle>& aVortexParticle,
                                          int ivp_self);

void setGradVeloVortexParticles(std::vector<CVortexParticle>& aVortexParticle);

class CGrid_Vortex{
public:
  class CDataVtx{
  public:
    delfem2::CVec3d circ;
    std::vector< std::pair<int, double> > aPairPtcleWeight;
  };
public:
  int nx;
  int ny;
  int nz;
  double h;
  delfem2::CVec3d cnt;
  std::vector<CDataVtx> aDataVtx; // (nx+1)*(ny+1)*(nz+1)
public:
  void drawBoundingBox() const;
};

void viscousityVortexParticleGrid(std::vector<CVortexParticle>& aVortexParticle,
                                  CGrid_Vortex& grid, double resolution);



// --------------------------------------

COMPLEX evaluateField_Helmholtz_Order0th(const std::vector<COMPLEX>& aSol,
                                         const delfem2::CVec3d& p,
                                         const delfem2::CVec3d& pos_source,
                                         double k, // wave number
                                         double Adm, // admittance
                                         const std::vector<unsigned int> &aTri,
                                         const std::vector<double>& aXYZ,
                                         bool is_inverted_norm);

COMPLEX evaluateField_Helmholtz_Order1st(const std::vector<COMPLEX>& aSol,
                                         const delfem2::CVec3d& p,
                                         const delfem2::CVec3d& pos_source,
                                         double k, // wave number
                                         double beta, // admittance
                                         const std::vector<int>& aTri,
                                         const std::vector<double>& aXYZ,
                                         bool is_inverted_norm,
                                         int ngauss);


#endif
