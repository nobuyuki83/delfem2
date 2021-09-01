/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

// (2020/12/26) DONE: change name to "mshprimitive.h" (2020/12/23)

#ifndef DFM2_MSHPRIMITIVE_H
#define DFM2_MSHPRIMITIVE_H

#include <math.h>
#include <vector>

#include "delfem2/isrf_sdf.h"
#include "delfem2/dfm2_inline.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846264338327950288
#endif

namespace delfem2 {

// -----------------------
// 3D primitives

/**
 * @brief a function to make a triangle mesh of sphere. The y axis is the pole.
 * @param radius  radius of  a sphere
 * @param nlong  number of subdivision in y-axis (longtitude)
 * @param nlat numbef of subdivision around the cross section of XZ plane (latitude)
 */
DFM2_INLINE void MeshTri3D_Sphere(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    double radius,
    int nlong, int nlat);

/**
 * @details y axis is the axis of cylinder
 */
DFM2_INLINE void MeshTri3D_CylinderOpen(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    double r, double l,
    int nr, int nl);

/**
 * @brief makae a mesh of 3D cylinder with closed ends
 * @details y axis is the axis of cylinder.
 * The first poit and the last points are at the center of the caps
 * @param nr number of division for circumference
 * @param nl number of division in axis direction
 */
template<typename T>
DFM2_INLINE void MeshTri3D_CylinderClosed(
    std::vector<T> &aXYZ,
    std::vector<unsigned int> &aTri,
    T r,
    T l,
    unsigned int nr,
    unsigned int nl);

template<typename T>
DFM2_INLINE void MeshTri3D_Cube(
    std::vector<T> &aXYZ,
    std::vector<unsigned int> &aTri,
    unsigned int n);

DFM2_INLINE void MeshTri3D_Disk(
    std::vector<double> &aXYZ, std::vector<unsigned int> &aTri,
    double r, int nr, int nth);

DFM2_INLINE void MeshTri3D_Icosahedron(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri);

// -----------------------
// below: cube

/**
 * @brief triangle mesh of torus.
 * @tparam REAL defined for "flaot" and "double" when DFM_HEADER_ONY is not defined.
 * @param len_latitude (in) length of the latitude.
 * @param len_merdician (in) length of the meridian. smaller length
 * @param ndiv_latitude (in) number of division in latitude
 * @param ndiv_meridian (in) number of division in meridian
 */
template<typename REAL>
DFM2_INLINE void MeshTri3_Torus(
    std::vector<REAL> &aXYZ,
    std::vector<unsigned int> &aTri,
    REAL len_latitutde,
    REAL len_meridian,
    unsigned int ndiv_latittude,
    unsigned int ndiv_meridian);

std::vector<unsigned int>
DFM2_INLINE MeshQuadTopo_CubeVox();

/**
 * @brief making a quad mesh of a cube. The order of the vertex is same as voxel element
 * @tparam REAL in the static compile, the REAL is defined for "float" and "double"
 */
template<typename REAL>
DFM2_INLINE void MeshQuad3_CubeVox(
    std::vector<REAL> &aXYZ,
    std::vector<unsigned int> &aQuad,
    const REAL bbmin[3],
    const REAL bbmax[3]);

// above: cube
// -----------
// below: capsule

template<typename T>
void MeshTri3_Capsule(
    std::vector<T> &aXYZ,
    std::vector<unsigned int> &aTri,
    T r,
    T l,
    unsigned int nc,
    unsigned int nr,
    unsigned int nl);

template<typename T>
DFM2_INLINE void MeshHex3_Grid(
    std::vector<T> &aXYZ,
    std::vector<unsigned int> &aHex,
    unsigned int nx,
    unsigned int ny,
    unsigned int nz,
    T elen);

// above: 3D primitives
// --------------------------------------------------------------------------------------------------------------
// below: 2D primitives

DFM2_INLINE void MeshQuad2D_Grid(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aQuad,
    unsigned int nx,
    unsigned int ny);

// =================================================================================
// primivive classes from here

template<typename REAL>
class CPlane : public CSDF3 {
 public:
  CPlane(const double norm[3], const double orgn[3]);
  double Projection(double n[3],
                    double px, double py, double pz) const; // normal
 public:
  double normal_[3] = {0, 0, 1};
  double origin_[3] = {0, 0, 0};
};

// ----------

template<typename REAL>
class CSphere : public CSDF3 {
 public:
  CSphere() {}
  CSphere(double rad, const std::vector<double> &c, bool is_out);
  // return penetration depth (inside is positive)
  virtual double Projection(double n[3],
                            double px, double py, double pz) const; // normal outward
  unsigned int FindInOut(double px, double py, double pz) const;
  bool IntersectionPoint(double p[3],
                         const double org[3], const double dir[3]) const;
 public:
  std::vector<double> cent_ = {0, 0, 0.0, 0.0};
  double radius_ = 1.0;
  bool is_out_ = true;    // true:normal points outward
};

// ------------

template<typename REAL>
class CCylinder : public CSDF3 {
 public:
  CCylinder() {};
  CCylinder(double r, const double cnt[3], const double dir[3], bool is_out);
  // return penetration depth (inside is positive)
  double Projection(double n[3],
                    double px, double py, double pz) const; // normal outward
  unsigned int FindInOut(double px, double py, double pz) const;
  bool IntersectionPoint(double p[3],
                         const double org[3], const double dir[3]) const;
  // ------
  void SetCenter(const double cnt[3]) {
    cent_[0] = cnt[0];
    cent_[1] = cnt[1];
    cent_[2] = cnt[2];
  }
  void SetDirection(const double dir[3]) {
    dir_[0] = dir[0];
    dir_[1] = dir[1];
    dir_[2] = dir[2];
  }
  void SetRadius(double r) { radius_ = r; }
 public:
  double cent_[3]{0.0, 0.0, 0.0};
  double dir_[3] = {1.0, 0.0, 0.0};
  double radius_ = 1.0;
  bool is_out_ = true;    // true:normal points outward
};

template<typename REAL>
class CTorus : public CSDF3 {
 public:
  CTorus() {}
  // return penetration depth (inside is positive)
  virtual double Projection(double n[3],
                            double px, double py, double pz) const; // normal outward
  unsigned int FindInOut(double px, double py, double pz) const;
 public:
  double cent_[3] = {0.0, 0.0, 0.0};
  double radius_ = 0.5;
  double radius_tube_ = 0.2;
};

template<typename REAL>
class CBox : public CSDF3 {
 public:
  virtual double Projection(double n[3],
                            double x, double y, double z) const {
    double len0, x0, y0, z0;
    {
      x0 = (x < 0) ? -hwx : +hwx;
      y0 = (fabs(y) < hwy) ? y : ((y < 0) ? -hwy : +hwy);
      z0 = (fabs(z) < hwz) ? z : ((z < 0) ? -hwz : +hwz);
      len0 = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0) + (z - z0) * (z - z0));
    }
    double len1, x1, y1, z1;
    {
      x1 = (fabs(x) < hwx) ? x : ((x < 0) ? -hwx : +hwx);
      y1 = (y < 0) ? -hwy : +hwy;
      z1 = (fabs(z) < hwz) ? z : ((z < 0) ? -hwz : +hwz);
      len1 = sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1));
    }
    double len2, x2, y2, z2;
    {
      x2 = (fabs(x) < hwx) ? x : ((x < 0) ? -hwx : +hwx);
      y2 = (fabs(y) < hwy) ? y : ((y < 0) ? -hwy : +hwy);
      z2 = (z < 0) ? -hwz : +hwz;
      len2 = sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2) + (z - z2) * (z - z2));
    }
    double len3, x3, y3, z3;
    if (len0 < len1) {
      len3 = len0;
      x3 = x0;
      y3 = y0;
      z3 = z0;
    }
    else {
      len3 = len1;
      x3 = x1;
      y3 = y1;
      z3 = z1;
    }
    if (len3 > len2) {
      len3 = len2;
      x3 = x2;
      y3 = y2;
      z3 = z2;
    }
    {
      n[0] = x - x3;
      n[1] = y - y3;
      n[2] = z - z3;
      double tmp = 1.0 / sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
      n[0] *= tmp;
      n[1] *= tmp;
      n[2] *= tmp;
    }
    if (fabs(x) > hwx || fabs(y) > hwy || fabs(z) > hwz) { // outside
      len3 = -len3;
      n[0] *= -1;
      n[1] *= -1;
      n[2] *= -1;
    }
    return len3;
  }
  bool IntersectionPoint(double p[3],
                         const double org[3], const double dir[3]) const { return true; }

 public:
  double hwx = 1; // half x width
  double hwy = 1; // half y width
  double hwz = 1; // half z width
};

/*
class CSDF3_Combine : public CSDF3
{
public:
	 CSDF3_Combine(){}
	~CSDF3_Combine(){ 
    for(unsigned int ipct=0;ipct<apCT.size();ipct++){ delete apCT[ipct]; }
  }
	virtual double Projection
	(double px, double py, double pz,
	 double n[3]) const; // normal
  virtual bool IntersectionPoint
  (double p[3], 
   const double org[3], const double dir[3]) const { return  true; }
private:
  std::vector<CSDF3*> apCT;    
};
 */

/*
class CSDF3_Transform : public CSDF3
{
public:
	CSDF3_Transform(CSDF3* pCT){ 
		phi = 0; theta = 0; psi = 0;		
		trans[0]=0; trans[1]=0;	trans[2]=0;
		this->pCT = pCT; 
	}
	~CSDF3_Transform(){ delete pCT; }
	void SetAngle_Bryant(double phi, double theta, double psi){
		this->phi = phi*3.1416/180.0;
		this->theta = theta*3.1416/180.0;
		this->psi = psi*3.1416/180.0;		
	}
	void SetTranslate(double x, double y, double z){
		trans[0] = x; trans[1] = y; trans[2] = z;		
	}
  
	virtual double Projection
	(double px, double py, double pz,
	 double n[3]) const; // normal
  virtual bool IntersectionPoint
  (double p[3], 
   const double org[3], const double dir[3]) const { return  true; }
private:
	double phi, theta, psi;	// Bryant Angle
	double trans[3];
	CSDF3* pCT;
};
 */


/*
////////////////////////////////////////////////////////////////////////////////
class CSpatialHash_Grid3D;
class CSDF3_Mesh : public CSDF3
{
public:
  CSDF3_Mesh();
  ~CSDF3_Mesh();
  // return penetration depth (inside is positive)
  virtual double Projection
  (double px, double py, double pz,
   double n[3]) const; // normal outward
  virtual bool IntersectionPoint(double p[3], const double org[3], const double dir[3]) const;
  
  ////////
  bool FindIntersectionTri(double psec[3], int& itri, double& r0, double& r1, const double org[3], const double dir[3]) const;
  double FindNearest(int& itri, double& r0, double& r1, double px, double py, double pz) const;
  void GetCenterWidth(double& cx, double& cy, double& cz,  double& wx, double& wy, double& wz);
  void Translate(double x, double y, double z);
  void BuildBoxel();
  void SetHole(bool is_hole){	this->is_hole = is_hole; }
  virtual void GetMesh(std::vector<unsigned int>& aTri, std::vector<double>& aXYZ, double elen) const;
  void SetMesh(const std::vector<unsigned int>& aTri, const std::vector<double>& aXYZ);
private:
  unsigned int FindInOut_IntersectionRay
  (double px, double py, double pz,
   const double dir[3]) const;
  unsigned int FindInOut_IntersectionRay_Boxel
  (double px, double py, double pz,
   const double dir[3]) const;
  virtual unsigned int FindInOut(double px, double py, double pz) const;
  virtual unsigned int FindInOut_Boxel
  (double px, double py, double pz) const;
  
public:
  bool is_hole;
  unsigned int nnode_;
  double* pXYZs_;
  unsigned int ntri_;
  unsigned int* aTri_;
  CSpatialHash_Grid3D* pBoxel_;
  mutable std::vector<unsigned int> aFlgTriUsed;
  mutable std::vector<unsigned int> aIndTriCand;
};
 */

// --------------------------------------------------------------------------------------------------------------------

} // end of delfem2


#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/mshprimitive.cpp"
#endif

#endif
