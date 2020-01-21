/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_VEC3_H
#define DFM2_VEC3_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#define NEARLY_ZERO 1.e-16


namespace delfem2 {

template <typename T>
T Distance3(const T p0[3], const T p1[3]);

template <typename T>
T SquareDistance3(const T p0[3], const T p1[3]);

template <typename T>
T Length3(const T v[3]);

template <typename T>
T SquareLength3(const T v[3]);
  
template <typename T>
T Dot3(const T a[3], const T b[3]);
  
template <typename T>
T Volume_Tet3(const T v1[3], const T v2[3], const T v3[3], const T v4[3]);
  
template <typename T>
void Cross3(T r[3],
            const T v1[3], const T v2[3]);

template <typename T>
void Normalize3(T v[3]);

template <typename T>
T Area_Tri3(const T v1[3], const T v2[3], const T v3[3]);
  
template <typename T>
void MatVec3(T y[3],
             const T m[9], const T x[3]);


double ScalarTripleProduct3D(const double a[], const double b[], const double c[]);

void UnitNormalAreaTri3D(double n[3], double& a,
                         const double v1[3], const double v2[3], const double v3[3]);
void NormalTri3D(double n[3],
                 const double v1[3], const double v2[3], const double v3[3]);
void InverseMat3(double Ainv[],
                 const double A[]);
void transposeMat3(double At[],
                   const double A[]);
void GetVertical2Vector3D(const double vec_n[3], double vec_x[3], double vec_y[3]);
void MatTransVec3(double y[3],
                  const double m[9], const double x[3]);
void VecMat3(double y[3],
             const double x[3], const double m[9]);
void GetRotMatrix_Rodrigues3D(double rot[9],
                              const double n[3], double theta);
void Mat4Vec3(double vo[3],
              const double M[16], const double vi[3]);

// -------------------------------------------------------------

class CVector3;
double Dot(const CVector3 &arg1, const CVector3 &arg2);
CVector3 Cross(const CVector3& arg1, const CVector3& arg2);
CVector3 operator + (const CVector3& lhs, const CVector3& rhs);
CVector3 operator - (const CVector3& lhs, const CVector3& rhs);
CVector3 operator * (double d, const CVector3& rhs);
CVector3 operator * (const CVector3& vec, double d);
double operator * (const CVector3& lhs, const CVector3& rhs);
CVector3 operator / (const CVector3& vec, double d);
CVector3 operator ^ (const CVector3& lhs, const CVector3& rhs);
std::ostream &operator<<(std::ostream &output, const CVector3& v);
std::istream &operator>>(std::istream &input, CVector3& v);

/**
 * @brief 3 dimentional vector class
 */
class CVector3  
{
public:
  CVector3(double vx, double vy, double vz) : p{vx,vy,vz} {}
  CVector3(): p{0.0, 0.0, 0.0} {}
	CVector3(const CVector3& rhs){ p[0] = rhs.p[0]; p[1] = rhs.p[1]; p[2] = rhs.p[2]; }
  CVector3(const double* prhs){ p[0] = prhs[0]; p[1] = prhs[1]; p[2] = prhs[2]; }
  CVector3(const float* prhs){ p[0] = prhs[0]; p[1] = prhs[1]; p[2] = prhs[2]; }
  CVector3(const std::vector<double>& v){ p[0] = v[0]; p[1] = v[1]; p[2] = v[2]; }
	virtual ~CVector3(){}

  std::vector<double> stlvec() const {
    std::vector<double> d = {p[0], p[1], p[2]};
    return d;
  }
	void SetVector(double vx, double vy, double vz){ p[0] = vx; p[1] = vy; p[2] = vz; }
  void CopyValueTo(double* v) const { v[0]=p[0]; v[1]=p[1]; v[2]=p[2]; }
  void CopyValueToScale(double* v, double s) const { v[0]=p[0]*s; v[1]=p[1]*s; v[2]=p[2]*s; }

	inline const CVector3 operator-() const{ return -1.0*(*this); }
	inline const CVector3 operator+() const{ return *this; }  
	inline CVector3& operator=(const CVector3& rhs){
		if( this != &rhs ){ p[0]= rhs.p[0]; p[1] = rhs.p[1]; p[2] = rhs.p[2]; }
		return *this;
	}
	inline CVector3& operator+=(const CVector3& rhs){
		p[0] += rhs.p[0];
    p[1] += rhs.p[1];
    p[2] += rhs.p[2];
		return *this;
	}
	inline CVector3& operator-=(const CVector3& rhs){
		p[0] -= rhs.p[0];
    p[1] -= rhs.p[1];
    p[2] -= rhs.p[2];
		return *this;
	}
	inline CVector3& operator*=(double d){
		p[0] *= d;
    p[1] *= d;
    p[2] *= d;
		return *this;
	}
	inline CVector3& operator/=(double d){
		if( fabs(d) < NEARLY_ZERO ){ return *this; }
		p[0] /= d;
    p[1] /= d;
    p[2] /= d;
		return *this;
	}
  inline double operator[](int i) const{
    if( i == 0 ) return p[0];
    if( i == 1 ) return p[1];
    if( i == 2 ) return p[2];
    return 0;
  }
  inline double& operator[](int i){
    if( i == 0 ) return p[0];
    if( i == 1 ) return p[1];
    if( i == 2 ) return p[2];
    assert(0);
    return p[0];
  }  
	inline CVector3 operator+(){ return *this; }
	inline CVector3 operator-(){ return CVector3(-p[0],-p[1],-p[2]); }
  inline double x() const { return p[0]; }
  inline double y() const { return p[1]; }
  inline double z() const { return p[2]; }

	friend bool operator==(const CVector3&, const CVector3&);
	friend bool operator!=(const CVector3&, const CVector3&);

	friend CVector3 Cross(const CVector3&, const CVector3&);
	friend double Dot(const CVector3&, const CVector3&);
  CVector3 Normalize() const {
    CVector3 r = (*this);
    r.SetNormalizedVector();
    return r;
  }
	inline double Length()  const{ return sqrt( p[0]*p[0]+p[1]*p[1]+p[2]*p[2] ); }
	inline double DLength() const{ return p[0]*p[0]+p[1]*p[1]+p[2]*p[2]; }
	void SetNormalizedVector();
	void SetZero();
  void Print() const {
    std::cout <<p[0]<< " " << p[1] << " " << p[2] << std::endl;
  }
  bool isNaN() const{
    double s=p[0]+p[1]+p[2];
    return !(s > s-1.0);
  }
  static CVector3 Axis(int idim){
    CVector3 r(0,0,0);
    if( idim >= 0 && idim < 3) { r[idim] = 1; }
    return r;
  }
public:
  double p[3];
  /*
	double x;	//!< x axis coordinate
	double y;	//!< y axis coordinate
	double z;	//!< z axis coordinate
   */
};


// --------------------------------------------------------------------------
// rule about naming, the method starts "Set" change it self (not const)
  
CVector3 mult_GlAffineMatrix(const float* m,
                             const CVector3& p);
CVector3 solve_GlAffineMatrix(const float* m,
                              const CVector3& p);
CVector3 solve_GlAffineMatrixDirection(const float* m,
                                       const CVector3& v);

CVector3 Mat3Vec(const double M[ 9], const CVector3& v);
CVector3 Mat4Vec(const double M[16], const CVector3& v);
CVector3 QuatVec(const double quat[4], const CVector3& v0);
CVector3 QuatConjVec(const double quat[4], const CVector3& v0);
CVector3 screenProjection(const CVector3& v,
                          const float* mMV, const float* mPj);
// opposite to the screen normal
CVector3 screenDepthDirection(const CVector3& v,
                              const float* mMV,
                              const float* mPj);
CVector3 screenUnProjection(const CVector3& v,
                            const float* mMV,
                            const float* mPj);
CVector3 screenUnProjectionDirection(const CVector3& v,
                                     const float* mMV,
                                     const float* mPj);
CVector3 screenDepthDirection(const CVector3& v,
                              const float* mMV,
                              const float* mPj);

double ScalarTripleProduct(const CVector3& a, const CVector3& b, const CVector3& c);
bool operator == (const CVector3& lhs, const CVector3& rhs);
bool operator != (const CVector3& lhs, const CVector3& rhs);
double Height(const CVector3& v1, const CVector3& v2, const CVector3& v3, const CVector3& v4);
void GetVertical2Vector (const CVector3& vec_n, CVector3& vec_x, CVector3& vec_y);

// --------------------------------------------------------------

CVector3 nearest_Line_Point(const CVector3& p, // point
                              const CVector3& s, // source
                              const CVector3& d);
CVector3 nearest_Line_Point(double& t,
                              const CVector3& p, // point
                              const CVector3& s, // source
                              const CVector3& d); // direction
CVector3 nearest_LineSeg_Point(const CVector3& p, // point
                                 const CVector3& s, // source
                                 const CVector3& e); // end
CVector3 nearest_LineSeg_Point(double& t,
                                 const CVector3& p, // point
                                 const CVector3& s,
                                 const CVector3& e); // direction
CVector3 nearest_Plane_Point(const CVector3& p, // point
                             const CVector3& o, // origin
                             const CVector3& n); // normal
CVector3 Nearest_Origin_Tri(double& r0,
                            double& r1,
                            const CVector3& q0,
                            const CVector3& q1,
                            const CVector3& q2);
CVector3 nearst_Origin_Quad(double& s0, double& s1,
                            const CVector3& p,
                            const CVector3& q0,
                            const CVector3& q1,
                            const CVector3& q2,
                            const CVector3& q3);

void nearest_Line_Line(double& D, CVector3& Da, CVector3& Db,
                         const CVector3& pa_, const CVector3& va,
                         const CVector3& pb_, const CVector3& vb);
void nearest_Line_Line(double& D, CVector3& Da, CVector3& Db,
                         double& ta, double& tb,
                         const CVector3& pa_, const CVector3& va,
                         const CVector3& pb_, const CVector3& vb);
void Nearest_Line_Circle(CVector3& p0,
                         CVector3& q0,
                         const CVector3& src,
                         const CVector3& dir,
                         const CVector3& org,
                         const CVector3& normal,
                         double rad);
  
CVector3 nearst_Origin_Quad(double& s0,
                            double& s1,
                            const CVector3& q0,
                            const CVector3& q1,
                            const CVector3& q2,
                            const CVector3& q3);
  
CVector3 nearest_Origin_LineSeg(const CVector3& s, // start
                                const CVector3& e); // end

// r0==0 -> p0==org
// r0==1 -> p1==org
CVector3 nearest_Origin_LineSeg(double& r0,
                                const CVector3& p0, // start
                                const CVector3& p1); // end

  
CVector3 Nearest_Orgin_PlaneTri(double& r0,
                                double& r1,
                                const CVector3& q0,
                                const CVector3& q1,
                                const CVector3& q2);
    
// -----------------------------------------

bool intersection_Plane_Line(CVector3& p0, double& r0, double& r1, double& r2,
                             double eps,
                             const CVector3& src, const CVector3& dir,
                             const CVector3& q0, const CVector3& q1, const CVector3& q2);
bool intersection_Point_Quad(CVector3& psec, double& s0, double& s1,
                             const CVector3& src, const CVector3& dir,
                             const CVector3& q0, const CVector3& q1, const CVector3& q2, const CVector3& q3);
CVector3 intersection_Plane_Line(const CVector3& o, // one point on plane
                                 const CVector3& n, // plane normal
                                 const CVector3& s, // one point on line
                                 const CVector3& d); // direction of line

// ----------------------------------------------------------------------------------

double DistanceFaceVertex(const CVector3& p0, const CVector3& p1, const CVector3& p2,
                          const CVector3& p3,
                          double& w0, double& w1);

double DistanceEdgeEdge(const CVector3& p0, const CVector3& p1,
                        const CVector3& q0, const CVector3& q1,
                        double& ratio_p, double& ratio_q);

double FindCoplanerInterp(const CVector3& s0, const CVector3& s1, const CVector3& s2, const CVector3& s3,
                          const CVector3& e0, const CVector3& e1, const CVector3& e2, const CVector3& e3);

bool IsContact_EE_Proximity(int ino0,        int ino1,        int jno0,        int jno1,
                            const CVector3& p0, const CVector3& p1, const CVector3& q0, const CVector3& q1,
                            const double delta);

bool IsContact_FV_CCD2(int ino0,        int ino1,        int ino2,        int ino3,
                       const CVector3& p0, const CVector3& p1, const CVector3& p2, const CVector3& p3,
                       const CVector3& q0, const CVector3& q1, const CVector3& q2, const CVector3& q3);

bool isIntersectTriPair(CVector3& P0, CVector3& P1,
                        int itri, int jtri,
                        const std::vector<unsigned int>& aTri,
                        const std::vector<double>& aXYZ);

// ---------------------------------------------------------

bool barycentricCoord_Origin_Tet(double& r0,
                                 double& r1,
                                 double& r2,
                                 const CVector3& p0,
                                 const CVector3& p1,
                                 const CVector3& p2,
                                 const CVector3& p3);
bool barycentricCoord_Origin_Pyramid(double& r0,
                                     double& r1,
                                     double& r2,
                                     const CVector3& p0,
                                     const CVector3& p1,
                                     const CVector3& p2,
                                     const CVector3& p3,
                                     const CVector3& p4);
bool barycentricCoord_Origin_Wedge(double& r0,
                                   double& r1,
                                   double& r2,
                                   const CVector3& p0,
                                   const CVector3& p1,
                                   const CVector3& p2,
                                   const CVector3& p3,
                                   const CVector3& p4,
                                   const CVector3& p5);
  
CVector3 positionBarycentricCoord_Pyramid(double r0,
                                          double r1,
                                          double r2,
                                          const CVector3& p0,
                                          const CVector3& p1,
                                          const CVector3& p2,
                                          const CVector3& p3,
                                          const CVector3& p4);
  
CVector3 positionBarycentricCoord_Wedge(double r0,
                                        double r1,
                                        double r2,
                                        const CVector3& p0,
                                        const CVector3& p1,
                                        const CVector3& p2,
                                        const CVector3& p3,
                                        const CVector3& p4,
                                        const CVector3& p5);
  
void iteration_barycentricCoord_Origin_Solid(double& r0,
                                             double& r1,
                                             double& r2,
                                             const CVector3& q, // q=positionBarycentricCoord_Wedge
                                             const CVector3& dpdr0,
                                             const CVector3& dpdr1,
                                             const CVector3& dpdr2,
                                             double damp=1.0);

  
// ---------------------------------------------------------

bool IsInside_Orgin_BoundingBoxPoint4(const CVector3& p0,
                                    const CVector3& p1,
                                    const CVector3& p2,
                                    const CVector3& p3);

bool IsInside_Orgin_BoundingBoxPoint5(const CVector3& p0,
                                    const CVector3& p1,
                                    const CVector3& p2,
                                    const CVector3& p3,
                                    const CVector3& p4);

bool IsInside_Orgin_BoundingBoxPoint6(const CVector3& p0,
                                    const CVector3& p1,
                                    const CVector3& p2,
                                    const CVector3& p3,
                                    const CVector3& p4,
                                    const CVector3& p5);

// -----------------------------------------------

double volume_OrgTet(const CVector3& v1,
                     const CVector3& v2,
                     const CVector3& v3 );
double Volume_Tet(const CVector3& v0,
                  const CVector3& v1,
                  const CVector3& v2,
                  const CVector3& v3 );
double volume_Pyramid(const CVector3& p0,
                      const CVector3& p1,
                      const CVector3& p2,
                      const CVector3& p3,
                      const CVector3& p4);
double volume_Wedge(const CVector3& p0,
                    const CVector3& p1,
                    const CVector3& p2,
                    const CVector3& p3,
                    const CVector3& p4,
                    const CVector3& p5);

// ---------------------------------------------

double SolidAngleTri(const CVector3& v1,
                     const CVector3& v2,
                     const CVector3& v3);

void Cross( CVector3& lhs, const CVector3& v1, const CVector3& v2 );
double Area_Tri(const CVector3& v1, const CVector3& v2, const CVector3& v3);
double SquareTriArea(const CVector3& v1, const CVector3& v2, const CVector3& v3);
double SquareDistance(const CVector3& ipo0, const CVector3& ipo1);
double SquareLength(const CVector3& point);
double Length(const CVector3& point);
double Distance(const CVector3& ipo0, const CVector3& ipo1);

double SqareLongestEdgeLength(const CVector3& ipo0,
                              const CVector3& ipo1,
                              const CVector3& ipo2,
                              const CVector3& ipo3 );

double LongestEdgeLength(const CVector3& ipo0,
                         const CVector3& ipo1,
                         const CVector3& ipo2,
                         const CVector3& ipo3 );

double SqareShortestEdgeLength(const CVector3& ipo0,
                               const CVector3& ipo1,
                               const CVector3& ipo2,
                               const CVector3& ipo3 );

double ShortestEdgeLength(const CVector3& ipo0,
                          const CVector3& ipo1,
                          const CVector3& ipo2,
                          const CVector3& ipo3 );

void Normal(CVector3& vnorm,
            const CVector3& v1,
            const CVector3& v2,
            const CVector3& v3);

CVector3 Normal(const CVector3& v1,
                const CVector3& v2,
                const CVector3& v3);

void UnitNormal(CVector3& vnorm,
                const CVector3& v1,
                const CVector3& v2,
                const CVector3& v3);
  
CVector3 UnitNormal(const CVector3& v1,
                    const CVector3& v2,
                    const CVector3& v3);
  
/**
 * @function check if Delaunay condition satisfied
 * @return
 * 0 : p3 is inside circum circle on the p0,p1,p2
 * 1 :       on
 * 2 :       outsdie
 */
int DetDelaunay(const CVector3& p0,
                const CVector3& p1,
                const CVector3& p2,
                const CVector3& p3);

double SquareCircumradius(const CVector3& ipo0,
                          const CVector3& ipo1,
                          const CVector3& ipo2,
                          const CVector3& ipo3);

CVector3 CircumCenter(const CVector3& ipo0,
                      const CVector3& ipo1,
                      const CVector3& ipo2,
                      const CVector3& ipo3);

double Circumradius(const CVector3& ipo0,
                    const CVector3& ipo1,
                    const CVector3& ipo2,
                    const CVector3& ipo3);
CVector3 RotateVector(const CVector3& vec0, const CVector3& rot );
CVector3 RandVector();
CVector3 RandUnitVector();
CVector3 RandGaussVector();

void MeanValueCoordinate(double w[3],
                         const CVector3& v0,
                         const CVector3& v1,
                         const CVector3& v2);

CVector3 ProjectPointOnTriangle(const CVector3 &p0,
                                const CVector3 &tri_p1, const CVector3 &tri_p2, const CVector3 &tri_p3);
bool isPointInsideTriangle(const CVector3 &p0,
                           const CVector3 &tri_p1, const CVector3 &tri_p2, const CVector3 &tri_p3);
bool isPointSameSide(const CVector3 &p0, const CVector3 &p1,
                     const CVector3 &line_p0, const CVector3 &line_p1);
bool isRayIntersectingTriangle(const CVector3 &line0, const CVector3 &line1,
                               const CVector3 &tri0, const CVector3 &tri1, const CVector3 &tri2,
                               CVector3 &intersectionPoint);

void GetConstConstDiff_Bend(double& C,
                            CVector3 dC[4],
                            // -----
                            const CVector3& p0,
                            const CVector3& p1,
                            const CVector3& p2,
                            const CVector3& p3);
  
void CheckConstDiff_Bend();
  
// ----------------------------------------------------------
// here starts std::vector<CVector3>

double TetVolume( int iv1, int iv2, int iv3, int iv4,
                 const std::vector<CVector3>& node);

double volume_Tet(int iv1, int iv2, int iv3, int iv4,
                  const std::vector<CVector3>& aPoint);

double TriArea(const int iv1, const int iv2, const int iv3,
               const std::vector<CVector3>& node );

bool IsOut(int itri, const CVector3& v,
           const std::vector<CVector3>& aXYZ,
           const std::vector<int>& aTri);
  
void ConvexHull(std::vector<int>& aTri,
                const std::vector<CVector3>& aXYZ);

inline CVector3 cg_Tri(unsigned int itri,
                       const std::vector<unsigned int>& aTri,
                       const std::vector<double>& aXYZ)
{
  CVector3 p;
  int i0 = aTri[itri*3+0];
  int i1 = aTri[itri*3+1];
  int i2 = aTri[itri*3+2];
  p.p[0] = (aXYZ[i0*3+0]+aXYZ[i1*3+0]+aXYZ[i2*3+0])/3.0;
  p.p[1] = (aXYZ[i0*3+1]+aXYZ[i1*3+1]+aXYZ[i2*3+1])/3.0;
  p.p[2] = (aXYZ[i0*3+2]+aXYZ[i1*3+2]+aXYZ[i2*3+2])/3.0;
  return p;
}

inline CVector3 normalTri(int itri,
                          const std::vector<int>& aTri,
                          const std::vector<double>& aXYZ)
{
  int i0 = aTri[itri*3+0];
  int i1 = aTri[itri*3+1];
  int i2 = aTri[itri*3+2];
  CVector3 p0(aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2]);
  CVector3 p1(aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2]);
  CVector3 p2(aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2]);
  return (p1-p0)^(p2-p0);
}

std::ostream &operator<<(std::ostream &output, const std::vector<CVector3>& aV);
std::istream &operator>>(std::istream &input, std::vector<CVector3>& aV);

} // end namespace delfem2


#endif // VEC3_H
