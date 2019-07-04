#ifndef CAD3D_H
#define CAD3D_H

#include "delfem2/vec3.h"
#include "delfem2/vec2.h"
#include "delfem2/paramgeo_v23.h"
#include "delfem2/mshsrch_v3.h"

#include "delfem2/gl_color.h"

class CCad3D_Vertex{
public:
  CCad3D_Vertex(){
    pos = CVector3(0,0,0);
    isConst[0] = false;
    isConst[1] = false;
    isConst[2] = false;
  }
  CCad3D_Vertex(const CVector3& p){
    pos = p;
    isConst[0] = false;
    isConst[1] = false;
    isConst[2] = false;
  }
  void Draw(bool is_selected, int ielem, double view_height) const;
  void WriteFile(std::ofstream& fout) const{
    fout << "  " << pos << std::endl;
    fout << "  " << norm << std::endl;
  }
  void ReadFile(std::ifstream& fin){
    fin >> pos;
    fin >> norm;
  }
public:
  CVector3 pos;
  CVector3 norm;
  bool isConst[3];
  int iq_right; // id of mesh point right
  int iq_left;  // id of mesh point left
};


class CCad3D_Edge{
public:
  CCad3D_Edge(){
    this->iv0 = -1;
    this->iv1 = -1;
    this->inorm = 0;
    this->is_sim = false;
    r0 = 0.35;
    r1 = 0.35;
  }
  CCad3D_Edge(int icp0, int icp1, bool is_sim, int inorm){
    this->iv0 = icp0;
    this->iv1 = icp1;
    this->inorm = inorm;
    this->is_sim = is_sim;
    r0 = 0.35;
    r1 = 0.35;
  }
  CVector3 getNorm() const {
    if( inorm >=0 && inorm < 3){
      CVector3 norm(0,0,0);
      norm[inorm] = 1.0;
      return norm;
    }
    return CVector3(0,1,0);
  }
  void MovePoints(const std::vector<CCad3D_Vertex>& aVertex){
    if( aP.size()<2 ) return;
    p0 = aVertex[iv0].pos;
    p1 = aVertex[iv1].pos;
    const double len01 = (p1-p0).Length();
    CVector3 n0 = aVertex[iv0].norm;
    CVector3 n1 = aVertex[iv1].norm;
    CVector3 en = this->getNorm();
    CVector3 v0 = Cross(+n0,en); //if( v0*(p1-p0)<0 ){ v0*=-1; }
    CVector3 v1 = Cross(-n1,en); //if( v1*(p0-p1)<0 ){ v1*=-1; }
    q0 = p0 + v0*len01*r0;
    q1 = p1 + v1*len01*r1;
    const int ndiv = (int)aP.size()-1;
    for(int ip=0;ip<ndiv+1;++ip){
      double t = (double)ip/ndiv;
      aP[ip] = getPointCubicBezierCurve(t, p0, q0, q1, p1);
    }
  }
  void Initialize(const std::vector<CCad3D_Vertex>& aVertex, double elen){
    p0 = aVertex[iv0].pos;
    p1 = aVertex[iv1].pos;
    double len = Distance(p0, p1);
    const int ndiv = len/elen;
    aP.resize((ndiv+1));
    MovePoints(aVertex);
  }
  bool isPick(double& ratio, const CVector2& sp0, const float mMV[16], const float mPj[16]) const;
  bool GetParameterIntersection(double& t, const CVector3& org, const CVector3& nrm) const {
    return getParameterCubicBezier_IntersectionWithPlane(t, org,nrm, p0,q0,q1,p1);
  }
  CVector3 GetPosInEdge(double t) const {
    return getPointCubicBezierCurve(t, p0, q0, q1, p1);
  }
  CVector3 GetTangentInEdge(double t) const {
    return getTangentCubicBezierCurve(t, p0, q0, q1, p1);
  }
  CVector3 getVtxPos(bool is_root, int ioff) const {
    if( is_root ){
      if(      ioff == 0 ){ return p0; }
      else if( ioff == 1 ){ return q0; }
      else if( ioff == 2 ){ return q1; }
      else if( ioff == 3 ){ return p1; }
    }
    else{
      if(      ioff == 0 ){ return p1; }
      else if( ioff == 1 ){ return q1; }
      else if( ioff == 2 ){ return q0; }
      else if( ioff == 3 ){ return p0; }
    }
    return CVector3(0,0,0);
  }
  void DrawLine(bool is_picked, double view_height) const;
  void DrawHandler(int ielem_picked, double view_height) const;
  
  void WriteFile(std::ofstream& fout) const{
    fout << "  " << iv0 << " " << iv1 << std::endl;
    fout << "  " << is_sim << std::endl;
    fout << "  " << inorm << std::endl;
    fout << "  " << r0 << " " << r1 << std::endl;
  }
  void ReadFile(std::ifstream& fin){
    fin >> iv0 >> iv1;
    fin >> is_sim;
    fin >> inorm;
    fin >> r0 >> r1;
  }
public:
  int iv0;
  int iv1;
  int inorm;
  double r0;
  double r1;
  ///
  CVector3 p0,p1,q0,q1;
  std::vector<CVector3> aP;
  bool is_sim;
  std::vector<int> aIQ_RightLeft;
};

class CCad3D_Face{
public:
  CCad3D_Face(){
  }
  CCad3D_Face(const std::vector< std::pair<int,bool> >& aIE0){
    aIE = aIE0;
  }
  int findIIE_CP(int icp_in, const std::vector<CCad3D_Edge>& aEdge) const {
    for(int iie=0;iie<aIE.size();++iie){
      int ie0 = aIE[iie].first;
      assert( ie0>=0 && ie0<aEdge.size() );
      const CCad3D_Edge& e0 = aEdge[ie0];
      const bool dir0 = aIE[iie].second;
      int icp0 = (dir0) ? e0.iv0 : e0.iv1;
      if( icp0 == icp_in ){ return iie; }
    }
    return -1;
  }
  void Initialize
  (const std::vector<CCad3D_Vertex>& aCtrlPoint,
   const std::vector<CCad3D_Edge>& aEdge,
   double elen);
  bool isPick(const CVector3& org, const CVector3& dir){
    CPointElemSurf pes = intersect_Ray_MeshTri3D(org,dir, aTri,aXYZ);
    if( pes.itri != -1 ){ return true; }
    return false;
  }
  void MovePoints(const std::vector<CCad3D_Vertex>& aVertex,
                  const std::vector<CCad3D_Edge>& aEdge);
  void DrawFace() const;
  void DrawBackFace() const;
  void DrawEdge() const;
  void WriteFile(std::ofstream& fout) const{
    fout << "  " << aIE.size() << std::endl;
    for(int iie=0;iie<aIE.size();++iie){
      fout << "   " << aIE[iie].first << " " << aIE[iie].second << std::endl;
    }
  }
  void ReadFile(std::ifstream& fin){
    int nIE;
    fin >> nIE;
    aIE.resize(nIE);
    for(int iie=0;iie<nIE;++iie){
      int ie0;
      bool flg;
      fin >> ie0 >> flg;
      aIE[iie].first = ie0;
      aIE[iie].second = flg;
    }
  }
public:
  class CFacePointInfo{
  public:
    int itype; // 0:vtx, 1:edge, 2:face
    int iv;
    int ie;
    int iep;
    CVector3 n;
    std::vector<double> aW0;
    std::vector<double> aW1;
    int iq_right;
    int iq_left;
  };
public:
  std::vector< std::pair<int,bool> > aIE; // index of edge, is this edge ccw?
  ///
  //CVector3 norm;
  std::vector<double> aXYZ;
  std::vector<double> aNorm;
  std::vector<unsigned int> aTri;
  std::vector<CFacePointInfo> aPInfo;
};

int AddPointEdge
(int ie_div, double ratio_edge,
 std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen);

void ConectEdge
(int iv0, int iv1, int iface_div, int inorm_new,
 std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen);

void MakeItSmooth
(std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace);

void findEdgeGroup
(std::vector< std::pair<int,bool> >& aIE,
 int iedge0,
 std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge);

void AddSphere_XSym
(std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen);
void AddSphere_ZSym
(std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen);
void AddTorus_XSym
(std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen);

void AddCube
(std::vector<CCad3D_Vertex>& aVertex,
std::vector<CCad3D_Edge>& aEdge,
std::vector<CCad3D_Face>& aFace,
double elen);

bool FindFittingPoint
(CVector2& p2d_near,
 CVector2& p2d_norm,
 const CVector2& p2d_org,
 const std::vector<CVector2>& aP2D,
 bool isConstX, bool isConstY,
 double half_view_height);

std::vector<int> getPointsInEdges
(const std::vector< std::pair<int,bool > >& aIE_picked,
 const std::vector<CCad3D_Edge>& aEdge);

bool MovePointsAlongSketch
(std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 const std::vector<CVector2>& aStroke,
 const std::vector< std::pair<int,bool > >& aIE_picked,
 const CVector3& plane_org, int inorm,
 float mMV[16], float mPj[16], double view_height);

void DivideFace
(int ifc,
 const CVector3& org, int inorm,
 std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen);

void BuildTriMesh
(std::vector<double>& aXYZ,
 std::vector<unsigned int>& aTri,
 std::vector<int>& aTriSurRel,
 std::vector<double>& aNorm,
 ////
 std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 int isym);

void UpdateTriMesh
(std::vector<double>& aXYZ,
 std::vector<double>& aNorm,
 const std::vector<int>& aTri,
 const std::vector<CCad3D_Vertex>& aVertex,
 const std::vector<CCad3D_Edge>& aEdge,
 const std::vector<CCad3D_Face>& aFace,
 int isym);


class CCad3D
{
public:
  CCad3D(){
    color_face = CColor(0.75164f, 0.60648f, 0.22648f, 1.0f);
    color_face_selected = CColor(1.0, 1.0, 0.0, 1);
    isym = 0;
    elen = 0.1;
    plane_inorm = -1; // if it is 0,1,2 it shouws plane
    plane_sizeX = 1.0; // norm[(inrom+1)%3]
    plane_sizeY = 1.0; // norm[(inrom+2)%3]
    imode_edit = EDIT_NONE;
    ivtx_picked = -1;
    iedge_picked = -1;
    ielem_edge_picked = 0; // 0:JustPicked, 1:PlaneEdgePickedd, 2:BezierHandlerA, 3:BezierHandlerB
    iface_picked = -1;
  }
  void Clear(){
    aVertex.clear();
    aEdge.clear();
    aFace.clear();
    ////
    aXYZ.clear();
    aTri.clear();
    aNorm.clear();
    aTriSurRel.clear();
  }
  void Initialize_Torus(){
    Clear();
    AddTorus_XSym(aVertex,aEdge,aFace,elen);
    isym = 0;
    BuildTriMesh(aXYZ,aTri,aTriSurRel,aNorm, aVertex,aEdge,aFace, isym);
  }
  void Initialize_Sphere(){
    Clear();
    AddSphere_XSym(aVertex,aEdge,aFace,elen);
    BuildTriMesh(aXYZ,aTri,aTriSurRel,aNorm, aVertex,aEdge,aFace, isym);
  }
  void Initialize_Cube(){
    Clear();
    AddCube(aVertex, aEdge, aFace, elen);
    BuildTriMesh(aXYZ, aTri, aTriSurRel, aNorm, aVertex, aEdge, aFace, isym);
  }
  void Pick
  (const CVector3& src_pick, const CVector3& dir_pick,
   CVector2 sp0, float mMV[16], float mPj[16],
   double view_height);
  
  void DrawFace_LeftRight() const;
  void DrawVtxEdgeHandler(double view_height) const;
  void DrawFace_RightSelected(bool is_edge) const;
  
  bool ReflectChangeForCurveAndSurface(std::vector<int>& aIsMoved_Edge,
                                       const std::vector<int>& aIsMoved_Vtx);
  void MouseDown(const CVector3& src_pick, const CVector3& dir_pick,
                 const CVector2& sp0, float mMV[16], float mPj[16],
                 double view_height);
  void MouseUp(float mMV[16], float mPj[16], double view_height);
  bool MouseMotion(const CVector3& src_pick, const CVector3& dir_pick,
                   const CVector2& sp0, const CVector2& sp1,
                   float mMV[16], float mPj[16]);
  
  void WriteFile(std::ofstream& fout) const;
  void ReadFile(std::ifstream& fin);
  void ReadFile(const std::string& fname){
    std::ifstream fin;
    fin.open(fname.c_str());
    this->ReadFile(fin);
  }
  bool isSym(int iv) const;
  void ChangeEdgeLength(double elen){
    this->elen = elen;
    for(int ie=0;ie<aEdge.size();ie++){
      aEdge[ie].Initialize(aVertex,elen); // ie0+0
    }
    for(int ifc=0;ifc<aFace.size();ifc++){
      aFace[ifc].Initialize(aVertex,aEdge, elen); // ie0+0
    }
    BuildTriMesh(aXYZ,aTri,aTriSurRel,aNorm, aVertex,aEdge,aFace, isym);
  }
public:
  ////////////////////////////////
  // fundamental data
  int isym;
  std::vector<CCad3D_Vertex> aVertex;
  std::vector<CCad3D_Edge> aEdge;
  std::vector<CCad3D_Face> aFace;

  ////////////////////////////////
  // aux data
  double elen;
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  std::vector<double> aNorm;
  std::vector<int> aTriSurRel;

  ////////////////////////////////
  // pick related
  int ivtx_picked;
  int ielem_vtx_picked;
  
  // edge related pick
  int iedge_picked;
  int ielem_edge_picked; // 0:JustPicked, 1:PlaneEdgePickedd, 2:BezierHandlerA, 3:BezierHandlerB
  double ratio_edge_picked;
  std::vector< std::pair<int,bool > > aIE_picked;
  
  int plane_inorm; // if it is 0,1,2 it shouws plane
  CVector3 plane_org;
  double plane_sizeX; // norm[(inrom+1)%3]
  double plane_sizeY; // norm[(inrom+2)%3]
  std::vector<CVector2> aStroke;
  
  // face related pick
  int iface_picked;
  
  enum EDIT_MODE {
    EDIT_NONE,
    EDIT_MOVE,
    EDIT_ADD_POINT_EDGE,
    EDIT_SKETCH,
    EDIT_ADD_CROSS_SECTION,
  } imode_edit;
  
  ////////////////////////////////
  // viz related
  CColor color_face;
  CColor color_face_selected;
};



#endif /* cad3d_h */
