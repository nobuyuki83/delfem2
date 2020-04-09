#include <iostream>
#include <fstream>
#include <math.h>
#include "delfem2/pgeo.h"
#include "delfem2/funcs.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// -----------------------------

int getEntity
(std::string& type,
 std::string& arg,
 const std::string& s)
{
  int is;
  for(is=0;is<s.size()-1;++is){
    if( s[is] == '#' ) break;
  }
  if( s[is] != '#' ) return -1;
  int ie;
  for(ie=is+1;ie<s.size();++ie){
    int iascii = (int)s[ie]-48;
    if( iascii < 0 || iascii > 9 ) break;
  }
  int ID0;
  {
    char sdist[256];
    strncpy(sdist, s.data()+is+1, ie-is-1);
    sdist[ie-is-1] = '\0';
    ID0 = std::stoi(sdist);
  }
  ////
  int ic;
  for(ic=ie;ic<s.size();++ic){
    if( s[ic] != ' ' && s[ic] != '=' ) break;
  }
  int id;
  for(id=ic;id<s.size();++id){
    int iascii = (int)s[id];
    if( !dfm2::isAlphabetUpper(iascii) && !dfm2::isNumeric(iascii) && (iascii != 95)  ) break;
//    if( s[id] == '(' ) break;
  }
  /*
  id--;
  for(;id>=0;--id){
    if( s[id] != ' ' ){ id++; break; }
  }
   */
  {
    char sdist[256];
    strncpy(sdist, s.data()+ic, id-ic);
    sdist[id-ic] = '\0';
//    std::cout << s << "  " << ic << " " << id << " ###" << sdist << "##" << std::endl;
    type = std::string(sdist);
  }
  
  arg = std::string(s.data()+id,s.data()+s.size());
  
  return ID0;
}

class CStep_Elem {
public:
  virtual void Init(int id, const std::string& stype, const std::string& arg) = 0;
  virtual void Draw() const = 0;
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind) = 0;
public:
  int id;
  std::string stype;
};

class CStep_CartesianPoint: public CStep_Elem
{
public:
  virtual void Init(int id, const std::string& stype, const std::string& arg){
    this->id = id;
    this->stype = stype;
    //
    std::string s = dfm2::Get_Parentheses(arg, "()");
    s =  dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s,',');
    p.p[0] = stod(aToken[0]);
    p.p[1] = stod(aToken[1]);
    p.p[2] = stod(aToken[2]);
//    std::cout << p.x << " " << p.y << " " << p.z << std::endl;
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind){}
  virtual void Draw() const {
    ::glColor3d(1,0,0);
    ::glPointSize(1);
    ::glBegin(GL_POINTS);
//    ::glVertex3d(p.x, p.y, p.z);
    ::glEnd();
  }
public:
  dfm2::CVec3d p;
};

class CStep_Direction: public CStep_Elem
{
public:
  void Init(int id, const std::string& stype, const std::string& arg){
    this->id = id;
    this->stype = stype;
    //
    std::string s =  dfm2::Get_Parentheses(arg, "()");
    s =  dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s,',');
    dir.p[0] = stod(aToken[0]);
    dir.p[1] = stod(aToken[1]);
    dir.p[2] = stod(aToken[2]);
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind){}
  virtual void Draw() const {}
public:
  dfm2::CVec3d dir;
};

class CStep_Vector: public CStep_Elem
{
  void Init(int id, const std::string& stype, const std::string& arg){
    this->id = id;
    this->stype = stype;
    //
    std::string s =  dfm2::Get_Parentheses(arg, "()");
    s = dfm2::RemoveSpace(s);
    std::vector<std::string> aToken = dfm2::Split(s,',');
    std::string s1(aToken[1].begin()+1,aToken[1].end());
    idDir = std::stoi(s1);
    len = stod(aToken[2]);
    pDir = 0;
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind){
    const int iDir = mapId2Ind[idDir];
    pDir = (CStep_Direction*)apStepElem[iDir];
  }
  virtual void Draw() const {}
public:
  double len;
  int idDir;
  ////
  const CStep_Direction* pDir;
};

class CStep_VertexPoint: public CStep_Elem
{
public:
  void Init(int id, const std::string& stype, const std::string& arg){
    this->id = id;
    this->stype = stype;
    //
    std::string s =  dfm2::Get_Parentheses(arg, "()");
    std::vector<std::string> aToken = dfm2::Split(s,',');
    s = dfm2::RemoveSpace(aToken[1]);
    s = std::string(s.begin()+1,s.end());
    this->idCP = std::stoi(s);
    pCP = 0;
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind){
    int icp = mapId2Ind[idCP];
    pCP = (CStep_CartesianPoint*)apStepElem[icp];
  }
  virtual void Draw() const {
    const dfm2::CVec3d& p0 = pCP->p;
    dfm2::opengl::myGlTranslate(+p0);
    ::glColor3d(1,0,1);
    dfm2::opengl::DrawSphere(16, 16);
    dfm2::opengl::myGlTranslate(-p0);
  }
public:
  int idCP;
  const CStep_CartesianPoint* pCP = 0;
};

class CStep_Axis2Placement3D: public CStep_Elem
{
public:
  void Init(int id, const std::string& stype, const std::string& arg){
    this->id = id;
    this->stype = stype;
    //
    std::string s = dfm2::RemoveSpace(arg);
    s =  dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s,',');
    assert( aToken.size() == 4 );
    {
      std::string s1 = aToken[1];
      s1 = std::string(s1.begin()+1,s1.end());
      this->idCP = std::stoi(s1);
      pCP = 0;
    }
    {
      std::string s2 = aToken[2];
      if( s2 == "$" ){ this->idDir1 = -1; }
      else{
        s2 = std::string(s2.begin()+1,s2.end());
        this->idDir1 = std::stoi(s2);
      }
      pDir1 = 0;
    }
    {
      std::string s3 = aToken[3];
      if( s3 == "$" ){ this->idDir2 = -1; }
      else{
        s3 = std::string(s3.begin()+1,s3.end());
        this->idDir2 = std::stoi(s3);
      }
      pDir2 = 0;
    }
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind){
    const int icp = mapId2Ind[idCP];
    const int idir1 = mapId2Ind[idDir1];
    const int idir2 = mapId2Ind[idDir2];
    pCP = (CStep_CartesianPoint*)apStepElem[icp];
    pDir1 = (CStep_Direction*)apStepElem[idir1];
    pDir2 = (CStep_Direction*)apStepElem[idir2];
  }
  virtual void Draw() const {
    /*
    const  dfm2::CVec3d& p0 = pCP->p;
    const  dfm2::CVec3d& d1 = pDir1->dir;
    const  dfm2::CVec3d& d2 = pDir2->dir;
    ::myGlTranslate(+p0);
    ::glColor3d(0,0,1);
    DrawSphere(16, 16);
    ::glColor3d(0,1,1);
    ::DrawArrow( dfm2::CVec3d(0,0,0), d1*5);
    ::glColor3d(1,0,1);
    ::DrawArrow( dfm2::CVec3d(0,0,0), d2*5);
    ::myGlTranslate(-p0);
     */
  }
public:
  int idCP, idDir1, idDir2;
  const CStep_CartesianPoint* pCP;
  const CStep_Direction* pDir1;
  const CStep_Direction* pDir2;
};

class CStep_Curve: public CStep_Elem
{
public:
  virtual double GetParameter(const dfm2::CVec3d& p) const = 0;
  virtual void SampleCurve(std::vector<dfm2::CVec3d>& polyLine, double r1, double r2, int nsmpl) const = 0;
};

class CStep_Line: public CStep_Curve
{
  void Init(int id, const std::string& stype, const std::string& arg){
    this->id = id;
    this->stype = stype;
    //
    std::string s = dfm2::RemoveSpace(arg);
    s =  dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s,',');
    assert( aToken.size() == 3 );
    {
      std::string s1 = std::string(aToken[1].begin()+1,aToken[1].end());
      this->idCP = std::stoi(s1);
      pCP = 0;
    }
    {
      std::string s2 = std::string(aToken[2].begin()+1,aToken[2].end());
      this->idVec = std::stoi(s2);
      pVec = 0;
    }
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind){
    const int iCP = mapId2Ind[idCP];
    assert( iCP != -1 );
    assert( apStepElem[iCP]->stype == "CARTESIAN_POINT" );
    pCP = (CStep_CartesianPoint*)apStepElem[iCP];
    ////
    const int iVec = mapId2Ind[idVec];
    assert( iVec != -1 );
    assert( apStepElem[iVec]->stype == "VECTOR" );
    pVec = (CStep_Vector*)apStepElem[iVec];
  }
  virtual void Draw() const {}
  //
  virtual double GetParameter(const dfm2::CVec3d& p) const{
    const dfm2::CVec3d& cp = pCP->p;
    const dfm2::CVec3d& d = pVec->pDir->dir;
    const double l = pVec->len;
    const dfm2::CVec3d v = l*d;
    return (p-cp)*v/v.DLength();
  }
  virtual void SampleCurve(std::vector<dfm2::CVec3d>& polyLine, double r1, double r2, int nsmpl) const{
    const dfm2::CVec3d& cp = pCP->p;
    const dfm2::CVec3d& d = pVec->pDir->dir;
    const double l = pVec->len;
    polyLine.clear();
    for(int is=0;is<nsmpl;++is){
      double r = (r2-r1)/(nsmpl-1)*is + r1;
      dfm2::CVec3d p = cp + r*l*d;
      polyLine.push_back(p);
    }
  }
public:
  int idCP;
  int idVec;
  ////
  const CStep_CartesianPoint* pCP;
  const CStep_Vector* pVec;
};

class CStep_Circle: public CStep_Curve
{
public:
  void Init(int id, const std::string& stype, const std::string& arg){
    this->id = id;
    this->stype = stype;
    //
    std::string s = dfm2::RemoveSpace(arg);
    s =  dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s,',');
    assert( aToken.size() == 3 );
    {
      s = std::string(aToken[1].begin()+1,aToken[1].end());
      this->idA2P3D = std::stoi(s);
      pA2P3D = 0;
    }
    {
      this->radius = std::stof(aToken[2]);
    }
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind){
    int iA2P3D = mapId2Ind[idA2P3D];
    pA2P3D = (CStep_Axis2Placement3D*)apStepElem[iA2P3D];
  }
  virtual void Draw() const {}
  /////
  virtual double GetParameter(const dfm2::CVec3d& p) const{
    const dfm2::CVec3d& c = pA2P3D->pCP->p;
    const dfm2::CVec3d& n = pA2P3D->pDir1->dir;
    const dfm2::CVec3d ex = (pA2P3D->pDir2->dir).Normalize();
    const dfm2::CVec3d ey = Cross(n,ex).Normalize();
    return atan2( (p-c)*ey, (p-c)*ex );
  }
  virtual void SampleCurve(std::vector<dfm2::CVec3d>& polyLine, double r1, double r2, int nsmpl) const{
    const dfm2::CVec3d& c = pA2P3D->pCP->p;
    const dfm2::CVec3d& n = pA2P3D->pDir1->dir;
    const dfm2::CVec3d ex = (pA2P3D->pDir2->dir).Normalize();
    const dfm2::CVec3d ey = Cross(n,ex).Normalize();
    if( r1 > r2 ){ r2 += M_PI*2; }
    polyLine.clear();
    for(int is=0;is<nsmpl;++is){
      double r = (r2-r1)/(nsmpl-1)*is + r1;
      dfm2::CVec3d p = c + radius*cos(r)*ex + radius*sin(r)*ey;
      polyLine.push_back(p);
    }
  }
public:
  int idA2P3D;
  double radius;
  ////
  const CStep_Axis2Placement3D* pA2P3D;
};

class CStep_BSplineCurveWithKnots: public CStep_Curve
{
public:
  void Init(int id, const std::string& stype, const std::string& arg)
  {
    this->id = id;
    this->stype = stype;
    /////
    aIdCP.clear();
    aKnotMulti.clear();
    std::string s =  dfm2::Get_Parentheses(arg, "()");
    s = dfm2::RemoveSpace(s);
    std::vector<std::string> aToken = dfm2::Split_Parentheses(s,',',"()");
    if( aToken.size() != 9 ){
      std::cout << "something is wrong" << std::endl;
    }
    {
      ndegree = stoi(aToken[1]);
    }
    {
      std::string sID = aToken[2];
      sID =  dfm2::Get_Parentheses(sID, "()");
      std::vector<std::string> at = dfm2::Split(sID,',');
      for(int i=0;i<at.size();++i){
        const int id = stoi(std::string(at[i].begin()+1,at[i].end()));
        aIdCP.push_back(id);
      }
    }
    {
      std::string sKM = aToken[6];
      sKM =  dfm2::Get_Parentheses(sKM, "()");
      std::vector<std::string> at = dfm2::Split(sKM,',');
      for(int i=0;i<at.size();++i){
        const int ikm = stoi(std::string(at[i].begin(),at[i].end()));
        aKnotMulti.push_back(ikm);
      }
    }
    {
      std::string sK = aToken[7];
      sK =  dfm2::Get_Parentheses(sK, "()");
      std::vector<std::string> at = dfm2::Split(sK,',');
      for(int i=0;i<at.size();++i){
        const double k = stod(std::string(at[i].begin(),at[i].end()));
        aKnot.push_back(k);
      }
    }
    dfm2::FlatKnot(aKnotFlat, aKnotMulti, aKnot);
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind){
    apCP.clear();
    for(int iicp=0;iicp<aIdCP.size();++iicp){
      int idcp = aIdCP[iicp];
      int icp = mapId2Ind[idcp];
      assert( apStepElem[icp]->stype == "CARTESIAN_POINT" );
      apCP.push_back( (CStep_CartesianPoint*)apStepElem[icp] );
    }
  }
  virtual void Draw() const {}
  /////
  virtual double GetParameter(const dfm2::CVec3d& p) const{ return 0; }
  virtual void SampleCurve(std::vector<dfm2::CVec3d>& polyLine, double r1, double r2, int nsmpl) const{
    std::vector<dfm2::CVec3d> aPosCP;
    for(int iicp=0;iicp<apCP.size();++iicp){
      aPosCP.push_back( apCP[iicp]->p );
    }
    SampleBSpline(polyLine, nsmpl, ndegree, aKnotFlat, aPosCP);
  }
public:
  std::vector<int> aIdCP;
  std::vector<int> aKnotMulti; // multiplicity
  std::vector<double> aKnot;
  int ndegree = 3;
  std::vector<double> aKnotFlat;
  ////
  std::vector<CStep_CartesianPoint*> apCP;
};

class CStep_EdgeCurve: public CStep_Elem
{
public:
  void Init(int id, const std::string& stype, const std::string& arg)
  {
    this->id = id;
    this->stype = stype;
    //
    std::string s =  dfm2::Get_Parentheses(arg, "()");
    s = dfm2::RemoveSpace(s);
    std::vector<std::string> aToken = dfm2::Split(s,',');
    assert(aToken.size()==5);
    {
      std::string scp1 = std::string(aToken[1].begin()+1,aToken[1].end());
      idCP1 = std::stoi(scp1);
      pCP1 = 0;
    }
    {
      std::string scp2 = std::string(aToken[2].begin()+1,aToken[2].end());
      idCP2 = std::stoi(scp2);
      pCP2 = 0;
    }
    {
      std::string scur = std::string(aToken[3].begin()+1,aToken[3].end());
      idCurve = std::stoi(scur);
      pCurve = 0;
    }
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind)
  {
    {
      const int iCP1 = mapId2Ind[idCP1];
      assert( iCP1 != -1 );
      assert( apStepElem[iCP1]->stype == "VERTEX_POINT" );
      pCP1 = (CStep_VertexPoint*)apStepElem[iCP1];
    }
    {
      const int iCP2 = mapId2Ind[idCP2];
      assert( iCP2 != -1 );
      assert( apStepElem[iCP2]->stype == "VERTEX_POINT" );
      pCP2 = (CStep_VertexPoint*)apStepElem[iCP2];
    }
    {
      const int iCurve = mapId2Ind[idCurve];
      if( iCurve != -1 ){ pCurve = (CStep_Curve*)apStepElem[iCurve]; }
      else{ pCurve = 0; }
    }
  }
  virtual void Draw() const {
    ::glColor3d(0,0,0);
    dfm2::opengl::drawPolyLine3D(polyLine);
  }
  /////
  void Sample(){
    dfm2::CVec3d p1(pCP1->pCP->p);
    dfm2::CVec3d p2(pCP2->pCP->p);
    polyLine.clear();
    if( pCurve != 0 ){
      double r1 = pCurve->GetParameter(p1);
      double r2 = pCurve->GetParameter(p2);
      pCurve->SampleCurve(polyLine,r1,r2,10);
    }
  }
public:
  int idCP1, idCP2, idCurve;
  const CStep_VertexPoint* pCP1;
  const CStep_VertexPoint* pCP2;
  const CStep_Curve* pCurve;
  ////
  std::vector<dfm2::CVec3d> polyLine;
};

class CStep_OrientedEdge: public CStep_Elem
{
public:
  void Init(int id, const std::string& stype, const std::string& arg)
  {
    this->id = id;
    this->stype = stype;
    //////
    std::string s = dfm2::RemoveSpace(arg);
    s =  dfm2::Get_Parentheses(s, "()" );
    std::vector<std::string> aToken = dfm2::Split(s,',');
    assert(aToken.size()==5);
    {
      std::string s1(aToken[3].begin()+1,aToken[3].end());
      idEdgeCurve = std::stoi(s1);
    }
    {
      if( aToken[4] == ".T." ){ isSameDir = true; }
      else{                     isSameDir = false; }
    }
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind){
    const int iEC = mapId2Ind[idEdgeCurve];
    pEC = (CStep_EdgeCurve*)apStepElem[iEC];
  }
  virtual void Draw() const {}
public:
  int idEdgeCurve;
  bool isSameDir;
  const CStep_EdgeCurve* pEC;
};

class CStep_EdgeLoop: public CStep_Elem
{
public:
  void Init(int id, const std::string& stype, const std::string& arg)
  {
    this->id = id;
    this->stype = stype;
    //////
    std::string s = dfm2::RemoveSpace(arg);
    s =  dfm2::Get_Parentheses(s, "()" );
    s =  dfm2::Get_Parentheses(s, "()" );
    std::vector<std::string> aToken = dfm2::Split(s,',');
    aIdOE.clear();
    for(int ioe=0;ioe<aIdOE.size();++ioe){
      int idOE = std::stoi( std::string(aToken[ioe].begin()+1,aToken[ioe].end()) );
      aIdOE.push_back(idOE);
    }
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind){
    apOE.clear();
    for(int iioe=0;iioe<aIdOE.size();++iioe){
      int idOE = aIdOE[iioe];
      const int iOE = mapId2Ind[idOE];
      assert( iOE != -1 );
      assert( apStepElem[iOE]->stype == "ORIENTED_EDGE" );
      apOE.push_back( (CStep_OrientedEdge*)apStepElem[iOE] );
    }
  }
  virtual void Draw() const {}
public:
  std::vector<int> aIdOE;
  std::vector<CStep_OrientedEdge*> apOE;
};

class CStep_FaceOuterBound: public CStep_Elem
{
public:
  void Init(int id, const std::string& stype, const std::string& arg)
  {
    this->id = id;
    this->stype = stype;
    //////
    std::string s = dfm2::RemoveSpace(arg);
    s =  dfm2::Get_Parentheses(s, "()" );
    std::vector<std::string> aToken = dfm2::Split(s,',');
    {
      std::string s1 = std::string(aToken[1].begin()+1,aToken[1].end());
      idEL = std::stoi(s1);
      pEL = 0;
    }
    {
      std::string s1 = aToken[2];
      if( s1 == ".T." ){ orientation = true; }
      else{              orientation = false; }
    }
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind){
    const int iEL = mapId2Ind[idEL];
    assert( iEL!=-1 );
    assert( apStepElem[iEL]->stype == "EDGE_LOOP" );
    pEL = (CStep_EdgeLoop*)apStepElem[iEL];
  }
  virtual void Draw() const {}
public:
  int idEL;
  bool orientation;
  ///
  const CStep_EdgeLoop* pEL;
};

class CStep_Surface: public CStep_Elem
{
  virtual void Project(double& r0, double& r1,
                       const dfm2::CVec3d& p0) const = 0;
public:
};

class CStep_Plane : public CStep_Surface
{
public:
  void Init(int id, const std::string& stype, const std::string& arg)
  {
    this->id = id;
    this->stype = stype;
    //////
    std::string s = dfm2::RemoveSpace(arg);
    s =  dfm2::Get_Parentheses(s, "()" );
    std::vector<std::string> aToken = dfm2::Split(s,',');
    assert( aToken.size() == 2 );
    {
      std::string s1 = std::string(aToken[1].begin()+1,aToken[1].end());
      idA2P3D = std::stoi(s1);
      pA2P3D = 0;
    }
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind){
    const int iA2P3D = mapId2Ind[idA2P3D];
    assert( iA2P3D!=-1 );
    assert( apStepElem[iA2P3D]->stype == "AXIS2_PLACEMENT_3D" );
    pA2P3D = (CStep_Axis2Placement3D*)apStepElem[iA2P3D];
  }
  virtual void Draw() const {
     dfm2::CVec3d org = pA2P3D->pCP->p;
     dfm2::CVec3d n = pA2P3D->pDir1->dir;
     dfm2::CVec3d ex = pA2P3D->pDir2->dir;
     dfm2::CVec3d ey = Cross(n,ex);
    double l = 100;
//    ::glEnable(GL_LIGHTING);
    ::glBegin(GL_LINE_LOOP);
    dfm2::opengl::myGlNormal(n);
    dfm2::opengl::myGlVertex(org-l*ex-l*ey);
    dfm2::opengl::myGlVertex(org+l*ex-l*ey);
    dfm2::opengl::myGlVertex(org+l*ex+l*ey);
    dfm2::opengl::myGlVertex(org-l*ex+l*ey);
    ::glEnd();
  }
  virtual void Project(double& r0, double& r1,
                       const  dfm2::CVec3d& p0) const
  {
    const  dfm2::CVec3d& org = pA2P3D->pCP->p;
    const  dfm2::CVec3d& n = pA2P3D->pDir1->dir;
    const  dfm2::CVec3d& ex = pA2P3D->pDir2->dir;
    const  dfm2::CVec3d& ey = Cross(n,ex);
    r0 = (p0-org)*ex;
    r1 = (p0-org)*ey;
  }
public:
  int idA2P3D;
  const CStep_Axis2Placement3D* pA2P3D;
};


class CStep_AdvancedFace: public CStep_Elem
{
public:
  void Init(int id, const std::string& stype, const std::string& arg)
  {
    this->id = id;
    this->stype = stype;
    //
    std::string s = dfm2::RemoveSpace(arg);
    std::vector<std::string> aToken;
    {
      std::string s0 = dfm2::Get_Parentheses(s, "()");
      aToken = dfm2::Split_Parentheses(s0,',', "()");
      for(int it=0;it<aToken.size();++it){
        std::cout << it << " " << aToken[it] << "   " << arg << " " << s << std::endl;
      }
    }
    assert( aToken.size() == 4 );
    {
      std::string s1 = std::string(aToken[2].begin()+1,aToken[2].end());
      idSurf = std::stoi(s1);
      pSurf = 0;
    }
    {
      s =  dfm2::Get_Parentheses(aToken[1], "()");
      aToken = dfm2::Split(s,',');
      aIdFOB.clear();
      for(int it=0;it<aToken.size();++it){
        s = std::string(aToken[it].begin()+1,aToken[it].end());
        int idFOB = std::stoi(s);
        aIdFOB.push_back(idFOB);
      }
    }
  }
  virtual void SetPtr(const std::vector<CStep_Elem*>& apStepElem,
                      const std::vector<int>& mapId2Ind){
    apFOB.clear();
    for(int ifob=0;ifob<aIdFOB.size();++ifob){
      int idFOB = aIdFOB[ifob];
      const int iFOB = mapId2Ind[idFOB];
      if( iFOB == -1 ){ apFOB.push_back(0); continue; }
      assert( iFOB!=-1 );
      assert( apStepElem[iFOB]->stype == "FACE_OUTER_BOUND" );
      apFOB.push_back( (CStep_FaceOuterBound*)apStepElem[iFOB] );
    }
  }
  virtual void Draw() const {}
  ///
  void SampleSurface() {
  }
public:
  std::vector<int> aIdFOB;
  int idSurf;
  CStep_Surface* pSurf;
  std::vector<CStep_FaceOuterBound*> apFOB;
  /////
  std::vector<double> aXYZ;
  std::vector<int> aTri;
};

void LoadStep
(const std::string& fname,
 std::vector<CStep_Elem*>& apStepElem)
{
  std::cout << "load Step" << std::endl;
  {
    apStepElem.clear();
    std::ifstream fin(fname.c_str());
    std::string s;
    s.reserve(1024*8);
    while(getline(fin,s,';')){
      std::string type, arg;
      int id = getEntity(type,arg,s);
      std::cout << id << " " << type << " " << arg << " " << s << std::endl;
      CStep_Elem* pSE = 0;
      if( type == "CARTESIAN_POINT"           ){ pSE = new CStep_CartesianPoint(); }
      if( type == "DIRECTION"                 ){ pSE = new CStep_Direction(); }
      if( type == "VERTEX_POINT"              ){ pSE = new CStep_VertexPoint(); }
      if( type == "VECTOR"                    ){ pSE = new CStep_Vector(); }
      if( type == "AXIS2_PLACEMENT_3D"        ){ pSE = new CStep_Axis2Placement3D(); }
      if( type == "LINE"                      ){ pSE = new CStep_Line(); }
      if( type == "CIRCLE"                    ){ pSE = new CStep_Circle(); }
      if( type == "B_SPLINE_CURVE_WITH_KNOTS" ){ pSE = new CStep_BSplineCurveWithKnots(); }
      if( type == "EDGE_CURVE"                ){ pSE = new CStep_EdgeCurve(); }
      if( type == "ORIENTED_EDGE"             ){ pSE = new CStep_OrientedEdge(); }
      if( type == "EDGE_LOOP"                 ){ pSE = new CStep_EdgeLoop(); }
      if( type == "FACE_OUTER_BOUND"          ){ pSE = new CStep_FaceOuterBound(); }
      if( type == "PLANE"                     ){ pSE = new CStep_Plane(); }
      if( type == "ADVANCED_FACE"             ){ pSE = new CStep_AdvancedFace(); }
      if( pSE != 0 ){
        pSE->Init(id,type,arg);
        apStepElem.push_back(pSE);
      }
    }
  }
  ////////////////////////////////////////////////
  { // id 2 ind
    std::vector<int> mapId2Ind;
    for(int ise=0;ise<apStepElem.size();++ise){
      int id = apStepElem[ise]->id;
      if( id >= mapId2Ind.size() ){
        mapId2Ind.resize(id+1,-1);
      }
      mapId2Ind[id] = ise;
    }
    for(int ise=0;ise<apStepElem.size();++ise){
      apStepElem[ise]->SetPtr(apStepElem,mapId2Ind);
    }
    for(int ise=0;ise<apStepElem.size();++ise){
      if( apStepElem[ise]->stype == "EDGE_CURVE" ){
        ((CStep_EdgeCurve*)apStepElem[ise])->Sample();
//        std::cout << "sample " << ise << std::endl;
      }
    }
  }
}

/*
void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
  {
    case 'q':
    case 'Q':
      break;
    case 'a': //
      break;
    case ' ':
    {
      static int ifile = 0;
      ifile++;
      if( ifile >= 3 ){ ifile=0; }
      ////
      std::string fpath;
      if(      ifile == 0 ){ fpath=std::string("grabcad/hook/Hook.STEP"); }
      else if( ifile == 1 ){ fpath=std::string("grabcad/bumper/BUMPER.stp"); }
      else if( ifile == 2 ){ fpath=std::string("grabcad/anchor/bolt.STEP");  }
      LoadStep(fpath,apStepElem);
    }
  }
  ::glutPostRedisplay();
}
 */

int main(int argc, char* argv[])
{
  std::vector<CStep_Elem*> apStepElem;
  //  LoadStep("Hook.STEP",aCP,aBSCWK);
  LoadStep(std::string(PATH_INPUT_DIR)+"/bolt.STEP",apStepElem);
  //  LoadStep("grabcad/bumper/BUMPER.stp",aCP,aBSCWK);
  //  LoadStep("bolt.STEP",aCP,aBSCWK);
  
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 100;
  
  while(!glfwWindowShouldClose(viewer.window)){
    viewer.DrawBegin_oldGL();
    for(int ipse=0;ipse<apStepElem.size();++ipse){
      apStepElem[ipse]->Draw();
    }
    viewer.DrawEnd_oldGL();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
