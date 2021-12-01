
#include <iostream>
#include <cmath>
#include <fstream>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/geo_bspline.h"
#include "delfem2/str.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/v3q.h"

#ifndef M_PI
#  define M_PI 3.141592653589793
#endif

namespace dfm2 = delfem2;

// -----------------------------

int getEntity(
    std::string &type,
    std::string &arg,
    const std::string &s) {
  int is;
  for (is = 0; is < (int) s.size() - 1; ++is) {
    if (s[is] == '#') break;
  }
  if (s[is] != '#') return -1;
  int ie;
  for (ie = is + 1; ie < (int) s.size(); ++ie) {
    int iascii = (int) s[ie] - 48;
    if (iascii < 0 || iascii > 9) break;
  }
  int ID0;
  {
    char sdist[256];
    strncpy(sdist, s.data() + is + 1, ie - is - 1);
    sdist[ie - is - 1] = '\0';
    ID0 = std::stoi(sdist);
  }
  //
  int ic;
  for (ic = ie; ic < (int) s.size(); ++ic) {
    if (s[ic] != ' ' && s[ic] != '=') break;
  }
  int id;
  for (id = ic; id < (int) s.size(); ++id) {
    int iascii = (int) s[id];
    if (!dfm2::isAlphabetUpper(iascii) && !dfm2::isNumeric(iascii) && (iascii != 95)) break;
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
    strncpy(sdist, s.data() + ic, id - ic);
    sdist[id - ic] = '\0';
//    std::cout << s << "  " << ic << " " << id << " ###" << sdist << "##" << std::endl;
    type = std::string(sdist);
  }

  arg = std::string(s.data() + id, s.data() + s.size());

  return ID0;
}

class CStep_Elem {
 public:
  virtual void Init(int id, const std::string &stype, const std::string &arg) = 0;
  virtual void Draw() const = 0;
  virtual void SetPtr(
      const std::vector<CStep_Elem *> &apStepElem,
      const std::vector<int> &mapId2Ind) = 0;
 public:
  int id{};
  std::string stype;
};

class CStep_CartesianPoint : public CStep_Elem {
 public:
  void Init(
	  int id0, 
	  const std::string &stype0, 
	  const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    //
    std::string s = dfm2::Get_Parentheses(arg0, "()");
    s = dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s, ',');
    p.p[0] = stod(aToken[0]);
    p.p[1] = stod(aToken[1]);
    p.p[2] = stod(aToken[2]);
//    std::cout << p.x << " " << p.y << " " << p.z << std::endl;
  }
  void SetPtr(
	  [[maybe_unused]] const std::vector<CStep_Elem *> &apStepElem,
	  [[maybe_unused]] const std::vector<int> &mapId2Ind) override {}
  void Draw() const override {
    ::glColor3d(1, 0, 0);
    ::glPointSize(1);
    ::glBegin(GL_POINTS);
//    ::glVertex3d(p.x, p.y, p.z);
    ::glEnd();
  }
 public:
  dfm2::CVec3d p;
};

class CStep_Direction : public CStep_Elem {
 public:
  void Init(
	  int id0, 
	  const std::string &stype0, 
	  const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    //
    std::string s = dfm2::Get_Parentheses(arg0, "()");
    s = dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s, ',');
    dir.p[0] = stod(aToken[0]);
    dir.p[1] = stod(aToken[1]);
    dir.p[2] = stod(aToken[2]);
  }
  void SetPtr(
      [[maybe_unused]] const std::vector<CStep_Elem *> &apStepElem,
      [[maybe_unused]] const std::vector<int> &mapId2Ind) override {}
  void Draw() const override {}
 public:
  dfm2::CVec3d dir;
};

class CStep_Vector : public CStep_Elem {
  void Init(
	  int id0, 
	  const std::string &stype0, 
	  const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    //
    std::string s = dfm2::Get_Parentheses(arg0, "()");
    s = dfm2::RemoveSpace(s);
    std::vector<std::string> aToken = dfm2::Split(s, ',');
    std::string s1(aToken[1].begin() + 1, aToken[1].end());
    idDir = std::stoi(s1);
    len = stod(aToken[2]);
    pDir = nullptr;
  }
  void SetPtr(
      const std::vector<CStep_Elem *> &apStepElem,
      const std::vector<int> &mapId2Ind) override {
    const int iDir = mapId2Ind[idDir];
    pDir = (CStep_Direction *) apStepElem[iDir];
  }
  void Draw() const override {}
 public:
  double len{};
  int idDir{};
  ////
  const CStep_Direction *pDir{};
};

class CStep_VertexPoint : public CStep_Elem {
 public:
  void Init(
	  int id0, 
	  const std::string &stype0, 
	  const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    //
    std::string s = dfm2::Get_Parentheses(arg0, "()");
    std::vector<std::string> aToken = dfm2::Split(s, ',');
    s = dfm2::RemoveSpace(aToken[1]);
    s = std::string(s.begin() + 1, s.end());
    this->idCP = std::stoi(s);
    pCP = nullptr;
  }
  void SetPtr(
      const std::vector<CStep_Elem *> &apStepElem,
      const std::vector<int> &mapId2Ind) override {
    int icp = mapId2Ind[idCP];
    pCP = (CStep_CartesianPoint *) apStepElem[icp];
  }
  void Draw() const override {
    const dfm2::CVec3d &p0 = pCP->p;
    dfm2::opengl::myGlTranslate(+p0);
    ::glColor3d(1, 0, 1);
    dfm2::opengl::DrawSphere(16, 16);
    dfm2::opengl::myGlTranslate(-p0);
  }
 public:
  int idCP{};
  const CStep_CartesianPoint *pCP = nullptr;
};

class CStep_Axis2Placement3D : public CStep_Elem {
 public:
  void Init(
      int id0,
      const std::string &stype0,
      const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    //
    std::string s = dfm2::RemoveSpace(arg0);
    s = dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s, ',');
    assert(aToken.size() == 4);
    {
      std::string s1 = aToken[1];
      s1 = std::string(s1.begin() + 1, s1.end());
      this->idCP = std::stoi(s1);
      pCP = nullptr;
    }
    {
      std::string s2 = aToken[2];
      if (s2 == "$") { this->idDir1 = -1; }
      else {
        s2 = std::string(s2.begin() + 1, s2.end());
        this->idDir1 = std::stoi(s2);
      }
      pDir1 = nullptr;
    }
    {
      std::string s3 = aToken[3];
      if (s3 == "$") { this->idDir2 = -1; }
      else {
        s3 = std::string(s3.begin() + 1, s3.end());
        this->idDir2 = std::stoi(s3);
      }
      pDir2 = nullptr;
    }
  }
  void SetPtr(
      const std::vector<CStep_Elem *> &apStepElem,
      const std::vector<int> &mapId2Ind) override {
    const int icp = mapId2Ind[idCP];
    const int idir1 = mapId2Ind[idDir1];
    const int idir2 = mapId2Ind[idDir2];
    pCP = (CStep_CartesianPoint *) apStepElem[icp];
    pDir1 = (CStep_Direction *) apStepElem[idir1];
    pDir2 = (CStep_Direction *) apStepElem[idir2];
  }
  void Draw() const override {
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
  int idCP{}, idDir1{}, idDir2{};
  const CStep_CartesianPoint *pCP{};
  const CStep_Direction *pDir1{};
  const CStep_Direction *pDir2{};
};

class CStep_Curve : public CStep_Elem {
 public:
  [[nodiscard]] virtual double GetParameter(const dfm2::CVec3d &p) const = 0;
  virtual void SampleCurve(std::vector<dfm2::CVec3d> &polyLine,
                           double r1, double r2,
                           unsigned int nsmpl) const = 0;
};

class CStep_Line : public CStep_Curve {
  void Init(
	  int id0, 
	  const std::string &stype0, 
	  const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    //
    std::string s = dfm2::RemoveSpace(arg0);
    s = dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s, ',');
    assert(aToken.size() == 3);
    {
      std::string s1 = std::string(aToken[1].begin() + 1, aToken[1].end());
      this->idCP = std::stoi(s1);
      pCP = nullptr;
    }
    {
      std::string s2 = std::string(aToken[2].begin() + 1, aToken[2].end());
      this->idVec = std::stoi(s2);
      pVec = nullptr;
    }
  }
  void SetPtr(
      const std::vector<CStep_Elem *> &apStepElem,
      const std::vector<int> &mapId2Ind) override {
    const int iCP = mapId2Ind[idCP];
    assert(iCP != -1);
    assert(apStepElem[iCP]->stype == "CARTESIAN_POINT");
    pCP = (CStep_CartesianPoint *) apStepElem[iCP];
    ////
    const int iVec = mapId2Ind[idVec];
    assert(iVec != -1);
    assert(apStepElem[iVec]->stype == "VECTOR");
    pVec = (CStep_Vector *) apStepElem[iVec];
  }
  void Draw() const override {}
  //
  [[nodiscard]] double GetParameter(const dfm2::CVec3d &p) const override {
    const dfm2::CVec3d &cp = pCP->p;
    const dfm2::CVec3d &d = pVec->pDir->dir;
    const double l = pVec->len;
    const dfm2::CVec3d v = l * d;
    return (p - cp).dot(v) / v.squaredNorm();
  }
  void SampleCurve(
      std::vector<dfm2::CVec3d> &polyLine,
      double r1, double r2,
      unsigned int nsmpl) const override {
    const dfm2::CVec3d &cp = pCP->p;
    const dfm2::CVec3d &d = pVec->pDir->dir;
    const double l = pVec->len;
    polyLine.clear();
    for (unsigned int is = 0; is < nsmpl; ++is) {
      double r = (r2 - r1) / (nsmpl - 1) * is + r1;
      dfm2::CVec3d p = cp + r * l * d;
      polyLine.push_back(p);
    }
  }
 public:
  int idCP{};
  int idVec{};
  //
  const CStep_CartesianPoint *pCP{};
  const CStep_Vector *pVec{};
};

class CStep_Circle : public CStep_Curve {
 public:
  void Init(
	  int id0, 
	  const std::string &stype0, 
	  const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    //
    std::string s = dfm2::RemoveSpace(arg0);
    s = dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s, ',');
    assert(aToken.size() == 3);
    {
      s = std::string(aToken[1].begin() + 1, aToken[1].end());
      this->idA2P3D = std::stoi(s);
      pA2P3D = nullptr;
    }
    {
      this->radius = std::stof(aToken[2]);
    }
  }
  void SetPtr(
      const std::vector<CStep_Elem *> &apStepElem,
      const std::vector<int> &mapId2Ind) override {
    int iA2P3D = mapId2Ind[idA2P3D];
    pA2P3D = (CStep_Axis2Placement3D *) apStepElem[iA2P3D];
  }
  void Draw() const override {}
  /////
  [[nodiscard]] double GetParameter(const dfm2::CVec3d &p) const override {
    const dfm2::CVec3d &c = pA2P3D->pCP->p;
    const dfm2::CVec3d &n = pA2P3D->pDir1->dir;
    const dfm2::CVec3d ex = (pA2P3D->pDir2->dir).normalized();
    const dfm2::CVec3d ey = Cross(n, ex).normalized();
    return atan2((p - c).dot(ey), (p - c).dot(ex));
  }
  void SampleCurve(
      std::vector<dfm2::CVec3d> &polyLine,
      double r1, double r2,
      unsigned int nsmpl) const override {
    const dfm2::CVec3d &c = pA2P3D->pCP->p;
    const dfm2::CVec3d &n = pA2P3D->pDir1->dir;
    const dfm2::CVec3d ex = (pA2P3D->pDir2->dir).normalized();
    const dfm2::CVec3d ey = Cross(n, ex).normalized();
    if (r1 > r2) { r2 += M_PI * 2; }
    polyLine.clear();
    for (unsigned int is = 0; is < nsmpl; ++is) {
      double r = (r2 - r1) / (nsmpl - 1) * is + r1;
      dfm2::CVec3d p = c + radius * cos(r) * ex + radius * sin(r) * ey;
      polyLine.push_back(p);
    }
  }
 public:
  int idA2P3D{};
  double radius{};
  ////
  const CStep_Axis2Placement3D *pA2P3D{};
};

class CStep_BSplineCurveWithKnots : public CStep_Curve {
 public:
  void Init(
	  int id0, 
	  const std::string &stype0, 
	  const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    // 
    aIdCP.clear();
    aKnotMulti.clear();
    std::string s = dfm2::Get_Parentheses(arg0, "()");
    s = dfm2::RemoveSpace(s);
    std::vector<std::string> aToken = dfm2::Split_Parentheses(s, ',', "()");
    if (aToken.size() != 9) {
      std::cout << "something is wrong" << std::endl;
    }
    {
      ndegree = stoi(aToken[1]);
    }
    {
      std::string sID = aToken[2];
      sID = dfm2::Get_Parentheses(sID, "()");
      std::vector<std::string> at = dfm2::Split(sID, ',');
      for (auto & i : at) {
        const int id1 = stoi(std::string(i.begin() + 1, i.end()));
        aIdCP.push_back(id1);
      }
    }
    {
      std::string sKM = aToken[6];
      sKM = dfm2::Get_Parentheses(sKM, "()");
      std::vector<std::string> at = dfm2::Split(sKM, ',');
      for (auto & i : at) {
        const int ikm = stoi(std::string(i.begin(), i.end()));
        aKnotMulti.push_back(ikm);
      }
    }
    {
      std::string sK = aToken[7];
      sK = dfm2::Get_Parentheses(sK, "()");
      std::vector<std::string> at = dfm2::Split(sK, ',');
      for (auto & i : at) {
        const double k = stod(std::string(i.begin(), i.end()));
        aKnot.push_back(k);
      }
    }
    dfm2::FlatKnot(aKnotFlat, aKnotMulti, aKnot);
  }
  void SetPtr(const std::vector<CStep_Elem *> &apStepElem,
                      const std::vector<int> &mapId2Ind) override {
    apCP.clear();
    for (int idcp : aIdCP) {
      int icp = mapId2Ind[idcp];
      assert(apStepElem[icp]->stype == "CARTESIAN_POINT");
      apCP.push_back((CStep_CartesianPoint *) apStepElem[icp]);
    }
  }
  void Draw() const override {}
  //
  [[nodiscard]] double GetParameter([[maybe_unused]] const dfm2::CVec3d &p) const override { return 0; }
  void SampleCurve(
      std::vector<dfm2::CVec3d> &polyLine,
      [[maybe_unused]] double r1,
      [[maybe_unused]] double r2,
      unsigned int nsmpl) const override {
    std::vector<dfm2::CVec3d> aPosCP;
    for (auto iicp : apCP) {
      aPosCP.push_back(iicp->p);
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
  std::vector<CStep_CartesianPoint *> apCP;
};

class CStep_EdgeCurve : public CStep_Elem {
 public:
  void Init(
      int id0,
      const std::string &stype0,
      const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    //
    std::string s = dfm2::Get_Parentheses(arg0, "()");
    s = dfm2::RemoveSpace(s);
    std::vector<std::string> aToken = dfm2::Split(s, ',');
    assert(aToken.size() == 5);
    {
      std::string scp1 = std::string(aToken[1].begin() + 1, aToken[1].end());
      idCP1 = std::stoi(scp1);
      pCP1 = nullptr;
    }
    {
      std::string scp2 = std::string(aToken[2].begin() + 1, aToken[2].end());
      idCP2 = std::stoi(scp2);
      pCP2 = nullptr;
    }
    {
      std::string scur = std::string(aToken[3].begin() + 1, aToken[3].end());
      idCurve = std::stoi(scur);
      pCurve = nullptr;
    }
  }
  void SetPtr(
      const std::vector<CStep_Elem *> &apStepElem,
      const std::vector<int> &mapId2Ind) override {
    {
      const int iCP1 = mapId2Ind[idCP1];
      assert(iCP1 != -1);
      assert(apStepElem[iCP1]->stype == "VERTEX_POINT");
      pCP1 = (CStep_VertexPoint *) apStepElem[iCP1];
    }
    {
      const int iCP2 = mapId2Ind[idCP2];
      assert(iCP2 != -1);
      assert(apStepElem[iCP2]->stype == "VERTEX_POINT");
      pCP2 = (CStep_VertexPoint *) apStepElem[iCP2];
    }
    {
      const int iCurve = mapId2Ind[idCurve];
      if (iCurve != -1) { pCurve = (CStep_Curve *) apStepElem[iCurve]; }
      else { pCurve = nullptr; }
    }
  }
  void Draw() const override {
    ::glColor3d(0, 0, 0);
    dfm2::opengl::drawPolyLine3D(polyLine);
  }
  /////
  void Sample() {
    dfm2::CVec3d p1(pCP1->pCP->p);
    dfm2::CVec3d p2(pCP2->pCP->p);
    polyLine.clear();
    if (pCurve != nullptr) {
      double r1 = pCurve->GetParameter(p1);
      double r2 = pCurve->GetParameter(p2);
      pCurve->SampleCurve(polyLine, r1, r2, 10);
    }
  }
 public:
  int idCP1{}, idCP2{}, idCurve{};
  const CStep_VertexPoint *pCP1{};
  const CStep_VertexPoint *pCP2{};
  const CStep_Curve *pCurve{};
  ////
  std::vector<dfm2::CVec3d> polyLine;
};

class CStep_OrientedEdge : public CStep_Elem {
 public:
  void Init(
	  int id0, 
	  const std::string &stype0, 
	  const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    //
    std::string s = dfm2::RemoveSpace(arg0);
    s = dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s, ',');
    assert(aToken.size() == 5);
    {
      std::string s1(aToken[3].begin() + 1, aToken[3].end());
      idEdgeCurve = std::stoi(s1);
    }
    {
      if (aToken[4] == ".T.") { isSameDir = true; }
      else { isSameDir = false; }
    }
  }
  void SetPtr(
      const std::vector<CStep_Elem *> &apStepElem,
      const std::vector<int> &mapId2Ind) override {
    const int iEC = mapId2Ind[idEdgeCurve];
    pEC = (CStep_EdgeCurve *) apStepElem[iEC];
  }
  void Draw() const override {}
 public:
  int idEdgeCurve{};
  bool isSameDir{};
  const CStep_EdgeCurve *pEC{};
};

class CStep_EdgeLoop : public CStep_Elem {
 public:
  void Init(
      int id0,
      const std::string &stype0,
      const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    //////
    std::string s = dfm2::RemoveSpace(arg0);
    s = dfm2::Get_Parentheses(s, "()");
    s = dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s, ',');
    aIdOE.clear();
    for (unsigned int ioe = 0; ioe < aIdOE.size(); ++ioe) {
      int idOE = std::stoi(std::string(aToken[ioe].begin() + 1, aToken[ioe].end()));
      aIdOE.push_back(idOE);
    }
  }
  void SetPtr(
      const std::vector<CStep_Elem *> &apStepElem,
      const std::vector<int> &mapId2Ind) override {
    apOE.clear();
    for (int idOE : aIdOE) {
      const int iOE = mapId2Ind[idOE];
      assert(iOE != -1);
      assert(apStepElem[iOE]->stype == "ORIENTED_EDGE");
      apOE.push_back((CStep_OrientedEdge *) apStepElem[iOE]);
    }
  }
  void Draw() const override {}
 public:
  std::vector<int> aIdOE;
  std::vector<CStep_OrientedEdge *> apOE;
};

class CStep_FaceOuterBound : public CStep_Elem {
 public:
  void Init(
	  int id0, 
	  const std::string &stype0, 
	  const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    //
    std::string s = dfm2::RemoveSpace(arg0);
    s = dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s, ',');
    {
      std::string s1 = std::string(aToken[1].begin() + 1, aToken[1].end());
      idEL = std::stoi(s1);
      pEL = nullptr;
    }
    {
      std::string s1 = aToken[2];
      if (s1 == ".T.") { orientation = true; }
      else { orientation = false; }
    }
  }
  void SetPtr(
      const std::vector<CStep_Elem *> &apStepElem,
      const std::vector<int> &mapId2Ind) override {
    const int iEL = mapId2Ind[idEL];
    assert(iEL != -1);
    assert(apStepElem[iEL]->stype == "EDGE_LOOP");
    pEL = (CStep_EdgeLoop *) apStepElem[iEL];
  }
  void Draw() const override {}
 public:
  int idEL{};
  bool orientation{};
  ///
  const CStep_EdgeLoop *pEL{};
};

class CStep_Surface : public CStep_Elem {
  virtual void Project(double &r0, double &r1,
                       const dfm2::CVec3d &p0) const = 0;
 public:
};

class CStep_Plane : public CStep_Surface {
 public:
  void Init(
	  int id0, 
	  const std::string &stype0, 
	  const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    //
    std::string s = dfm2::RemoveSpace(arg0);
    s = dfm2::Get_Parentheses(s, "()");
    std::vector<std::string> aToken = dfm2::Split(s, ',');
    assert(aToken.size() == 2);
    {
      std::string s1 = std::string(aToken[1].begin() + 1, aToken[1].end());
      idA2P3D = std::stoi(s1);
      pA2P3D = nullptr;
    }
  }
  void SetPtr(const std::vector<CStep_Elem *> &apStepElem,
                      const std::vector<int> &mapId2Ind) override {
    const int iA2P3D = mapId2Ind[idA2P3D];
    assert(iA2P3D != -1);
    assert(apStepElem[iA2P3D]->stype == "AXIS2_PLACEMENT_3D");
    pA2P3D = (CStep_Axis2Placement3D *) apStepElem[iA2P3D];
  }
  void Draw() const override {
    dfm2::CVec3d org = pA2P3D->pCP->p;
    dfm2::CVec3d n = pA2P3D->pDir1->dir;
    dfm2::CVec3d ex = pA2P3D->pDir2->dir;
    dfm2::CVec3d ey = Cross(n, ex);
    double l = 100;
//    ::glEnable(GL_LIGHTING);
    ::glBegin(GL_LINE_LOOP);
    dfm2::opengl::myGlNormal(n);
    dfm2::opengl::myGlVertex(org - l * ex - l * ey);
    dfm2::opengl::myGlVertex(org + l * ex - l * ey);
    dfm2::opengl::myGlVertex(org + l * ex + l * ey);
    dfm2::opengl::myGlVertex(org - l * ex + l * ey);
    ::glEnd();
  }
  void Project(
      double &r0, double &r1,
      const dfm2::CVec3d &p0) const override {
    const dfm2::CVec3d &org = pA2P3D->pCP->p;
    const dfm2::CVec3d &n = pA2P3D->pDir1->dir;
    const dfm2::CVec3d &ex = pA2P3D->pDir2->dir;
    const dfm2::CVec3d &ey = Cross(n, ex);
    r0 = (p0 - org).dot(ex);
    r1 = (p0 - org).dot(ey);
  }
 public:
  int idA2P3D{};
  const CStep_Axis2Placement3D *pA2P3D{};
};

class CStep_AdvancedFace : public CStep_Elem {
 public:
  void Init(
	  int id0, 
	  const std::string &stype0,
	  const std::string &arg0) override {
    this->id = id0;
    this->stype = stype0;
    //
    std::string s = dfm2::RemoveSpace(arg0);
    std::vector<std::string> aToken;
    {
      std::string s0 = dfm2::Get_Parentheses(s, "()");
      aToken = dfm2::Split_Parentheses(s0, ',', "()");
      for (unsigned int it = 0; it < aToken.size(); ++it) {
        std::cout << it << " " << aToken[it] << "   " << arg0 << " " << s << std::endl;
      }
    }
    assert(aToken.size() == 4);
    {
      std::string s1 = std::string(aToken[2].begin() + 1, aToken[2].end());
      idSurf = std::stoi(s1);
      pSurf = nullptr;
    }
    {
      s = dfm2::Get_Parentheses(aToken[1], "()");
      aToken = dfm2::Split(s, ',');
      aIdFOB.clear();
      for (auto & it : aToken) {
        s = std::string(it.begin() + 1, it.end());
        int idFOB = std::stoi(s);
        aIdFOB.push_back(idFOB);
      }
    }
  }
  void SetPtr(
      const std::vector<CStep_Elem *> &apStepElem,
      const std::vector<int> &mapId2Ind) override {
    apFOB.clear();
    for (int idFOB : aIdFOB) {
      const int iFOB = mapId2Ind[idFOB];
      if (iFOB == -1) {
        apFOB.push_back(nullptr);
        continue;
      }
      assert(iFOB != -1);
      assert(apStepElem[iFOB]->stype == "FACE_OUTER_BOUND");
      apFOB.push_back((CStep_FaceOuterBound *) apStepElem[iFOB]);
    }
  }
  void Draw() const override {}
  ///
  void SampleSurface() {
  }
 public:
  std::vector<int> aIdFOB;
  int idSurf{};
  CStep_Surface *pSurf{};
  std::vector<CStep_FaceOuterBound *> apFOB;
  /////
  std::vector<double> aXYZ;
  std::vector<int> aTri;
};

void LoadStep
    (const std::string &fname,
     std::vector<CStep_Elem *> &apStepElem) {
  std::cout << "load Step" << std::endl;
  {
    apStepElem.clear();
    std::ifstream fin(fname.c_str());
    std::string s;
    s.reserve(1024 * 8);
    while (getline(fin, s, ';')) {
      std::string type, arg;
      int id = getEntity(type, arg, s);
      std::cout << id << " " << type << " " << arg << " " << s << std::endl;
      CStep_Elem *pSE = 0;
      if (type == "CARTESIAN_POINT") { pSE = new CStep_CartesianPoint(); }
      if (type == "DIRECTION") { pSE = new CStep_Direction(); }
      if (type == "VERTEX_POINT") { pSE = new CStep_VertexPoint(); }
      if (type == "VECTOR") { pSE = new CStep_Vector(); }
      if (type == "AXIS2_PLACEMENT_3D") { pSE = new CStep_Axis2Placement3D(); }
      if (type == "LINE") { pSE = new CStep_Line(); }
      if (type == "CIRCLE") { pSE = new CStep_Circle(); }
      if (type == "B_SPLINE_CURVE_WITH_KNOTS") { pSE = new CStep_BSplineCurveWithKnots(); }
      if (type == "EDGE_CURVE") { pSE = new CStep_EdgeCurve(); }
      if (type == "ORIENTED_EDGE") { pSE = new CStep_OrientedEdge(); }
      if (type == "EDGE_LOOP") { pSE = new CStep_EdgeLoop(); }
      if (type == "FACE_OUTER_BOUND") { pSE = new CStep_FaceOuterBound(); }
      if (type == "PLANE") { pSE = new CStep_Plane(); }
      if (type == "ADVANCED_FACE") { pSE = new CStep_AdvancedFace(); }
      if (pSE != nullptr) {
        pSE->Init(id, type, arg);
        apStepElem.push_back(pSE);
      }
    }
  }
  // ---------------------------------------
  { // id 2 ind
    std::vector<int> mapId2Ind;
    for (unsigned int ise = 0; ise < apStepElem.size(); ++ise) {
      int id = apStepElem[ise]->id;
      if (id >= (int) mapId2Ind.size()) {
        mapId2Ind.resize(id + 1, -1);
      }
      mapId2Ind[id] = ise;
    }
    for (unsigned int ise = 0; ise < apStepElem.size(); ++ise) {
      apStepElem[ise]->SetPtr(apStepElem, mapId2Ind);
    }
    for (auto & ise : apStepElem) {
      if (ise->stype == "EDGE_CURVE") {
        ((CStep_EdgeCurve *) ise)->Sample();
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

int main() {
  std::vector<CStep_Elem *> apStepElem;
  //  LoadStep("Hook.STEP",aCP,aBSCWK);
  LoadStep(std::string(PATH_INPUT_DIR) + "/bolt.STEP", apStepElem);
  //  LoadStep("grabcad/bumper/BUMPER.stp",aCP,aBSCWK);
  //  LoadStep("bolt.STEP",aCP,aBSCWK);

  delfem2::glfw::CViewer3 viewer(100);
  //
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  //
  while (!glfwWindowShouldClose(viewer.window)) {
    viewer.DrawBegin_oldGL();
    for (auto & ipse : apStepElem) {
      ipse->Draw();
    }
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
