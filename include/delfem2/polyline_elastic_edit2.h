
#ifndef DFM2_POLYLINE_ELASTIC_EDIT_H
#define DFM2_POLYLINE_ELASTIC_EDIT_H

#include <cmath>
#include <vector>
#include <climits>

#include "delfem2/ls_tridiagonal.h"

class PolylineElasticEdit2 {
 public:
  PolylineElasticEdit2();
  ~PolylineElasticEdit2();

  template<class ARRAY2D>
  void Initialize(
      const ARRAY2D &v2s,
      unsigned int ip_drag,
      double X,
      double Y);

  template<class ARRAY2D>
  void Drag(
      ARRAY2D &v2s,
      double x,
      double y);

 private:
  void SolveLinearStatic();

  void ClearMemory();

  [[nodiscard]] unsigned int GetSizeNode() const { return ut.size() / 3; }

  void ProjectPoint(
      double x_in, double y_in, int &idiv_min,
      double &alpha, double &ndist, double &norm_x, double &norm_y);

 private:

  static double FindNearestPointParam_Line_Point(
      const double pc[2],
      const double ps[2],
      const double pe[2]) {
    const double es[2] = {pe[0] - ps[0], pe[1] - ps[1]};
    const double sc[2] = {ps[0] - pc[0], ps[1] - pc[1]};
    const double a = es[0] * es[0] + es[1] * es[1];
    const double b = es[0] * sc[0] + es[1] * sc[1];
    return -b / a;
  }

 public:
  double X0 = 0.;
  double Y0 = 0.;
  unsigned int ip_drag_ = UINT_MAX;
  // --------------
  std::vector<double> ut;   // displacement [nno*3]
  std::vector<double> ini_x; // initial position
  std::vector<int> bc_flag;
  // --------------
  double EI;
  double ARho;
  double AE;
  double rho;
 public:
  std::vector<double> dut;
  std::vector<double> Res;
  BlockTriDiagonalMatrix m_mat;
};

PolylineElasticEdit2::PolylineElasticEdit2() {
  this->EI = 30.0;
  this->ARho = 1;
  this->AE = 100000.0;
}

template<class ARRAY2D>
void PolylineElasticEdit2::Initialize(
    const ARRAY2D &v2s,
    unsigned int ip_drag,
    double X,
    double Y) {
  this->ClearMemory();
  this->X0 = X;
  this->Y0 = Y;
  this->ip_drag_ = ip_drag;
  const size_t nno = v2s.size();
  this->m_mat.Initialize(nno);
  this->ini_x.resize(nno * 2);
  this->ut.resize(nno * 3);
  this->bc_flag.resize(nno * 3);
  this->Res.resize(nno * 3);
  this->dut.resize(nno * 3);
  // ------------
  for (unsigned int ino = 0; ino < nno; ino++) {
    ini_x[ino * 2 + 0] = v2s[ino][0];
    ini_x[ino * 2 + 1] = v2s[ino][1];
  }
  ut.assign(nno * 3, 0.);
  bc_flag.assign(nno * 3, 0);
  // ----------------------------
  bc_flag[0 * 3 + 0] = 1;
  bc_flag[0 * 3 + 1] = 1;
  bc_flag[0 * 3 + 2] = 1;
//  bc_flag[(nno-1)*3+0] = 1;
//  bc_flag[(nno-1)*3+1] = 1;
//  bc_flag[(nno-1)*3+2] = 1;
  bc_flag[ip_drag * 3 + 0] = 1;
  bc_flag[ip_drag * 3 + 1] = 1;
//  bc_flag[ip_drag * 3 + 2] = 1;
}

template<class ARRAY2D>
void PolylineElasticEdit2::Drag(
    ARRAY2D &v2s,
    double x,
    double y) {
  const size_t nno = ut.size() / 3;
  v2s.resize(nno);
  ut[ip_drag_ * 3 + 0] = (x - X0);//-ini_x[ip_drag*2+0];
  ut[ip_drag_ * 3 + 1] = (y - Y0);//-ini_x[ip_drag*2+1];
  this->SolveLinearStatic();
  for (unsigned int ino = 0; ino < nno; ino++) {
    v2s[ino][0] = ut[ino * 3 + 0] + ini_x[ino * 2 + 0];
    v2s[ino][1] = ut[ino * 3 + 1] + ini_x[ino * 2 + 1];
  }
  /*
  ut[ip_drag_*3+0] = x-ini_x[ip_drag_*2+0];
  ut[ip_drag_*3+1] = y-ini_x[ip_drag_*2+1];
  this->SolveLinearStatic();
  for(unsigned int ino=0;ino<nno;ino++){
    v2s[ino][0] = ut[ino*3+0]+ini_x[ino*2+0];
    v2s[ino][1] = ut[ino*3+1]+ini_x[ino*2+1];
    ini_x[ino*2+0] += ut[ino*3+0];
    ini_x[ino*2+1] += ut[ino*3+1];
    ut[ino*3+0] = 0;
    ut[ino*3+1] = 0;
    ut[ino*3+2] = 0;
  }
   */
}

void PolylineElasticEdit2::ClearMemory() {
  ini_x.clear();
  ut.clear();
  bc_flag.clear();
  Res.clear();
  dut.clear();
  m_mat.Clear();
}

PolylineElasticEdit2::~PolylineElasticEdit2() {
  this->ClearMemory();
}

void PolylineElasticEdit2::ProjectPoint(
    double x_in, double y_in, int &idiv_min,
    double &alpha, double &ndist, double &norm_x, double &norm_y) {
  const size_t nno = ut.size() / 3;
  const double x0[2] = {x_in, y_in};
  double dist = std::sqrt(
      (ini_x[0] + ut[0] - x0[0]) * (ini_x[0] + ut[0] - x0[0]) +
          (ini_x[1] + ut[1] - x0[1]) * (ini_x[1] + ut[1] - x0[1]));
  idiv_min = -1;
  const unsigned int ndiv = nno - 1;
  for (unsigned int idiv = 0; idiv < ndiv; idiv++) {
    double x1[2] = {ini_x[idiv * 2 + 0] + ut[idiv * 3 + 0], ini_x[idiv * 2 + 1] + ut[idiv * 3 + 1]};
    double x2[2] = {ini_x[idiv * 2 + 2] + ut[idiv * 3 + 3], ini_x[idiv * 2 + 3] + ut[idiv * 3 + 4]};
    double t = FindNearestPointParam_Line_Point(x0, x1, x2);
    if (t < -0.001 || t > 1.001) { continue; }
    double x3[2] = {x1[0] * (1 - t) + x2[0] * t, x1[1] * (1 - t) + x2[1] * t};
    double d[2] = {x0[0] - x3[0], x0[1] - x3[1]};
    double elen = std::sqrt((x1[0] - x2[0]) * (x1[0] - x2[0]) + (x1[1] - x2[1]) * (x1[1] - x2[1]));
    double n[2] = {(x2[1] - x1[1]) / elen, (x1[0] - x2[0]) / elen};
    double dlen = std::sqrt(d[0] * d[0] + d[1] * d[1]);
    if (dlen < dist) {
      idiv_min = idiv;
      dist = dlen;
      alpha = t;
      ndist = (n[0] * d[0] + n[1] * d[1]) / elen;
      norm_x = n[0];
      norm_y = n[1];
    }
  }
}

void dWddW_Beam2(
    double dW[2][3],
    double ddW[2][2][3][3],
    //
    double EI,
    double AE,
    const double x[2][2],
    const double u[2][3]) {
  const double eLen = sqrt((x[1][0] - x[0][0]) * (x[1][0] - x[0][0]) + (x[1][1] - x[0][1]) * (x[1][1] - x[0][1]));

  double eC[2][2][3][3];
  {
    std::fill_n(&eC[0][0][0][0], 36, 0.0);
    const double tmp1 = EI / (eLen * eLen * eLen);  // bending stiffness
    const double tmp2 = AE / eLen;  // axal stiffness x direction
    /*
    eC[0][0] = eC[3][3] =  tmp2;
    eC[3][0] = eC[0][3] = -tmp2;
    eC[1][1] = eC[4][4] =  tmp1*12;
    eC[4][1] = eC[1][4] = -tmp1*12;
    eC[1][2] = eC[2][1] =  eC[5][1] = eC[1][5] =  tmp1*eLen*6;
    eC[2][4] = eC[4][2] =  eC[4][5] = eC[5][4] = -tmp1*eLen*6;
    eC[2][2] = eC[5][5] =  tmp1*eLen*eLen*4;
    eC[5][2] = eC[2][5] =  tmp1*eLen*eLen*2;
     */
    eC[0][0][0][0] = eC[1][1][0][0] = tmp2;
    eC[1][0][0][0] = eC[0][1][0][0] = -tmp2;
    eC[0][0][1][1] = eC[1][1][1][1] = tmp1 * 12;
    eC[1][0][1][1] = eC[0][1][1][1] = -tmp1 * 12;
    eC[0][0][1][2] = eC[0][0][2][1] = eC[1][0][2][1] = eC[0][1][1][2] = tmp1 * eLen * 6;
    eC[0][1][2][1] = eC[1][0][1][2] = eC[1][1][1][2] = eC[1][1][2][1] = -tmp1 * eLen * 6;
    eC[0][0][2][2] = eC[1][1][2][2] = tmp1 * eLen * eLen * 4;
    eC[1][0][2][2] = eC[0][1][2][2] = tmp1 * eLen * eLen * 2;
  }

  const double inv_eLen = 1.0 / eLen;
  const double T[2] = {
      (x[1][0] - x[0][0]) * inv_eLen,
      (x[1][1] - x[0][1]) * inv_eLen};  // cosine and sine
  // rotation for the of the stiffness matrix
  const double eR[3][3] = {
      {T[0], -T[1], 0},
      {T[1], T[0], 0},
      {0, 0, 1} };
  std::fill_n(&ddW[0][0][0][0], 36, 0.0);
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int n = 0; n < 3; n++) {
        for (unsigned int m = 0; m < 3; m++) {
          ddW[0][0][i][j] += eR[i][n] * eC[0][0][n][m] * eR[j][m];
          ddW[0][1][i][j] += eR[i][n] * eC[0][1][n][m] * eR[j][m];
          ddW[1][0][i][j] += eR[i][n] * eC[1][0][n][m] * eR[j][m];
          ddW[1][1][i][j] += eR[i][n] * eC[1][1][n][m] * eR[j][m];
        }
      }
    }
  }
  for (unsigned int n = 0; n < 2; n++) {
    for (unsigned int i = 0; i < 3; i++) {
      dW[n][i]
          = ddW[n][0][i][0] * u[0][0] + ddW[n][0][i][1] * u[0][1] + ddW[n][0][i][2] * u[0][2]
          + ddW[n][1][i][0] * u[1][0] + ddW[n][1][i][1] * u[1][1] + ddW[n][1][i][2] * u[1][2];
    }
  }
}

void PolylineElasticEdit2::SolveLinearStatic() {
  const size_t nno = ut.size() / 3;
  for (unsigned int i = 0; i < nno * 3; i++) { Res[i] = 0; }
  m_mat.SetZero();
  //
  const unsigned int ndiv = nno - 1;
  for (unsigned int idiv = 0; idiv < ndiv; idiv++) {
    const double x1[2][2] = {
        {ini_x[idiv * 2 + 0], ini_x[idiv * 2 + 1]},
        {ini_x[idiv * 2 + 2], ini_x[idiv * 2 + 3]}};
    const double u1[2][3] = {
        {ut[idiv * 3 + 0], ut[idiv * 3 + 1], ut[idiv * 3 + 2]},
        {ut[idiv * 3 + 3], ut[idiv * 3 + 4], ut[idiv * 3 + 5]}};
    double ddW[2][2][3][3], dW[2][3];
    dWddW_Beam2(
        dW, ddW,
        EI, AE,
        x1, u1);

    /*
    eC[0][0][0][0] += 10000;
    eC[0][0][1][1] += 10000;
    eC[1][1][0][0] += 10000;
    eC[1][1][1][1] += 10000;
     */
    m_mat.Merge(idiv, ddW);
    for (unsigned int i = 0; i < 3; i++) {
      Res[idiv * 3 + i] -= dW[0][i];
      Res[idiv * 3 + 3 + i] -= dW[1][i];
    }
  }
  for (unsigned int ino = 0; ino < nno; ino++) {
    for (unsigned int idof = 0; idof < 3; idof++) {
      if (bc_flag[ino * 3 + idof] == 0) continue;
      m_mat.FixBC(ino, idof);
      Res[ino * 3 + idof] = 0;
    }
  }
  {
    double nres = 0;
    for (unsigned int i = 0; i < nno * 3; i++) { nres += Res[i] * Res[i]; }
  }
  m_mat.Decompose();
  {
    for (unsigned int i = 0; i < nno * 3; i++) { dut[i] = Res[i]; }
    m_mat.Solve(dut);
  }
  for (unsigned int ino = 0; ino < nno; ino++) {
    ut[ino * 3 + 0] += dut[ino * 3 + 0];
    ut[ino * 3 + 1] += dut[ino * 3 + 1];
    ut[ino * 3 + 2] += dut[ino * 3 + 2];
  }
}

#endif 
