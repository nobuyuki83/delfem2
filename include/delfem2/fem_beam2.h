//
// Created by Nobuyuki Umetani on 2021-11-07.
//

#ifndef FEM_BEAM2_H_
#define FEM_BEAM2_H_

namespace delfem2 {

template <typename T>
void dWddW_Beam2(
    T dW[2][3],
    T ddW[2][2][3][3],
    //
    T EI,
    T AE,
    const T x[2][2],
    const T u[2][3]) {
  const T eLen = std::sqrt((x[1][0] - x[0][0]) * (x[1][0] - x[0][0]) + (x[1][1] - x[0][1]) * (x[1][1] - x[0][1]));

  T eC[2][2][3][3];
  {
    std::fill_n(&eC[0][0][0][0], 36, 0.0);
    const T tmp1 = EI / (eLen * eLen * eLen);  // bending stiffness
    const T tmp2 = AE / eLen;  // axal stiffness x direction
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

  const T inv_eLen = 1.0 / eLen;
  const T cs[2] = {
      (x[1][0] - x[0][0]) * inv_eLen,
      (x[1][1] - x[0][1]) * inv_eLen};  // cosine and sine
  // rotation for the of the stiffness matrix
  const T eR[3][3] = {
      {cs[0], -cs[1], 0},
      {cs[1], cs[0], 0},
      {0, 0, 1}};

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

}  // namespace delfem2


#endif //FEM_BEAM2_H_
