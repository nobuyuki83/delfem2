/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEMCLOTH_H
#define DFM2_FEMCLOTH_H

#include <vector>

#include "delfem2/dfm2_inline.h"
#include "delfem2/femutil.h"

#ifdef DFM2_STATIC_LIBRARY
// Merge use explicitly use the template so for static library we need to include the template itself.
#  include "delfem2/lsmats.h"
#endif

namespace delfem2 {

/**
 * @brief compute energy and its 1st and 2nd derivative for cloth bending
 * @param[out] W strain energy
 * @param[out] dW 1st derivative of energy
 * @param[out] ddW 2nd derivative of energy
 * @param[in] C undeformed quad vertex positions
 * @param[in] c deformed quad vertex positions
 */
void WdWddW_Bend(
    double& W,
    double dW[4][3],
    double ddW[4][4][3][3],
    //
    const double C[4][3],
    const double c[4][3],
    double stiff);

void WdWddW_CST(
    double& W, // (out) energy
    double dW[3][3], // (out) 1st derivative of energy
    double ddW[3][3][3][3], // (out) 2nd derivative of energy
    //
    const double C[3][3], // (in) undeformed triangle vertex positions
    const double c[3][3], // (in) deformed triangle vertex positions
    const double lambda, // (in) Lame's 1st parameter
    const double myu);   // (in) Lame's 2nd parameter

/**
 * @brief compute energy and its 1st and 2nd derivative for contact against object
 */
class CInput_Contact
{
public:
  virtual double penetrationNormal(
      double& nx, double& ny, double& nz,
      double px, double py, double pz) const = 0;
};

void WdWddW_Contact(
    double& W,  // (out) energy
    double dW[3], // (out) 1st derivative of energy
    double ddW[3][3], // (out) 2nd derivative of energy
    //
    const double c[3], // (in) deformed triangle vertex positions
    double stiff_contact,
    double contact_clearance,
    const CInput_Contact& input);


// compute total energy and its first and second derivatives
template <class MAT>
double MergeLinSys_Cloth(
    MAT& ddW, // (out) second derivative of energy
    double* dW, // (out) first derivative of energy
    //
    double lambda, // (in) Lame's 1st parameter
    double myu,  // (in) Lame's 2nd parameter
    double stiff_bend, // (in) bending stiffness
    const double* aPosIni,
    unsigned int np,
    unsigned int ndim,
    const unsigned int* aTri,
    unsigned int nTri, // (in) triangle index
    const unsigned int* aQuad,
    unsigned int nQuad, // (in) index of 4 vertices required for bending
    const double* aXYZ)
{
  assert( ndim == 2 || ndim == 3 );
  double W = 0;
  std::vector<unsigned int> tmp_buffer(np,UINT_MAX);
  // marge element in-plane strain energy
  for(unsigned int itri=0;itri<nTri;itri++){
    const unsigned int aIP[3] = { aTri[itri*3+0], aTri[itri*3+1], aTri[itri*3+2] };
    double C[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    double c[3][3];
    for(int ino=0;ino<3;ino++){
      const unsigned int ip = aIP[ino];
      for(unsigned int i=0;i<ndim;i++){ C[ino][i] = aPosIni[ip*ndim+i]; }
      for(int i=0;i<3;i++){ c[ino][i] = aXYZ[ip*3+i]; }
    }
    double e, de[3][3], dde[3][3][3][3];
    WdWddW_CST( e,de,dde, C,c, lambda,myu );
    W += e;  // marge energy
    // marge de
    for(int ino=0;ino<3;ino++){
      const unsigned int ip = aIP[ino];
      for(int i =0;i<3;i++){ dW[ip*3+i] += de[ino][i]; }
    }
    // marge dde
//    ddW.Mearge(3, aIP, 3, aIP, 9, &dde[0][0][0][0], tmp_buffer);
    Merge<3,3,3,3,double>(ddW,aIP,aIP,dde,tmp_buffer);
  }
//  std::cout << "cst:" << W << std::endl;
  // marge element bending energy
  for(unsigned int iq=0;iq<nQuad;iq++){
    const unsigned int aIP[4] = { aQuad[iq*4+0], aQuad[iq*4+1], aQuad[iq*4+2], aQuad[iq*4+3] };
    double C[4][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
    double c[4][3];
    for(int ino=0;ino<4;ino++){
      const unsigned int ip = aIP[ino];
      for(unsigned int i=0;i<ndim;i++){ C[ino][i] = aPosIni[ip*ndim+i]; }
      for(int i=0;i<3;i++){ c[ino][i] = aXYZ [ip*3+i]; }
    }
    double e, de[4][3], dde[4][4][3][3];
    WdWddW_Bend( e,de,dde, C,c, stiff_bend );
    W += e;  // marge energy
    // marge de
    for(int ino=0;ino<4;ino++){
      const unsigned int ip = aIP[ino];
      for(int i =0;i<3;i++){ dW[ip*3+i] += de[ino][i]; }
    }
    // marge dde
//    ddW.Mearge(4, aIP, 4, aIP, 9, &dde[0][0][0][0], tmp_buffer);
    Merge<4,4,3,3,double>(ddW,aIP,aIP,dde,tmp_buffer);
  }
  return W;
}

template <class MAT>
double MergeLinSys_Contact(
    MAT& ddW,
    double* dW,
    //
    double stiff_contact,
    double contact_clearance,
    const CInput_Contact& input,
    const double* aXYZ,
    int nXYZ)
{
  const unsigned int np = nXYZ;
  std::vector<unsigned int> tmp_buffer(np,UINT_MAX);
  double W = 0;
  for(unsigned int ip=0;ip<np;ip++){
    double c[3] = { aXYZ[ip*3+0], aXYZ[ip*3+1], aXYZ[ip*3+2] };
    double e, de[3], dde[1][1][3][3];
    WdWddW_Contact( e,de,dde[0][0], c, stiff_contact,contact_clearance, input );
    W += e;  // marge energy
    for(int i =0;i<3;i++){ dW[ip*3+i] += de[i]; }     // marge de
//    ddW.Mearge(1, &ip, 1, &ip, 9, &dde[0][0], tmp_buffer);
    Merge<1,1,3,3,double>(ddW,&ip,&ip,dde,tmp_buffer);     // marge dde
  }
  return W;
}



} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/femcloth.cpp"
#endif
  
#endif /* fem_ematrix_h */
