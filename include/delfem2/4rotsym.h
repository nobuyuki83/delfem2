//
// Created by Nobuyuki Umetani on 2020-09-27.
//

/**
 * @brief implementation of 4 rotatoinal symetry field
 * @details implementation is based on the paper
 * "Wenzel Jakob, Marco Tarini, Daniele Panozzo, and Olga Sorkine-Hornung.
 * Instant field-aligned meshes. Siggraph Asia 2015"
 */

#ifndef DFM2_4ROTSYM_H
#define DFM2_4ROTSYM_H

#include "delfem2/vec3.h"

namespace delfem2
{

void FindNearestOrientation(
    CVec3d& d0,
    CVec3d& d1,
    const CVec3d& o0,
    const CVec3d& p0,
    const CVec3d& o1,
    const CVec3d& p1)
{
  const CVec3d ad0[4] = { o0, p0, -o0, -p0 };
  const CVec3d ad1[2] = { o1, p1 };
  double dot_max = -2;
  for(const auto & j : ad1){
    for(const auto & i : ad0){
      const double dot = i*j;
      if( dot < dot_max ){ continue; }
      d0 = i;
      d1 = j;
      dot_max = dot;
    }
  }
}

void Smooth4RotSym
    (std::vector<double>& aOdir,
     const std::vector<double>& aNorm,
     const std::vector<unsigned int>& psup_ind,
     const std::vector<unsigned int>& psup)
{
  for(unsigned int iip=0;iip<aOdir.size()/3;++iip){
    const unsigned int ip0 = iip;
    assert( ip0 < psup_ind.size() );
    const unsigned int npj = psup_ind[ip0+1] - psup_ind[ip0+0];
    const CVec3d n0 = CVec3d(aNorm.data()+ip0*3);
    CVec3d o_new = CVec3d(aOdir.data()+ip0*3);
    double weight = 0.0;
    for(unsigned int jjp=0;jjp<npj;++jjp){
      unsigned int jp1 = psup[psup_ind[ip0]+jjp];
      const CVec3d n1 = CVec3d(aNorm.data()+jp1*3);
      const CVec3d o1 = CVec3d(aOdir.data()+jp1*3);
      CVec3d d0, d1;
      FindNearestOrientation(d0,d1,
                             o_new,n0^o_new, o1,n1^o1);
      o_new = d0 * weight + d1;
      o_new = (o_new - (o_new*n0)*n0).Normalize();
      weight += 1.0;
    }
    o_new.CopyTo(aOdir.data()+ip0*3);
  }
}

void Smooth4RotSym_RandomPermutation(
    std::vector<double>& aOdir,
    std::vector<unsigned int>& permutation0,
    const std::vector<double>& aNorm,
    const std::vector<unsigned int>& psup_ind,
    const std::vector<unsigned int>& psup)
{
  std::vector<unsigned int> permutation1(10);
  std::random_shuffle( permutation0.begin(), permutation0.end() );
  for(int iip=0;iip<aOdir.size()/3;++iip){
    const unsigned int ip0 = permutation0[iip];
    assert( ip0 < psup_ind.size() );
    const unsigned int npj = psup_ind[ip0+1] - psup_ind[ip0+0];
    const CVec3d n0 = CVec3d(aNorm.data()+ip0*3);
    CVec3d o_new = CVec3d(aOdir.data()+ip0*3);
    double weight = 0.0;
    {
      permutation1.resize(npj);
      for(int jjp=0;jjp<npj;++jjp){ permutation1[jjp] = jjp; }
      std::random_shuffle( permutation1.begin(), permutation1.end() );
    }
    for(int jjp=0;jjp<npj;++jjp){
      unsigned int jp1 = psup[psup_ind[ip0]+permutation1[jjp]];
      const CVec3d n1 = CVec3d(aNorm.data()+jp1*3);
      const CVec3d o1 = CVec3d(aOdir.data()+jp1*3);
      CVec3d d0, d1;
      FindNearestOrientation(d0,d1,
                             o_new,n0^o_new, o1,n1^o1);
      o_new = d0 * weight + d1;
      o_new = (o_new - (o_new*n0)*n0).Normalize();
      weight += 1.0;
    }
    o_new.CopyTo(aOdir.data()+ip0*3);
  }
}

}


#endif //EXAMPLES_GLFWOLD_HDRONLY_4ROTSYM_H
