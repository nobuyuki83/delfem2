//
// Created by Nobuyuki Umetani on 2020-09-27.
//

/**
 * @brief implementation of 4 rotatoinal symetry field
 * @details implementation is based on the paper:
 * "Wenzel Jakob, Marco Tarini, Daniele Panozzo, and Olga Sorkine-Hornung.
 * Instant field-aligned meshes. Siggraph Asia 2015"
 */

#ifndef DFM2_4ROTSYM_H
#define DFM2_4ROTSYM_H

#include <climits>
#include <algorithm> // shuffle
#include <random>
#include "delfem2/vec3.h"

namespace delfem2
{

CVec3d FindNearestOrientation(
    unsigned int& idiff,
    const CVec3d& o0,
    const CVec3d& n0,
    const CVec3d& o1,
    const CVec3d& n1,
    double weight)
{
  const auto p0 =n0^o0;
  const auto p1 =n1^o1;
  const CVec3d ad0[2] = { o0, p0 };
  const CVec3d ad1[4] = { o1, p1, -o1, -p1 };
  double dot_max = -2;
  unsigned int id0_near = UINT_MAX;
  unsigned int id1_near = UINT_MAX;
  for(unsigned int id0=0;id0<2;++id0){
    for(unsigned int id1=0;id1<4;++id1){
      const double dot = ad0[id0].dot(ad1[id1]);
      if( dot < dot_max ){ continue; }
      id0_near = id0;
      id1_near = id1;
      dot_max = dot;
    }
  }
  idiff = (id0_near+4-id1_near)%4;
  CVec3d d0new = ad0[id0_near] * weight + ad1[id1_near];
  d0new = (d0new - (d0new.dot(n0)) * n0).normalized();
  if( id0_near == 1 ){ d0new = d0new ^ n0; }
  return d0new;
}

void Smooth4RotSym(
    std::vector<double>& aOdir,
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
      unsigned int idiff;
      o_new = FindNearestOrientation(idiff,
          o_new,n0, o1,n1, weight);
      weight += 1.0;
    }
    o_new.CopyTo(aOdir.data()+ip0*3);
  }
}

void Smooth4RotSym2(
    double* aOdir,
    unsigned int np,
    const unsigned int* psup_ind,
    const unsigned int* psup)
{
  for(unsigned int ip=0;ip<np;++ip){
    const unsigned int npj = psup_ind[ip+1] - psup_ind[ip+0];
    const CVec3d N = CVec3d(0,0,1);
    CVec3d o_new = CVec3d(aOdir[ip*2+0],aOdir[ip*2+1],0);
    double weight = 1.0;
    for(unsigned int jjp=0;jjp<npj;++jjp){
      unsigned int jp1 = psup[psup_ind[ip]+jjp];
      const CVec3d o1 = CVec3d(aOdir[jp1*2+0],aOdir[jp1*2+1],0);
      unsigned int idiff;
      o_new = FindNearestOrientation(idiff,
          o_new,N, o1, N, weight);
      weight += 1.0;
    }
    aOdir[ip*2+0] = o_new.x;
    aOdir[ip*2+1] = o_new.y;
  }
}

void Smooth4RotSym_RandomPermutation(
    std::vector<double>& aOdir,
    std::vector<unsigned int>& permutation0,
    const std::vector<double>& aNorm,
    const std::vector<unsigned int>& psup_ind,
    const std::vector<unsigned int>& psup)
{
  std::random_device rd;
  std::mt19937 g(rd());
  std::vector<unsigned int> permutation1(10);
  std::shuffle( permutation0.begin(), permutation0.end(), g);
  for(unsigned int iip=0;iip<aOdir.size()/3;++iip){
    const unsigned int ip0 = permutation0[iip];
    assert( ip0 < psup_ind.size() );
    const unsigned int npj = psup_ind[ip0+1] - psup_ind[ip0+0];
    const CVec3d n0 = CVec3d(aNorm.data()+ip0*3);
    CVec3d o_new = CVec3d(aOdir.data()+ip0*3);
    double weight = 0.0;
    {
      permutation1.resize(npj);
      for(unsigned int jjp=0;jjp<npj;++jjp){ permutation1[jjp] = jjp; }
      std::shuffle( permutation1.begin(), permutation1.end(), g);
    }
    for(unsigned int jjp=0;jjp<npj;++jjp){
      unsigned int jp1 = psup[psup_ind[ip0]+permutation1[jjp]];
      const CVec3d n1 = CVec3d(aNorm.data()+jp1*3);
      const CVec3d o1 = CVec3d(aOdir.data()+jp1*3);
      unsigned int idiff;
      FindNearestOrientation(idiff,
          o_new,n0, o1,n1,weight);
      weight += 1.0;
    }
    o_new.CopyTo(aOdir.data()+ip0*3);
  }
}

}


#endif //EXAMPLES_GLFWOLD_HDRONLY_4ROTSYM_H
