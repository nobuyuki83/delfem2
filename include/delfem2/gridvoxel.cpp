/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/vec3.h"
#include "delfem2/gridvoxel.h"

namespace dfm2 = delfem2;

const unsigned int noelElemFace_Vox[6][4] = {
  { 0, 4, 6, 2 }, // -x
  { 1, 3, 7, 5 }, // +x
  { 0, 1, 5, 4 }, // -y
  { 2, 6, 7, 3 }, // +y
  { 0, 2, 3, 1 }, // -z
  { 4, 5, 7, 6 } }; // +z

const dfm2::CVec3d normalHexFace[6] = {
  dfm2::CVec3d(-1, 0, 0),
  dfm2::CVec3d(+1, 0, 0),
  dfm2::CVec3d( 0,-1, 0),
  dfm2::CVec3d( 0,+1, 0),
  dfm2::CVec3d( 0, 0,-1),
  dfm2::CVec3d( 0, 0,+1)
};

bool dfm2::IsInclude_AABB(const int aabb[8], int igvx, int igvy, int igvz)
{
  if( igvx < aabb[0] || igvx >= aabb[1] ){ return false; }
  if( igvy < aabb[2] || igvy >= aabb[3] ){ return false; }
  if( igvz < aabb[4] || igvz >= aabb[5] ){ return false; }
  return true;
}

void dfm2::Add_AABB(int aabb[8], int ivx, int ivy, int ivz)
{
  const int ipx0 = ivx+0;  const int ipx1 = ivx+1;
  const int ipy0 = ivy+0;  const int ipy1 = ivy+1;
  const int ipz0 = ivz+0;  const int ipz1 = ivz+1;
  if( aabb[0]>aabb[1] ){
    aabb[0] = ipx0;  aabb[1] = ipx1;
    aabb[2] = ipy0;  aabb[3] = ipy1;
    aabb[4] = ipz0;  aabb[5] = ipz1;
  }
  else{
    if( ipx0 < aabb[0] ){ aabb[0] = ipx0; }
    if( ipx1 > aabb[1] ){ aabb[1] = ipx1; }
    if( ipy0 < aabb[2] ){ aabb[2] = ipy0; }
    if( ipy1 > aabb[3] ){ aabb[3] = ipy1; }
    if( ipz0 < aabb[4] ){ aabb[4] = ipz0; }
    if( ipz1 > aabb[5] ){ aabb[5] = ipz1; }
  }
}


void dfm2::MeshQuad3D_VoxelGrid(
    std::vector<double>& aXYZ,
    std::vector<unsigned int>& aQuad,
    unsigned int ndivx,
    unsigned int ndivy,
    unsigned int ndivz,
    const std::vector<int>& aIsVox)
{
  int ioffx=0;
  int ioffy=0;
  int ioffz=0;
  aQuad.clear();
  aXYZ.clear();
  //
  const unsigned int mdivx = ndivx+1;
  const unsigned int mdivy = ndivy+1;
  const unsigned int mdivz = ndivz+1;
  for(unsigned int igpx=0;igpx<mdivx;++igpx){
    for(unsigned int igpy=0;igpy<mdivy;++igpy){
      for(unsigned int igpz=0;igpz<mdivz;++igpz){
        aXYZ.push_back( igpx+ioffx );
        aXYZ.push_back( igpy+ioffy );
        aXYZ.push_back( igpz+ioffz );
      }
    }
  }
  //
  assert( aIsVox.size() == ndivx*ndivy*ndivz );
  for(unsigned int igvx=0;igvx<ndivx;++igvx){
    for(unsigned int igvy=0;igvy<ndivy;++igvy){
      for(unsigned int igvz=0;igvz<ndivz;++igvz){
        const unsigned int ivoxel = igvx*(ndivy*ndivz)+igvy*ndivz+igvz;
        assert( ivoxel < aIsVox.size() );
        if( aIsVox[ivoxel] == 0 ){ continue; }
        //
        int aIGP_Vox[8] = {0,0,0,0, 0,0,0,0};
        {
          aIGP_Vox[0] = (igvx+0)*(mdivy*mdivz)+(igvy+0)*mdivz+(igvz+0);
          aIGP_Vox[1] = (igvx+1)*(mdivy*mdivz)+(igvy+0)*mdivz+(igvz+0);
          aIGP_Vox[2] = (igvx+0)*(mdivy*mdivz)+(igvy+1)*mdivz+(igvz+0);
          aIGP_Vox[3] = (igvx+1)*(mdivy*mdivz)+(igvy+1)*mdivz+(igvz+0);
          aIGP_Vox[4] = (igvx+0)*(mdivy*mdivz)+(igvy+0)*mdivz+(igvz+1);
          aIGP_Vox[5] = (igvx+1)*(mdivy*mdivz)+(igvy+0)*mdivz+(igvz+1);
          aIGP_Vox[6] = (igvx+0)*(mdivy*mdivz)+(igvy+1)*mdivz+(igvz+1);
          aIGP_Vox[7] = (igvx+1)*(mdivy*mdivz)+(igvy+1)*mdivz+(igvz+1);
        }
        for(int iface=0;iface<6;++iface){
          const int jgv0 = Adj_Grid(ivoxel, iface, ndivx, ndivy, ndivz);
          if( jgv0 >= 0 ){
            assert( jgv0 < (int)aIsVox.size() );
            if( aIsVox[jgv0] == 1 ){ continue; } // facing to adjacent voxel -> no outward face.
          }
          //
          const int aIGP0 = aIGP_Vox[ noelElemFace_Vox[iface][0] ];
          const int aIGP1 = aIGP_Vox[ noelElemFace_Vox[iface][1] ];
          const int aIGP2 = aIGP_Vox[ noelElemFace_Vox[iface][2] ];
          const int aIGP3 = aIGP_Vox[ noelElemFace_Vox[iface][3] ];
          aQuad.push_back(aIGP0);
          aQuad.push_back(aIGP1);
          aQuad.push_back(aIGP2);
          aQuad.push_back(aIGP3);
        }
      }
    }
  }
}

void dfm2::MeshHex3D_VoxelGrid(
    std::vector<double>& aXYZ,
    std::vector<int>& aHex,
    unsigned int ndivx,
    unsigned int ndivy,
    unsigned int ndivz,
    const std::vector<int>& aIsVox)
{
  int ioffx=0;
  int ioffy=0;
  int ioffz=0;
  aHex.clear();
  aXYZ.clear();
  //
  const unsigned int mdivx = ndivx+1;
  const unsigned int mdivy = ndivy+1;
  const unsigned int mdivz = ndivz+1;
  for(unsigned int igpx=0;igpx<mdivx;++igpx){
    for(unsigned int igpy=0;igpy<mdivy;++igpy){
      for(unsigned int igpz=0;igpz<mdivz;++igpz){
        aXYZ.push_back( igpx+ioffx );
        aXYZ.push_back( igpy+ioffy );
        aXYZ.push_back( igpz+ioffz );
      }
    }
  }
  //////
  assert( aIsVox.size() == ndivx*ndivy*ndivz );
  for(unsigned int igvx=0;igvx<ndivx;++igvx){
    for(unsigned int igvy=0;igvy<ndivy;++igvy){
      for(unsigned int igvz=0;igvz<ndivz;++igvz){
        const unsigned int ivoxel = igvx*(ndivy*ndivz)+igvy*ndivz+igvz;
        assert( ivoxel < aIsVox.size() );
        if( aIsVox[ivoxel] == 0 ){ continue; }
        //
        int aIGP_Hex[8] = {0,0,0,0, 0,0,0,0};
        {
          aIGP_Hex[0] = (igvx+0)*(mdivy*mdivz)+(igvy+0)*mdivz+(igvz+0);
          aIGP_Hex[1] = (igvx+1)*(mdivy*mdivz)+(igvy+0)*mdivz+(igvz+0);
          aIGP_Hex[2] = (igvx+1)*(mdivy*mdivz)+(igvy+1)*mdivz+(igvz+0);
          aIGP_Hex[3] = (igvx+0)*(mdivy*mdivz)+(igvy+1)*mdivz+(igvz+0);
          aIGP_Hex[4] = (igvx+0)*(mdivy*mdivz)+(igvy+0)*mdivz+(igvz+1);
          aIGP_Hex[5] = (igvx+1)*(mdivy*mdivz)+(igvy+0)*mdivz+(igvz+1);
          aIGP_Hex[6] = (igvx+1)*(mdivy*mdivz)+(igvy+1)*mdivz+(igvz+1);
          aIGP_Hex[7] = (igvx+0)*(mdivy*mdivz)+(igvy+1)*mdivz+(igvz+1);
        }
        aHex.push_back(aIGP_Hex[0]);
        aHex.push_back(aIGP_Hex[1]);
        aHex.push_back(aIGP_Hex[2]);
        aHex.push_back(aIGP_Hex[3]);
        aHex.push_back(aIGP_Hex[4]);
        aHex.push_back(aIGP_Hex[5]);
        aHex.push_back(aIGP_Hex[6]);
        aHex.push_back(aIGP_Hex[7]);
      }
    }
  }
}



void dfm2::MeshTet3D_VoxelGrid(
    std::vector<double>& aXYZ,
    std::vector<int>& aTet,
    unsigned int ndivx,
    unsigned int ndivy,
    unsigned int ndivz,
    const std::vector<int>& aIsVox)
{
  int ioffx=0;
  int ioffy=0;
  int ioffz=0;
  aTet.clear();
  aXYZ.clear();
  //
  const unsigned int mdivx = ndivx+1;
  const unsigned int mdivy = ndivy+1;
  const unsigned int mdivz = ndivz+1;
  for(unsigned int igpx=0;igpx<mdivx;++igpx){
    for(unsigned int igpy=0;igpy<mdivy;++igpy){
      for(unsigned int igpz=0;igpz<mdivz;++igpz){
        aXYZ.push_back( igpx+ioffx );
        aXYZ.push_back( igpy+ioffy );
        aXYZ.push_back( igpz+ioffz );
      }
    }
  }
  //
  assert( aIsVox.size() == ndivx*ndivy*ndivz );
  for(unsigned int igvx=0;igvx<ndivx;++igvx){
    for(unsigned int igvy=0;igvy<ndivy;++igvy){
      for(unsigned int igvz=0;igvz<ndivz;++igvz){
        const unsigned int ivoxel = igvx*(ndivy*ndivz)+igvy*ndivz+igvz;
        assert( ivoxel < aIsVox.size() );
        if( aIsVox[ivoxel] == 0 ){ continue; }
        //
        int aIP[8] = {0,0,0,0, 0,0,0,0};
        {
          aIP[0] = (igvx+0)*(mdivy*mdivz)+(igvy+0)*mdivz+(igvz+0);
          aIP[1] = (igvx+1)*(mdivy*mdivz)+(igvy+0)*mdivz+(igvz+0);
          aIP[2] = (igvx+1)*(mdivy*mdivz)+(igvy+1)*mdivz+(igvz+0);
          aIP[3] = (igvx+0)*(mdivy*mdivz)+(igvy+1)*mdivz+(igvz+0);
          aIP[4] = (igvx+0)*(mdivy*mdivz)+(igvy+0)*mdivz+(igvz+1);
          aIP[5] = (igvx+1)*(mdivy*mdivz)+(igvy+0)*mdivz+(igvz+1);
          aIP[6] = (igvx+1)*(mdivy*mdivz)+(igvy+1)*mdivz+(igvz+1);
          aIP[7] = (igvx+0)*(mdivy*mdivz)+(igvy+1)*mdivz+(igvz+1);
        }
        aTet.push_back(aIP[0]);
        aTet.push_back(aIP[1]);
        aTet.push_back(aIP[2]);
        aTet.push_back(aIP[6]);
        ////
        aTet.push_back(aIP[0]);
        aTet.push_back(aIP[1]);
        aTet.push_back(aIP[6]);
        aTet.push_back(aIP[5]);
        ////
        aTet.push_back(aIP[0]);
        aTet.push_back(aIP[4]);
        aTet.push_back(aIP[5]);
        aTet.push_back(aIP[6]);
        ////
        aTet.push_back(aIP[0]);
        aTet.push_back(aIP[2]);
        aTet.push_back(aIP[3]);
        aTet.push_back(aIP[6]);
        ////
        aTet.push_back(aIP[0]);
        aTet.push_back(aIP[3]);
        aTet.push_back(aIP[7]);
        aTet.push_back(aIP[6]);
        ////
        aTet.push_back(aIP[0]);
        aTet.push_back(aIP[4]);
        aTet.push_back(aIP[6]);
        aTet.push_back(aIP[7]);
      }
    }
  }
}

int dfm2::Adj_Grid
(int igridvox, int iface,
 int ndivx, int ndivy, int ndivz)
{
  int ivx0 = igridvox/(ndivy*ndivz);
  int ivy0 = (igridvox-ivx0*(ndivy*ndivz))/ndivz;
  int ivz0 = igridvox-ivx0*(ndivy*ndivz)-ivy0*ndivz;
  if( iface == 0 ){ ivx0 -= 1; }
  if( iface == 1 ){ ivx0 += 1; }
  if( iface == 2 ){ ivy0 -= 1; }
  if( iface == 3 ){ ivy0 += 1; }
  if( iface == 4 ){ ivz0 -= 1; }
  if( iface == 5 ){ ivz0 += 1; }
  if( ivx0 < 0 || ivx0 >= ndivx ){ return -1; }
  if( ivy0 < 0 || ivy0 >= ndivy ){ return -1; }
  if( ivz0 < 0 || ivz0 >= ndivz ){ return -1; }
  int igv1 = ivx0*(ndivy*ndivz)+ivy0*ndivz+ivz0;
  assert( igv1 >= 0 && igv1 < ndivx*ndivy*ndivz );
  return igv1;
}

// ---------------------------------------------------------------------


void delfem2::Grid3Voxel_Dilation
 (dfm2::CGrid3<int>& grid)
{
  const int nx = (int)grid.ndivx;
  const int ny = (int)grid.ndivy;
  const int nz = (int)grid.ndivz;
  for(int iz1=0;iz1<nz;++iz1){
    for(int iy1=0;iy1<ny;++iy1){
      for(int ix1=0;ix1<nx;++ix1){
        int flgx0 = ix1-1>=0 ? grid.aVal[iz1*ny*nx + iy1*nx + (ix1-1)] : 0;
        int flgx2 = ix1+1<nx ? grid.aVal[iz1*ny*nx + iy1*nx + (ix1+1)] : 0;
        int flgy0 = iy1-1>=0 ? grid.aVal[iz1*ny*nx + (iy1-1)*nx + ix1] : 0;
        int flgy2 = iy1+1<ny ? grid.aVal[iz1*ny*nx + (iy1+1)*nx + ix1] : 0;
        int flgz0 = iz1-1>=0 ? grid.aVal[(iz1-1)*ny*nx + iy1*nx + ix1] : 0;
        int flgz2 = iz1+1<nz ? grid.aVal[(iz1+1)*ny*nx + iy1*nx + ix1] : 0;
        const int ivox = iz1*ny*nx + iy1*nx + ix1;
        if( grid.aVal[ivox] != 0 ) continue;
        if( flgx0 == 1 ){ grid.aVal[ivox] = 2; }
        if( flgx2 == 1 ){ grid.aVal[ivox] = 2; }
        if( flgy0 == 1 ){ grid.aVal[ivox] = 2; }
        if( flgy2 == 1 ){ grid.aVal[ivox] = 2; }
        if( flgz0 == 1 ){ grid.aVal[ivox] = 2; }
        if( flgz2 == 1 ){ grid.aVal[ivox] = 2; }
      }
    }
  }
  const int nvox = nx*ny*nz;
  for(int ivox=0;ivox<nvox;++ivox){
    if( grid.aVal[ivox] == 0 ){ continue; }
    grid.aVal[ivox] = 1;
  }
}

void delfem2::Grid3Voxel_Erosion
(dfm2::CGrid3<int>& grid)
{
  const int nx = (int)grid.ndivx;
  const int ny = (int)grid.ndivy;
  const int nz = (int)grid.ndivz;
  for(int iz1=0;iz1<nz;++iz1){
    for(int iy1=0;iy1<ny;++iy1){
      for(int ix1=0;ix1<nx;++ix1){
        int flgx0 = ix1-1>=0 ? grid.aVal[iz1*ny*nx + iy1*nx + (ix1-1)] : 0;
        int flgx2 = ix1+1<nx ? grid.aVal[iz1*ny*nx + iy1*nx + (ix1+1)] : 0;
        int flgy0 = iy1-1>=0 ? grid.aVal[iz1*ny*nx + (iy1-1)*nx + ix1] : 0;
        int flgy2 = iy1+1<ny ? grid.aVal[iz1*ny*nx + (iy1+1)*nx + ix1] : 0;
        int flgz0 = iz1-1>=0 ? grid.aVal[(iz1-1)*ny*nx + iy1*nx + ix1] : 0;
        int flgz2 = iz1+1<nz ? grid.aVal[(iz1+1)*ny*nx + iy1*nx + ix1] : 0;
        const int ivox = iz1*ny*nx + iy1*nx + ix1;
        if( grid.aVal[ivox] == 0 ) continue;
        if( flgx0 == 0 ){ grid.aVal[ivox] = 2; }
        if( flgx2 == 0 ){ grid.aVal[ivox] = 2; }
        if( flgy0 == 0 ){ grid.aVal[ivox] = 2; }
        if( flgy2 == 0 ){ grid.aVal[ivox] = 2; }
        if( flgz0 == 0 ){ grid.aVal[ivox] = 2; }
        if( flgz2 == 0 ){ grid.aVal[ivox] = 2; }
      }
    }
  }
  const int nvox = nx*ny*nz;
  for(int ivox=0;ivox<nvox;++ivox){
    if( grid.aVal[ivox] != 2 ){ continue; }
    grid.aVal[ivox] = 0;
  }
}
