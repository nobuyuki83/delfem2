/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include "delfem2/dyntri_v3.h"
#include "delfem2/mshio.h"
#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"
#include "delfem2/mshtopoio.h"

///////////////////////////////////

// probably std::stroi is safer to use but it is only for C++11
static int myStoi(const std::string& str){
  char* e;
  long d = std::strtol(str.c_str(),&e,0);
  return (int)d;
}

static double myStof(const std::string& str){
  char* e;
  float fval = std::strtof(str.c_str(),&e);
  return fval;
}

////////////////////////////////////


void MeshTri3D_GeodesicPolyhedron
(std::vector<double>& aXYZ1,
 std::vector<unsigned int>& aTri1)
{
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aTri0;
  MeshTri3D_Icosahedron(aXYZ0, aTri0);
  ////
  const int np0 = aXYZ0.size()/3;
  std::vector<int> elsup_ind, elsup;
  JArrayElemSurPoint_MeshElem(elsup_ind, elsup,
                              aTri0.data(), aTri0.size()/3, 3, np0);
  ////
  std::vector<int> psup_ind, psup;
  JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                                              aTri0.data(),
                                              elsup_ind, elsup,
                                              3, np0);
  //  std::cout << "psup" << std::endl;
  //  Print_IndexedArray(psup_ind, psup);
  /////
  std::vector<int> edge_ind, edge;
  JArrayEdgeUnidir_PointSurPoint(edge_ind, edge,
                             psup_ind,psup);
  //  std::cout << "edge" << std::endl;
  //  Print_IndexedArray(edge_ind, edge);
  ////
  double r0 = sqrt((5+sqrt(5))*0.5);
  aXYZ1 = aXYZ0;
  for(int ip=0;ip<np0;++ip){
    for(int iedge=edge_ind[ip];iedge<edge_ind[ip+1];++iedge){
      const int ip0 = edge[iedge];
      const double x1 = (aXYZ1[ip*3+0] + aXYZ1[ip0*3+0])*0.5;
      const double y1 = (aXYZ1[ip*3+1] + aXYZ1[ip0*3+1])*0.5;
      const double z1 = (aXYZ1[ip*3+2] + aXYZ1[ip0*3+2])*0.5;
      double mag = r0/sqrt(x1*x1+y1*y1+z1*z1);
      aXYZ1.push_back(x1*mag);
      aXYZ1.push_back(y1*mag);
      aXYZ1.push_back(z1*mag);
    }
  }
  aTri1.clear();
  aTri1.reserve(aTri0.size()*3);
  for(unsigned int itri=0;itri<aTri0.size()/3;++itri){
    const int ip0 = aTri0[itri*3+0];
    const int ip1 = aTri0[itri*3+1];
    const int ip2 = aTri0[itri*3+2];
    int iedge01,iedge12,iedge20;
    {
      if( ip0 < ip1 ){ iedge01 = findEdge(ip0,ip1, edge_ind,edge); }
      else {           iedge01 = findEdge(ip1,ip0, edge_ind,edge); }
      if( ip1 < ip2 ){ iedge12 = findEdge(ip1,ip2, edge_ind,edge); }
      else {           iedge12 = findEdge(ip2,ip1, edge_ind,edge); }
      if( ip2 < ip0 ){ iedge20 = findEdge(ip2,ip0, edge_ind,edge); }
      else {           iedge20 = findEdge(ip0,ip2, edge_ind,edge); }
    }
    aTri1.push_back(ip0); aTri1.push_back(iedge01+np0); aTri1.push_back(iedge20+np0);
    aTri1.push_back(ip1); aTri1.push_back(iedge12+np0); aTri1.push_back(iedge01+np0);
    aTri1.push_back(ip2); aTri1.push_back(iedge20+np0); aTri1.push_back(iedge12+np0);
    aTri1.push_back(iedge01+np0); aTri1.push_back(iedge12+np0); aTri1.push_back(iedge20+np0);
  }
}

////////

void CMeshMultiElem::ReadObj(const std::string& path_obj)
{
  std::string fname_mtl;
  Load_Obj(path_obj,
           fname_mtl, aXYZ, aNorm, aObjGroupTri);
  std::string path_dir = std::string(path_obj.begin(),path_obj.begin()+path_obj.rfind("/"));
  Load_Mtl(path_dir+"/"+fname_mtl,
           aMaterial);
//  std::cout << aObjGroupTri.size() << " " << aMaterial.size() << std::endl;
  { //
    std::map< std::string, int > mapMtlName2Ind;
    for(int imtl=0;imtl<(int)aMaterial.size();++imtl){
      mapMtlName2Ind.insert( std::make_pair(aMaterial[imtl].name_mtl, imtl) );
    }
    for(int iogt=0;iogt<(int)aObjGroupTri.size();++iogt){
      std::string name_mtl = aObjGroupTri[iogt].name_mtl;
      auto itr = mapMtlName2Ind.find(name_mtl);
      if( name_mtl.empty() || itr == mapMtlName2Ind.end() ){
        aObjGroupTri[iogt].imtl = -1;
        continue;
      }
      aObjGroupTri[iogt].imtl = itr->second;
    }
  }
}

std::vector<double> CMeshMultiElem::AABB3_MinMax() const
{
  double cw[6]; GetCenterWidth(cw, aXYZ);
  std::vector<double> aabb(6);
  aabb[0] = cw[0]-0.5*cw[3];
  aabb[1] = cw[0]+0.5*cw[3];
  aabb[2] = cw[1]-0.5*cw[4];
  aabb[3] = cw[1]+0.5*cw[4];
  aabb[4] = cw[2]-0.5*cw[5];
  aabb[5] = cw[2]+0.5*cw[5];
  return aabb;
}

void CMeshMultiElem::ScaleXYZ(double s)
{
  Scale(s,aXYZ);
}

void CMeshMultiElem::TranslateXYZ(double x, double y, double z)
{
  Translate(x,y,z, aXYZ);
}

void Load_Mtl
(const std::string& fname,
 std::vector<CMaterial>& aMtl)
{
  std::ifstream fin;
  fin.open(fname.c_str());
  if (fin.fail()){
    std::cout<<"File Read Fail"<<std::endl;
    return;
  }
  aMtl.clear();
  const int BUFF_SIZE = 256;
  char buff[BUFF_SIZE];
  while (fin.getline(buff, BUFF_SIZE)){
    if (buff[0]=='#'){ continue; }
    if (buff[0]=='\n'){ continue; }
    std::stringstream ss(buff);
    std::string str0, str1, str2, str3, str4;
    ss >> str0;
    if( str0 == "newmtl" ){
      aMtl.resize(aMtl.size()+1);
      const int imtl0 = aMtl.size()-1;
      ss >> str1;
      aMtl[imtl0].name_mtl = str1;
    }
    if( str0 == "Kd" ){
      const int imtl0 = aMtl.size()-1;
      ss >> str1 >> str2 >> str3;
      aMtl[imtl0].Kd[0] = myStof(str1);
      aMtl[imtl0].Kd[1] = myStof(str2);
      aMtl[imtl0].Kd[2] = myStof(str3);
      aMtl[imtl0].Kd[3] = 1.0;
    }
    if( str0 == "Ka" ){
      const int imtl0 = aMtl.size()-1;
      ss >> str1 >> str2 >> str3;
      aMtl[imtl0].Ka[0] = myStof(str1);
      aMtl[imtl0].Ka[1] = myStof(str2);
      aMtl[imtl0].Ka[2] = myStof(str3);
      aMtl[imtl0].Ka[3] = 1.0;
    }
    if( str0 == "Ks" ){
      const int imtl0 = aMtl.size()-1;
      ss >> str1 >> str2 >> str3;
      aMtl[imtl0].Ks[0] = myStof(str1);
      aMtl[imtl0].Ks[1] = myStof(str2);
      aMtl[imtl0].Ks[2] = myStof(str3);
      aMtl[imtl0].Ks[3] = 1.0;
    }
    if( str0 == "Ke" ){
      const int imtl0 = aMtl.size()-1;
      ss >> str1 >> str2 >> str3;
      aMtl[imtl0].Ke[0] = myStof(str1);
      aMtl[imtl0].Ke[1] = myStof(str2);
      aMtl[imtl0].Ke[2] = myStof(str3);
      aMtl[imtl0].Ke[3] = myStof(str3);
    }
    if( str0 == "Ns" ){
      const int imtl0 = aMtl.size()-1;
      ss >> str1;
      aMtl[imtl0].Ns = myStof(str1);
    }
    if( str0 == "illum" ){
      const int imtl0 = aMtl.size()-1;
      ss >> str1;
      aMtl[imtl0].illum = myStoi(str1);
    }
    if( str0 == "map_Kd" ){
      const int imtl0 = aMtl.size()-1;
      ss >> str1;
      aMtl[imtl0].map_Kd = str1;
    }
  }
}



