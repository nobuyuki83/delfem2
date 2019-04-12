#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#if defined(__APPLE__) && defined(__MACH__) // Mac
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/glu.h>
#elif defined(_WIN32) // windows
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#else // linux
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "delfem2/funcs_gl.h"
#include "delfem2/dyntri_v3.h"

#include "delfem2/mshio.h"
#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"

#include "delfem2/mshtopoio_gl.h"


void MeshTri3D_GeodesicPolyhedron
(std::vector<double>& aXYZ1,
 std::vector<int>& aTri1)
{
  std::vector<double> aXYZ0;
  std::vector<int> aTri0;
  MeshTri3D_Icosahedron(aXYZ0, aTri0);
  ////
  const int np0 = aXYZ0.size()/3;
  std::vector<int> elsup_ind, elsup;
  makeElemSurroundingPoint(elsup_ind, elsup,
                           aTri0.data(), aTri0.size()/3, 3, np0);
  ////
  std::vector<int> psup_ind, psup;
  makeOneRingNeighborhood(psup_ind, psup,
                          aTri0.data(),
                          elsup_ind, elsup,
                          3, np0);
  //  std::cout << "psup" << std::endl;
  //  Print_IndexedArray(psup_ind, psup);
  /////
  std::vector<int> edge_ind, edge;
  makeEdge(edge_ind, edge,
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

void CMeshElem::DrawFace_ElemWiseNorm() const
{
  if( elem_type == MESHELEM_TRI ){
    if( ndim == 3 ){ DrawMeshTri3D_FaceNorm(aPos, aElem); }
    if( ndim == 2 ){ DrawMeshTri2D_Face(aElem,aPos); }
  }
  else if( elem_type == MESHELEM_QUAD ){
    if( ndim == 3 ){ DrawMeshQuad3D_FaceNorm(aPos, aElem); }
    if( ndim == 2 ){ DrawMeshQuad2D_Edge(aPos, aElem); }
  }
}

void CMeshElem::DrawEdge() const {
  if( elem_type == MESHELEM_TRI ){
    if( ndim == 3 ){ DrawMeshTri3D_Edge(aPos, aElem); }
    if( ndim == 2 ){ DrawMeshTri2D_Edge(aElem, aPos); }
  }
  else if( elem_type == MESHELEM_QUAD ){
    if( ndim == 3 ){ DrawMeshQuad3D_Edge(aPos, aElem); }
  }
}

void CMeshElem::Draw() const {
  if(      color_face.size() == 4 ){
    glColor4d(color_face[0], color_face[1], color_face[2], color_face[4]);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color_face.data());
  }
  else if( color_face.size() == 3 ){
    const float color[4] = {color_face[0], color_face[1], color_face[2], 1.0};
    glColor4fv(color);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
  }
  /////
  this->DrawFace_ElemWiseNorm();
  if( this->is_draw_edge ) {
    glDisable(GL_LIGHTING);
    glLineWidth(1);
    this->DrawEdge();
  }
}

CMeshElem Read_MeshTri3D_Nas_CMeshElem(const std::string& fpath){
  CMeshElem em;
  em.elem_type = MESHELEM_TRI;
  em.ndim = 3;
  Read_MeshTri3D_Nas(em.aPos, em.aElem, fpath.c_str());
  return em;
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

void CMeshMultiElem::Draw() const
{
  for(int iogt=0;iogt<(int)aObjGroupTri.size();++iogt){
    int imtl = aObjGroupTri[iogt].imtl;
    ::glEnable(GL_LIGHTING);
    if( imtl>=0 && imtl<(int)aMaterial.size() ){ aMaterial[imtl].GL(); }
    else{ ::myGlColorDiffuse(CColor::White()); }
    DrawMeshTri3D_FaceNorm(aXYZ, aObjGroupTri[iogt].aTriVtx,
                          aNorm, aObjGroupTri[iogt].aTriNrm);
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

