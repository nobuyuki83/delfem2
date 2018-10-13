#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#if defined(__APPLE__) && defined(__MACH__) // Mac
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/glu.h>
#elif defined(WIN32) // windows
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

void CMeshElem::DrawFace_ElemWiseNorm() const
{
  if( elem_type == MESHELEM_TRI ){
    if( ndim == 3 ){ DrawMeshTri3D_FaceNorm(aPos, aElem); }
    if( ndim == 2 ){ DrawMeshTri2D_Face(aElem,aPos); }
  }
  else if( elem_type == MESHELEM_QUAD ){
    if( ndim == 3 ){ DrawMeshQuad3D_FaceNorm(aPos, aElem); }
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
  glDisable(GL_LIGHTING);
  glLineWidth(1);
  this->DrawEdge();
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

/////////

CTriangulationOutput Triangulation
(const std::vector<double>& aXY,
 double edge_length)
{
  CTriangulationOutput out;
  std::vector< std::vector<double> > aaXY;
  aaXY.push_back(aXY);
  GenerateTesselation2(out.me.aElem, out.me.aPos,
                       out.aPtrVtxInd, out.aVtxInd,
                       edge_length, true, aaXY);
  out.me.ndim = 2;
  out.me.elem_type = MESHELEM_TRI;
  return out;
}

