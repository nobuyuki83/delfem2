/*
 * Copyright (c) 2019-2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/color_glold.h"
#include "delfem2/points.h"
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/color.h"
#include <GLFW/glfw3.h>
#include <vector>
#include <string>
#include <cassert>
#include <cstdlib>
#include <random>
#include <set>
#include <queue>

namespace dfm2 = delfem2;

class CNode_MeshGeodesic {
public:
  CNode_MeshGeodesic(unsigned int itri_, unsigned int idist_)
      : itri(itri_), idist(idist_) {}
  bool operator < (const CNode_MeshGeodesic& lhs) const {
    return this->idist > lhs.idist; // smaller distance have more priority
  }
public:
  unsigned int itri;
  unsigned int idist;
};

void MeshGeodesic(
    std::vector<unsigned int>& aDist,
    //
    unsigned int itri_ker,
    const std::vector<unsigned int>& aTriSuTri)
{
  const unsigned int nelem = aDist.size();
  aDist.assign(nelem,UINT_MAX);
  aDist[itri_ker] = 0;
  const unsigned int nedge = aTriSuTri.size()/nelem;
  std::priority_queue<CNode_MeshGeodesic> que;
  que.push(CNode_MeshGeodesic(itri_ker,0));
  while(!que.empty()) {
    unsigned int itri_fix = que.top().itri;
    unsigned int idist_fix = que.top().idist;
    que.pop();
    assert( aDist[itri_fix] == idist_fix );
    for(unsigned int iedtri=0;iedtri<nedge;++iedtri){
      unsigned int itri1 = aTriSuTri[itri_fix*nedge+iedtri];
      if( itri1 == UINT_MAX){ continue; }
      if( idist_fix+1 >= aDist[itri1] ){ continue; }
      // Found the shortest path ever examined.
      aDist[itri1] = idist_fix+1;
      // put this in the que because this is a candidate for the shortest path
      que.push(CNode_MeshGeodesic(itri1,idist_fix+1));
    }
  }
}

void MeshClustering(
    std::vector<unsigned int>& aFlgElm,
    //
    unsigned int ncluster,
    const std::vector<unsigned int>& aTriSuTri,
    unsigned int ntri)
{
  std::vector<unsigned int> aDist0(ntri,UINT_MAX);
  std::random_device rd;
  std::mt19937 rdeng(rd());
  std::uniform_int_distribution<unsigned int> dist0(0,ntri-1);
  const unsigned int itri_ker = dist0(rdeng); assert(itri_ker<ntri);
  aFlgElm.assign(ntri,0);
  MeshGeodesic(aDist0,
               itri_ker,aTriSuTri);
  for(unsigned int icluster=1;icluster<ncluster;++icluster){
    unsigned int itri_maxdist;
    { // find triangle with maximum distance
      double idist_max = 0;
      for(unsigned int it=0;it<ntri;++it){
        if( aDist0[it] <= idist_max ){ continue; }
        idist_max = aDist0[it];
        itri_maxdist = it;
      }
    }
    std::vector<unsigned int> aDist1(ntri,UINT_MAX);
    MeshGeodesic(aDist1,
                 itri_maxdist,aTriSuTri);
    for(unsigned int it=0;it<ntri;++it) {
      if( aDist1[it] < aDist0[it] ){
        aDist0[it] = aDist1[it];
        aFlgElm[it] = icluster;
      }
    }
  }
}

// ---------------------------

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;

  delfem2::Read_Ply(
//      std::string(PATH_INPUT_DIR)+"/bunny_34k.ply",
      std::string(PATH_INPUT_DIR)+"/arm_16k.ply",
      aXYZ,aTri);
  delfem2::Normalize_Points3(aXYZ);
  std::vector<unsigned int> aTriSuTri;
  ElSuEl_MeshElem(aTriSuTri,
      aTri.data(), aTri.size()/3,
      delfem2::MESHELEM_TRI,
      aXYZ.size()/3);

  // ------

  std::random_device rd;
  std::mt19937 rdeng(rd());
  std::uniform_int_distribution<unsigned int> ncluster_gen(1,100);

  // above: data preparation
  // -----------------------
  // below: view

  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  delfem2::opengl::setSomeLighting();
  viewer.nav.camera.view_height = 0.5;
  viewer.nav.camera.camera_rot_mode  = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  
  while (!glfwWindowShouldClose(viewer.window))
  {
    unsigned int ncluster = ncluster_gen(rdeng);
    std::vector< std::pair<int,dfm2::CColor> > aColor;
    for(unsigned int ic=0;ic<ncluster;++ic){
      dfm2::CColor c;
      c.setRandomVividColor();
      aColor.emplace_back(2,c);
    }
    std::vector<unsigned int> aFlgElm;
    MeshClustering(aFlgElm,ncluster,aTriSuTri,aTri.size()/3);
    //
    for(unsigned int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      ::glDisable(GL_LIGHTING);
      ::glColor3d(0, 0, 0);
      delfem2::opengl::DrawMeshTri3D_Edge(aXYZ, aTri);
      delfem2::opengl::DrawMeshTri3DFlag_FaceNorm(aXYZ, aTri, aFlgElm, aColor);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
    }
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
