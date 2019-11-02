#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <fstream>
#include <time.h>
#include "delfem2/adf.h"
#include "delfem2/msh.h"
#include "delfem2/mshio.h"
#include "delfem2/primitive.h"
#include "delfem2/bv.h"
#include "delfem2/bvh.h"
#include "delfem2/srch_v3bvhmshtopo.h"

// ----------------

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/gl2_color.h"
#include "delfem2/opengl/gl2_funcs.h"

// -----------------------

CADF3 adf;
std::vector<unsigned int> aTri;
std::vector<double> aXYZ;

// ---------------

void SetProblem(int iprob)
{
  if( iprob == 0 )
  {
    class CInSphere : public CInput_ADF3
    {
    public:
      virtual double sdf(double x, double y, double z) const {
        double n[3];
        return obj.Projection(n,
                              x, y, z);
      }
    public:
      CSphere obj;
    };
    CInSphere sphere;
    sphere.obj.radius_ = 0.5;
//    sphere.sdf.GetMesh(aTri, aXYZ, 0.01);
    double bb[6] = { -1, 1, -1, 1, -1,1 };
    adf.SetUp(sphere, bb);
    adf.BuildIsoSurface_MarchingCube();
    adf.BuildMarchingCubeEdge();
  }
  else if( iprob == 1 ){
    class CInTorus : public CInput_ADF3
    {
    public:
      virtual double sdf(double x, double y, double z) const {
        double n[3];
        return obj.Projection(n,
                              x, y, z);
      }
    public:
      CTorus obj;
    };
    CInTorus torus;
    torus.obj.radius_ = 0.5;
    torus.obj.radius_tube_ = 0.2;
//    torus.sdf.GetMesh(aTri, aXYZ, 0.01);
    double bb[6] = { -1, 1, -1, 1, -1,1 };
    adf.SetUp(torus, bb);
    adf.BuildIsoSurface_MarchingCube();
    adf.BuildMarchingCubeEdge();
  }
  else if( iprob == 2 ){
    class CMesh : public CInput_ADF3
    {
    public:
      virtual double sdf(double x, double y, double z) const {
        CVector3 n0;
        double sdf0 = obj.SignedDistanceFunction(n0,
                                                 CVector3(x,y,z),
                                                 aXYZ, aTri, aNorm);
        return sdf0;
      }
    public:
      std::vector<double> aNorm;
      CBVH_MeshTri3D<CBV3D_Sphere> obj;
    };
    CMesh mesh;
    {
      std::cout << PATH_INPUT_DIR << std::endl;
      Read_Ply(std::string(PATH_INPUT_DIR)+"/bunny_1k.ply", aXYZ, aTri);
      Normalize(aXYZ,1.7);
      mesh.obj.Init(aXYZ.data(), aXYZ.size()/3,
                    aTri.data(), aTri.size()/3,
                    0.0);
      mesh.aNorm.resize(aXYZ.size());
      Normal_MeshTri3D(mesh.aNorm.data(),
                       aXYZ.data(), aXYZ.size()/3, aTri.data(), aTri.size()/3);
    }
    double bb[6] = { -1, 1, -1, 1, -1,1 };
    adf.SetUp(mesh, bb);
    adf.BuildIsoSurface_MarchingCube();
    adf.BuildMarchingCubeEdge();
  }
}

// ------------------------------------------------

int main(int argc,char* argv[])
{
  CViewer_GLFW viewer;
  viewer.Init_GLold();
  
  viewer.nav.camera.view_height = 2.0;
  viewer.nav.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  opengl::setSomeLighting();
  
  while(!glfwWindowShouldClose(viewer.window)){
    int iproblem  = 0;
    {
      static int iframe = 0;
      iproblem = iframe/1000;
      if( iframe % 1000  == 0 ){ SetProblem(iproblem); }
      iframe = (iframe+1)%3000;
    }
    // --------------------
    viewer.DrawBegin_Glold();
    if( iproblem == 0 ){
      adf.SetShowCage(false);
      adf.Draw();
    }
    else if( iproblem == 1 ){
      adf.SetShowCage(true);
      adf.Draw();
    }
    else if( iproblem == 2 ){
//      opengl::DrawMeshTri3D_FaceNorm(aXYZ,aTri);
      adf.Draw();
    }
    viewer.DrawEnd_oldGL();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
