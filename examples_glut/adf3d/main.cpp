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

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/opengl/gl2_color.h"
#include "delfem2/opengl/gl2_funcs.h"
#include "../glut_cam.h"

// -----------------------

CADF3 adf;
CNav3D_GLUT nav;
bool is_animation;
double cur_time = 0;

int imode_display = 0;
std::vector<unsigned int> aTri;
std::vector<double> aXYZ;

// ---------------

void SetProblem()
{
  const unsigned int nprob = 3;	// number of problems
  static int iprob = 0;
  
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
  iprob++;
  if( iprob == nprob ){ iprob = 0; }
}

///////////////////////////////

void myGlutResize(int w, int h)
{
  ::glViewport(0, 0, w, h);
}

void myGlutMotion( int x, int y )
{
  nav.glutMotion(x,y);
}

void myGlutMouse(int button, int state, int x, int y)
{
  nav.glutMouse(button,state,x,y);
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
  {
    case 'q':
    case 'Q':
    case '\033':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
    case 'a':
      is_animation = !is_animation;
      break;
    case 'd':
      imode_display = (imode_display+1)%3;
      break;
    case ' ':	
      SetProblem();
      break;
  }
}

void myGlutSpecial(int Key, int x, int y)
{
  nav.glutSpecial(Key,x,y);
}

void myGlutIdle(){
  ::glutPostRedisplay();
}

void myGlutDisplay(void)
{
  ::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  
  nav.SetGL_Camera();
  
  if( imode_display == 0 ){
    adf.SetShowCage(false);
    adf.Draw();
  }
  else if( imode_display == 1 ){
    adf.SetShowCage(true);
    adf.Draw();
  }
  else if( imode_display == 2 ){
    //    DrawTri3D_SurfaceNorm(aXYZ, aTri);
    DrawMeshTri3D_FaceNorm(aXYZ,aTri);
    //    Draw_Edge(aXYZ, aTri);
  }
  
  ShowFPS();
  ::glutSwapBuffers();
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int argc,char* argv[])
{	
  // Initialize GLUT
  glutInitWindowPosition(200,200);
  glutInitWindowSize(400, 300);
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  glutCreateWindow("FEM View");
  
  // Setting call back function
  glutDisplayFunc(myGlutDisplay);
  glutReshapeFunc(myGlutResize);
  glutMotionFunc(myGlutMotion);
  glutMouseFunc(myGlutMouse);
  glutKeyboardFunc(myGlutKeyboard);
  glutSpecialFunc(myGlutSpecial);
  glutIdleFunc(myGlutIdle);
  
  nav.camera.view_height = 2.0;
  nav.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  SetProblem();
  setSomeLighting();
  glutMainLoop();
  return 0;
}
