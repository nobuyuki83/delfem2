#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits>
#include <vector>
#include <set>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/msh.h" // Rotate, NormalMeshTri3D
#include "delfem2/mshtopo.h"
#include "delfem2/dyntri.h"
#include "delfem2/bv.h"
#include "delfem2/bvh.h"
#include "delfem2/primitive.h"

#include "delfem2/objfunc_v23.h"
#include "delfem2/objfunc_v23dtri.h"
#include "delfem2/v23m3q.h"
#include "delfem2/dyntri_v2.h"
#include "delfem2/cad2d.h"
#include "delfem2/srch_v3bvhmshtopo.h"

#include "delfem2/gl_cad_dyntri_v23.h"
#include "delfem2/gl_funcs.h"
#include "delfem2/glut_funcs.h"




////////////////////////////////////////////////////////////////////////////////////

std::vector<CEPo2> aPo2D;
std::vector<CVector2> aVec2;
std::vector<ETri> aETri;
std::vector<unsigned int> aLine;
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aXYZt;
std::vector<double> aUVW; // deformed vertex velocity
std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)
std::vector<CInfoNearest> aInfoNearest;


std::vector<double> aXYZ_Contact;
std::vector<unsigned int> aTri_Contact;
std::vector<double> aNorm_Contact(aXYZ.size());
CBVH_MeshTri3D<CBV3D_Sphere> bvh;

const double dt = 0.01;
const double gravity[3] = {0.0, 0.0, 0.0};
const double contact_clearance = 0.0001;

CGlutWindowManager win;
bool is_animation = false;

//////////////////////////////////

void StepTime()
{
  PBD_Pre3D(aXYZt,
            dt, gravity, aXYZ, aUVW, aBCFlag);
  PBD_TriStrain(aXYZt.data(),
                aXYZt.size()/3,aETri, aVec2);
  PBD_Bend(aXYZt.data(),
           aXYZt.size()/3, aETri, aVec2);
  PBD_Seam(aXYZt.data(),
           aXYZt.size()/3, aLine.data(), aLine.size()/2);
  Project_PointsIncludedInBVH_Outside_Cache(aXYZt.data(), aInfoNearest,
                                            aXYZt.size()/3,
                                            contact_clearance,bvh,
                                            aXYZ_Contact.data(), aXYZ_Contact.size()/3,
                                            aTri_Contact.data(), aTri_Contact.size()/3,
                                            aNorm_Contact.data(), 0.1);
  PBD_Post(aXYZ, aUVW,
           dt, aXYZt, aBCFlag);

}

//////////////////////////////////

void myGlutResize(int w, int h)
{
  glViewport(0, 0, w, h);
  ::glutPostRedisplay();
}



void myGlutDisplay(void)
{
  ::glClearColor(1.0, 1.0, 1.0, 1.0);
  //  ::glClearColor(0.0, .0, 0.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1, 4.0 );
  
  win.SetGL_Camera();
  
  ::glPointSize(5);
  ::glLineWidth(1);
  {
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0.8, 0.8, 0.8);
    /*
    float color[4] = {200.0/256.0, 200.0/256.0, 200.0/256.0,1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
    ::glEnable(GL_DEPTH_TEST);
     */
//    DrawMeshTri3D_FaceNorm(aXYZ, aTri);
  }
  
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  DrawMeshDynTri3D_Edge(aXYZ, aETri);

  ::glColor3d(1,0,0);
  DrawMeshTri3D_Edge(aXYZ_Contact.data(), aXYZ_Contact.size()/3,
                     aTri_Contact.data(), aTri_Contact.size()/3);
//  DrawSphere_Edge(rad0);

  
  ShowFPS();
  
  glutSwapBuffers();
}


void myGlutIdle(){
  if( is_animation ){
    StepTime();
  }
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y )
{
  win.glutMotion(x,y);
}

void myGlutMouse(int button, int state, int x, int y)
{
  win.glutMouse(button, state, x,y);
}

void myGlutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
    case 'q':
    case 'Q':
    case '\033':  /* '\033' ÇÕ ESC ÇÃ ASCII ÉRÅ[Éh */
      exit(0);
      break;
    case 'a':
    {
      is_animation = !is_animation;
      break;
    }
  }
}

void myGlutSpecial(int key, int x, int y){
  win.glutSpecial(key,x,y);
}


class CRigidTrans_2DTo3D
{
public:
  CVector2 org2;
  CVector3 org3;
  CMatrix3 R;
};

int main(int argc,char* argv[])
{
  // Initailze GLUT
  ::glutInitWindowPosition(200,200);
  ::glutInitWindowSize(400, 300);
  ::glutInit(&argc, argv);
  ::glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  ::glutCreateWindow("Cad View");
  
  // Set callback function
  ::glutMotionFunc(myGlutMotion);
  ::glutMouseFunc(myGlutMouse);
  ::glutDisplayFunc(myGlutDisplay);
  ::glutReshapeFunc(myGlutResize);
  ::glutKeyboardFunc(myGlutKeyboard);
  ::glutSpecialFunc(myGlutSpecial);
  ::glutIdleFunc(myGlutIdle);
  
  CCad2D cad;
  {
    const double xys0[8] = { -0.0,-0.0,  +1.0,-0.0, +1.0,+1.0, -0.0,+1.0, };
    const double xys1[8] = { +2.0,-0.0,  +3.0,-0.0, +3.0,+1.0, +2.0,+1.0, };
    cad.AddPolygon(std::vector<double>(xys0,xys0+8));
    cad.AddPolygon(std::vector<double>(xys1,xys1+8));
  }
  CMesher_Cad2D mesher;
  mesher.edge_length = 0.04;
  CMeshDynTri2D dmesh;
  mesher.Meshing(dmesh,
                 cad);
  dmesh.Check();
  aPo2D = dmesh.aEPo;
  aETri = dmesh.aETri;
  aVec2 = dmesh.aVec2;

  /////////
  const int np = aPo2D.size();
  aUVW.resize(np*3,0.0);
  aBCFlag.resize(np,0);
  aXYZ.resize(np*3);
  {
    CRigidTrans_2DTo3D rt23;
    rt23.org2 = CVector2(2.5,0.5);
    rt23.org3 = CVector3(0.0, 0.0, 0.5);
    rt23.R.SetRotMatrix_Cartesian(0.0, 3.1415, 0.0);
    std::vector<int> aIP = mesher.IndPoint_IndFaceArray(std::vector<int>(1,1), cad);
    for(int iip=0;iip<aIP.size();++iip){
      const int ip = aIP[iip];
      CVector3 p0(aVec2[ip].x-rt23.org2.x, aVec2[ip].y-rt23.org2.y,0.0);
      CVector3 p1 = rt23.org3+rt23.R*p0;
      aXYZ[ip*3+0] = p1.x;
      aXYZ[ip*3+1] = p1.y;
      aXYZ[ip*3+2] = p1.z;
    }
    {
      CRigidTrans_2DTo3D rt23;
      rt23.org2 = CVector2(0.5,0.5);
      rt23.org3 = CVector3(0.0, 0.0, -0.5);
      rt23.R.SetIdentity();
      std::vector<int> aIP = mesher.IndPoint_IndFaceArray(std::vector<int>(1,0), cad);
      for(int iip=0;iip<aIP.size();++iip){
        const int ip = aIP[iip];
        CVector3 p0(aVec2[ip].x-rt23.org2.x, aVec2[ip].y-rt23.org2.y,0.0);
        CVector3 p1 = rt23.org3+rt23.R*p0;
        aXYZ[ip*3+0] = p1.x;
        aXYZ[ip*3+1] = p1.y;
        aXYZ[ip*3+2] = p1.z;
      }
    }
    aLine.clear();
    {
      std::vector<unsigned int> aIP0 = mesher.IndPoint_IndEdge(1, true, cad);
      std::vector<unsigned int> aIP1 = mesher.IndPoint_IndEdge(7, true, cad);
      const int npe = aIP0.size();
      assert( aIP1.size() == npe );
      for(int iip=0;iip<npe;++iip){
        int ip0 = aIP0[iip];
        int ip1 = aIP1[npe-iip-1];
        aLine.push_back(ip0);
        aLine.push_back(ip1);
      }
    }
    {
      std::vector<unsigned int> aIP0 = mesher.IndPoint_IndEdge(3, true, cad);
      std::vector<unsigned int> aIP1 = mesher.IndPoint_IndEdge(5, true, cad);
      const int npe = aIP0.size();
      assert( aIP1.size() == npe );
      for(int iip=0;iip<npe;++iip){
        int ip0 = aIP0[iip];
        int ip1 = aIP1[npe-iip-1];
        aLine.push_back(ip0);
        aLine.push_back(ip1);
      }
    }
  }
  aXYZt = aXYZ;
  
  { // make a unit sphere
    MeshTri3D_Sphere(aXYZ_Contact, aTri_Contact, 0.3, 32, 32);
    Rotate(aXYZ_Contact, 0.2, 0.3, 0.4);
    aNorm_Contact.resize(aXYZ_Contact.size());
    Normal_MeshTri3D(aNorm_Contact.data(),
                     aXYZ_Contact.data(), aXYZ_Contact.size()/3,
                     aTri_Contact.data(), aTri_Contact.size()/3);
    bvh.Init(aXYZ_Contact.data(), aXYZ_Contact.size()/3,
             aTri_Contact.data(), aTri_Contact.size()/3,
             0.01);
  }

  
  win.camera.view_height = 1.0;
  win.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  setSomeLighting();
  // Enter main loop
  ::glutMainLoop();
  return 0;
}
