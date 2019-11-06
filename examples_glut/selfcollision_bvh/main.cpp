#include <iostream>
#include <vector>
#include "delfem2/bv.h"
#include "delfem2/bvh.h"
#include "delfem2/primitive.h"
#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"
#include "delfem2/srchbi_v3bvh.h"

// -----

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/opengl/gl2_funcs.h"
#include "../glut_cam.h"

/* ------------------------------------------------------------------------ */
// input parameter for simulation
std::vector<double> aXYZ0; // undeformed vertex positions
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aUVW; // vertex deformation modes
std::vector<unsigned int> aTri;  // index of triangles

// variables for self-collision
int iroot_bvh; // index BVH root node
std::vector<CNodeBVH> aNodeBVH; // array of BVH node
std::vector<CBV3D_Sphere> aBB_BVH; // array of AABB same size as aNodeBVH

std::vector<CIntersectTriPair> aITP;

// data for camera
bool is_animation;
double cur_time = 0;
CNav3D_GLUT nav;
int imode_draw = 0;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void myGlutDisplay(void)
{
  //	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  
  nav.SetGL_Camera();
  
  bool is_lighting = glIsEnabled(GL_LIGHTING);
  
  //  Draw_SurfaceMeshNorm(aXYZ, aTri, aNormal);
  delfem2::opengl::DrawMeshTri3D_Edge(aXYZ,aTri);
  
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(2);
  ::glColor3d(1,0,0);
  ::glBegin(GL_LINES);
  for(int iitp=0;iitp<aITP.size();++iitp){
    const CIntersectTriPair& itp = aITP[iitp];
    glVertex3d(itp.P[0].x, itp.P[0].y, itp.P[0].z);
    glVertex3d(itp.P[1].x, itp.P[1].y, itp.P[1].z);
  }
  ::glEnd();
  
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }
  else{              ::glDisable(GL_LIGHTING); }
  
  glColor3d(0,0,0);
  ShowFPS();
  
  ::glutSwapBuffers();
}

void myGlutIdle(){
  
  if( is_animation ){
    cur_time += 0.02;
    double d = sin(cur_time);
    for(int ip=0;ip<(int)aXYZ.size()/3;ip++){
      aXYZ[ip*3+0] =  aXYZ0[ip*3+0] + aUVW[ip*3+0]*d;
      aXYZ[ip*3+1] =  aXYZ0[ip*3+1] + aUVW[ip*3+1]*d;
      aXYZ[ip*3+2] =  aXYZ0[ip*3+2] + aUVW[ip*3+2]*d;
    }
    ////
    BVH_BuildBVHGeometry(iroot_bvh,
                         1.0e-5,
                         aXYZ.data(),aXYZ.size()/3,
                         aTri.data(),3,aTri.size()/3,
                         aNodeBVH,aBB_BVH);
    aITP.clear();
    GetIntersectTriPairs(aITP,
                         aXYZ,aTri,
                         iroot_bvh,
                         aNodeBVH,aBB_BVH); // output
    std::cout << aITP.size() << std::endl;
  }
  
  ::glutPostRedisplay();
}


void myGlutResize(int w, int h)
{
  glViewport(0, 0, w, h);
  ::glutPostRedisplay();
}

void myGlutSpecial(int Key, int x, int y)
{
  nav.glutSpecial(Key, x, y);
}

void myGlutMotion( int x, int y )
{
  nav.glutMotion(x, y);
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  nav.glutMouse(button, state, x, y);
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
    case 'd': // change draw mode
      imode_draw++;
      if( imode_draw >= 2 ){
        imode_draw = 0;
      }
      break;
    case 't':
      //      StepTime();
      break;
    case ' ':
      //      imode_contact++;
      aXYZ = aXYZ0;
      aUVW.assign(aUVW.size(),0.0);
      break;
  }
  ::glutPostRedisplay();
}

int main(int argc,char* argv[])
{
  {
    delfem2::MeshTri3D_Sphere(aXYZ0, aTri, 1.0, 16, 16);
    delfem2::Rotate(aXYZ0, 0.2, 0.3, 0.4);
    aXYZ = aXYZ0;
    {
      const int ntri = aTri.size()/3;
      std::vector<double> aElemCenter(ntri*3);
      for(int itri=0;itri<ntri;++itri){
        CVector3 p0 = cg_Tri(itri, aTri, aXYZ);
        aElemCenter[itri*3+0] = p0.x;
        aElemCenter[itri*3+1] = p0.y;
        aElemCenter[itri*3+2] = p0.z;
      }
      std::vector<int> aTriSurRel;
      makeSurroundingRelationship(aTriSurRel,
                                  aTri.data(), aTri.size()/3, 
                                  MESHELEM_TRI, aXYZ.size()/3);
      iroot_bvh = BVH_MakeTreeTopology(aNodeBVH,
                                       3,aTriSurRel,
                                       aElemCenter);
      std::cout << "aNodeBVH.size(): " << aNodeBVH.size() << std::endl;
    }
    //    aEdge.SetEdgeOfElem(aTri,(int)aTri.size()/3,3, aXYZ.size()/3,false);
  }
  {
    aUVW.assign(aXYZ.size(),0.0);
    for(int ixyz=0;ixyz<aXYZ.size()/3;++ixyz){
      double x0 = aXYZ[ixyz*3+0];
      aUVW[ixyz*3+0] = -3*x0*x0*x0*x0*x0;
    }
  }
  
  ///////////////////////////
  glutInit(&argc, argv);
  glutInitWindowPosition(200,200);
  glutInitWindowSize(400, 300);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  glutCreateWindow("3D View");
  glutDisplayFunc(myGlutDisplay);
  glutIdleFunc(myGlutIdle);
  glutReshapeFunc(myGlutResize);
  glutMotionFunc(myGlutMotion);
  glutMouseFunc(myGlutMouse);
  glutKeyboardFunc(myGlutKeyboard);
  glutSpecialFunc(myGlutSpecial);
  // ----------------
  
  delfem2::opengl::setSomeLighting();
  nav.camera.view_height = 1.5;
  
  glutMainLoop();
  return 0;
}


