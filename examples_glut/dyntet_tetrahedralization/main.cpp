#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <set>
#include <math.h>
#include <time.h>
#include "delfem2/vec3.h"
#include "delfem2/msh.h"
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"
#include "delfem2/dtet_v3.h"

// ----------------

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/opengl/gl2_color.h"
#include "../glut_cam.h"

// --------------------------------------------


void FindEdge
(int i0, int i1,
 std::vector<CEPo3D>& aPo3D,
 std::vector<CETet>& aSTet)
{
  int itet0 = aPo3D[i0].e;
  if( itet0 == -1 ) return;
  int inotet0 = aPo3D[i0].poel;
}

void Inactivate
(int it,
 std::vector<CETet>& aSTet)
{
  for(int ift=0;ift<4;++ift){
    int jt0 = aSTet[it].s[ift];
    if( jt0 == -1 ){
    }
    else{
      int jft0 = tetRel[ aSTet[it].f[ift] ][ift];
      assert( aSTet[jt0].s[jft0] == it );
      aSTet[jt0].s[jft0] = -1;
      aSTet[it].s[ift] = -1;
    }
  }
  aSTet[it].v[0] = -1;
  aSTet[it].v[1] = -1;
  aSTet[it].v[2] = -1;
  aSTet[it].v[3] = -1;
}

//////////////////////////////////////////////////////////////////////////////////////

std::vector<CEPo3D> aPo3D;
std::vector<CETet> aSTet;

int imode_display = 0;
bool is_animation = true;
CNav3D_GLUT nav;

int it_wrong = -1;

//////////////////////////////////////////////////////////////////////////////////////


static void myGlVertex3d(const CVector3& v)
{
  ::glVertex3d(v.x, v.y, v.z);
}

static void myGlVertex3d(int i, const std::vector<CVector3>& aV)
{
  const CVector3& v = aV[i];
  ::glVertex3d(v.x, v.y, v.z);
}


void myGlutDisplay(void)
{
  //	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);

  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);

  nav.SetGL_Camera();

  DrawBackground();

  ::glEnable(GL_BLEND);
  ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  ::glEnable(GL_CULL_FACE);
  ::glCullFace(GL_BACK);
  
  
  if (imode_display==1){
    ::glColor3d(0, 0, 0);
    ::glPointSize(10);
    ::glBegin(GL_POINTS);
    myGlVertex3d(aPo3D[631].p);
    myGlVertex3d(aPo3D[732].p);
    myGlVertex3d(aPo3D[964].p);
    ::glEnd();
  }
  
//  else if(imode_display ==1){
  {
    ::glColor3d(0, 0, 0);
    ::glPointSize(1);
    ::glBegin(GL_POINTS);
    for (int ip = 0; ip<(int)aPo3D.size(); ip++){
      glVertex3d(aPo3D[ip].p.x, aPo3D[ip].p.y, aPo3D[ip].p.z);
    }
    ::glEnd();
    ///
    ::glBegin(GL_LINES);
    for (int it = 0; it<aSTet.size(); ++it){
      int ip0 = aSTet[it].v[0];
      if( ip0 == -1 ) continue;
      int ip1 = aSTet[it].v[1];
      int ip2 = aSTet[it].v[2];
      int ip3 = aSTet[it].v[3];
      myGlVertex3d(aPo3D[ip0].p);
      myGlVertex3d(aPo3D[ip1].p);
      myGlVertex3d(aPo3D[ip0].p);
      myGlVertex3d(aPo3D[ip2].p);
      myGlVertex3d(aPo3D[ip0].p);
      myGlVertex3d(aPo3D[ip3].p);
      myGlVertex3d(aPo3D[ip1].p);
      myGlVertex3d(aPo3D[ip2].p);
      myGlVertex3d(aPo3D[ip1].p);
      myGlVertex3d(aPo3D[ip3].p);
      myGlVertex3d(aPo3D[ip2].p);
      myGlVertex3d(aPo3D[ip3].p);
    }
    ::glEnd();
  }
  
  ShowFPS();

  ::glutSwapBuffers();
}

void myGlutResize(int w, int h)
{
  ::glViewport(0, 0, w, h);
  ::glutPostRedisplay();
}

void myGlutMotion(int x, int y)
{
  nav.glutMotion(x, y);
  ::glutPostRedisplay();
}

void myGlutMouse(int ibutton, int state, int x, int y)
{
  nav.glutMouse(ibutton, state, x, y);
}

void Initialize()
{
  aPo3D.clear();
  aSTet.clear();
  {
    double len = 3.0;
    aPo3D.resize(4);
    aPo3D[0].p = CVector3(-len, +len, -len); aPo3D[0].e = 0; aPo3D[0].poel = 0;
    aPo3D[1].p = CVector3(+len, -len, -len); aPo3D[1].e = 0; aPo3D[1].poel = 1;
    aPo3D[2].p = CVector3(+len, +len, +len); aPo3D[2].e = 0; aPo3D[2].poel = 2;
    aPo3D[3].p = CVector3(-len, -len, +len); aPo3D[3].e = 0; aPo3D[3].poel = 3;
    aSTet.resize(1);
    aSTet[0].v[0] = 0;
    aSTet[0].v[1] = 1;
    aSTet[0].v[2] = 2;
    aSTet[0].v[3] = 3;
    aSTet[0].s[0] = -1;
    aSTet[0].s[1] = -1;
    aSTet[0].s[2] = -1;
    aSTet[0].s[3] = -1;
    aSTet[0].setCircumCenter(aPo3D);
  }

}

void AddRandomPoint()
{
  std::vector<int> tmp_buffer;
  for(int itr=0;itr<100;itr++){
    double x0 = 2.0*(double)rand()/(RAND_MAX+1.0)-1.0;
    double y0 = 2.0*(double)rand()/(RAND_MAX+1.0)-1.0;
    double z0 = 2.0*(double)rand()/(RAND_MAX+1.0)-1.0;
    int ip_ins = (int)aPo3D.size();
    aPo3D.resize(ip_ins+1);
    aPo3D[ip_ins].p = CVector3(x0,y0,z0);
    aPo3D[ip_ins].e = -1;
    aPo3D[ip_ins].poel = -1;
    int itet_ins = -1;
    {
      for (int it = 0; it<aSTet.size(); ++it){
        int j0 = aSTet[it].v[0];
        int j1 = aSTet[it].v[1];
        int j2 = aSTet[it].v[2];
        int j3 = aSTet[it].v[3];
        if( j0 == -1 ) continue; // floating tet
        double v0 = TetVolume(ip_ins, j1, j2, j3, aPo3D);
        double v1 = TetVolume(j0, ip_ins, j2, j3, aPo3D);
        double v2 = TetVolume(j0, j1, ip_ins, j3, aPo3D);
        double v3 = TetVolume(j0, j1, j2, ip_ins, aPo3D);
        //    double v4 = TetVolume(j0, j1, j2, j3, aPo3D);
        if (v0>0&&v1>0&&v2>0&&v3>0){ itet_ins = it; break; }
      }
    }
    if (itet_ins==-1){ return; }
    AddPointTetDelaunay(ip_ins,itet_ins, aPo3D, aSTet, tmp_buffer);
#ifndef NDEBUG
    CheckTet(aSTet, aPo3D);
#endif
  }
  std::cout << aPo3D.size() << std::endl;
}


void myGlutIdle()
{
  if (!is_animation){
    AddRandomPoint();
    ::glutPostRedisplay();
    return;
  }
  ::glutPostRedisplay();
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch (Key)
  {
    case 'q':
    case 'Q':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
      break;
    case '\033':
      break;
    case 'a':
      is_animation = !is_animation;
      break;
    case 'd':
      imode_display = (imode_display+1)%2;
      break;
    case ' ':
    {
      //      ConvexHull(aTri,aXYZ);
      it_wrong = -1;
      Initialize();
      break;
    }
    case 'v':
    {
      AddRandomPoint();
      break;
    }
    case '1':
    {
      std::vector<double> aXYZ;
      std::vector<unsigned int> aTri;
      Read_Obj("models/bunny2k.obj",aXYZ, aTri);
      Normalize(aXYZ);
      Scale(2.0,aXYZ);
      ////
      Initialize();
      std::vector<int> tmp_buffer;
      for(int ixyz=0;ixyz<aXYZ.size()/3;ixyz++){
        double x0 = aXYZ[ixyz*3+0];
        double y0 = aXYZ[ixyz*3+1];
        double z0 = aXYZ[ixyz*3+2];
        int ip_ins = (int)aPo3D.size();
        aPo3D.resize(ip_ins+1);
        aPo3D[ip_ins].p = CVector3(x0,y0,z0);
        aPo3D[ip_ins].e = -1;
        aPo3D[ip_ins].poel = -1;
        int itet_ins = -1;
        { // find tet
          for (int it = 0; it<aSTet.size(); ++it){
            int j0 = aSTet[it].v[0];
            int j1 = aSTet[it].v[1];
            int j2 = aSTet[it].v[2];
            int j3 = aSTet[it].v[3];
            if( j0 == -1 ) continue; // floating tet
            double v0 = TetVolume(ip_ins, j1, j2, j3, aPo3D);
            double v1 = TetVolume(j0, ip_ins, j2, j3, aPo3D);
            double v2 = TetVolume(j0, j1, ip_ins, j3, aPo3D);
            double v3 = TetVolume(j0, j1, j2, ip_ins, aPo3D);
            //    double v4 = TetVolume(j0, j1, j2, j3, aPo3D);
            if (v0>-0.0000001&&v1>-0.00000001&&v2>-0.000000001&&v3>-0.00000001){ itet_ins = it; break; }
          }
        }
        if (itet_ins==-1){ continue; }
        AddPointTetDelaunay(ip_ins,itet_ins, aPo3D, aSTet, tmp_buffer);
#ifndef NDEBUG
        CheckTet(aSTet, aPo3D);
#endif
      }
      { // remove SuperTet verteces
        for(int it=0;it<aSTet.size();++it){
          int i0 = aSTet[it].v[0];
          int i1 = aSTet[it].v[1];
          int i2 = aSTet[it].v[2];
          int i3 = aSTet[it].v[3];
          if( i0 < 4 || i1 < 4 || i2 < 4 || i3 < 4 ){
            Inactivate(it,aSTet);
          }
        }
        aPo3D[0].e = -1;
        aPo3D[1].e = -1;
        aPo3D[2].e = -1;
        aPo3D[3].e = -1;
        for(int it=0;it<aSTet.size();++it){
          for(int ift=0;ift<4;++ift){
            int iv = aSTet[it].v[ift];
            if( iv == -1 ) continue;
            aPo3D[iv].e = it;
            aPo3D[iv].poel = ift;
          }
        }
      }
      { // edge recovery
        std::vector<int> psup_ind;
        std::vector<int> psup;
        JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                                                    aTri.data(),aTri.size()/3,3,(int)aXYZ.size()/3);
        std::vector<int> edge_ind;
        std::vector<int> edge;
        JArrayEdgeUnidir_PointSurPoint(edge_ind, edge,
                                psup_ind, psup);
//        CJaggedArray edge;
//        edge.SetEdgeOfElem(aTri, (int)aTri.size()/3, 3, (int)aXYZ.size()/3, false);
        for(int ixyz=0;ixyz<(int)aXYZ.size()/3;++ixyz){
          int ip0 = ixyz+4;
          ElemAroundPoint elarpo;
          {
            int itet0 = aPo3D[ip0].e;
            if( itet0 == -1 ){
              std::cout << ip0 << " " << aPo3D[ip0].e << std::endl;
              continue;
            }
            MakeElemAroundPoint(elarpo, itet0, aPo3D[ip0].poel, aSTet);
          }
          for(int ie=edge_ind[ixyz];ie<edge_ind[ixyz+1];++ie){
            int jxyz = edge[ie];
            int jp0 = jxyz+3;
          }
        }
      }
      break;
    }

  }
  ::glutPostRedisplay();
}


void myGlutSpecial(int Key, int x, int y)
{
  nav.glutSpecial(Key, x, y);
  ::glutPostRedisplay();
}

int main(int argc, char* argv[])
{
  ::glutInit(&argc, argv);

  // Initialize GLUT
  glutInitWindowPosition(200, 200);
  glutInitWindowSize(400, 300);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  glutCreateWindow("Initial");

  // Define callback functions
  glutDisplayFunc(myGlutDisplay);
  glutReshapeFunc(myGlutResize);
  glutMotionFunc(myGlutMotion);
  glutMouseFunc(myGlutMouse);
  glutKeyboardFunc(myGlutKeyboard);
  glutSpecialFunc(myGlutSpecial);
  glutIdleFunc(myGlutIdle);

  nav.camera.view_height = 2.5;
  
  Initialize();
  AddRandomPoint();

  glutMainLoop();
  return 0;
}
