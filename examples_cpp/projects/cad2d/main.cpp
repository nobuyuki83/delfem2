//
//  main.cpp
//  rigid_joint
//
//  Created by Nobuyuki Umetani on 11/3/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <math.h>

#include <complex>
#include <set>
#include <stack>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/mat3.h"

#include "delfem2/dyntri_v3.h"

#include "delfem2/funcs_glut.h"
#include "delfem2/funcs_gl.h"
#include "delfem2/v23_gl.h"

#ifndef M_PI
#define M_PI 3.141592653589793
#endif


class CCadTopo
{
public:
  void Clear(){
    nVertex = 0;
    aEdge.clear();
    aFace.clear();
  }
  void AddPolygon(int np){
    const int iv0 = nVertex;
    nVertex += np;
    const int ie0 = aEdge.size();
    for(int iie=0;iie<np;++iie){
      CEdge edge0;
      edge0.iv0 = iv0 + (iie+0)%np;
      edge0.iv1 = iv0 + (iie+1)%np;
      aEdge.push_back(edge0);
    }
    CFace face0;
    for(int iie=0;iie<np;++iie){
      face0.aIE.push_back( std::make_pair(ie0+iie,true ) );
    }
    aFace.push_back(face0);
  }
public:
  class CEdge{
  public:
    int iv0,iv1;
  };
  class CFace{
  public:
    std::vector< std::pair<int,bool> > aIE; // index of edge, is this edge ccw?
  };
public:
  int nVertex;
  std::vector<CEdge> aEdge;
  std::vector<CFace> aFace;
};


//////////////////

class CCad2D
{
public:
  void Clear(){
    aVtx.clear();
    aEdge.clear();
    aFace.clear();
    aTopo.Clear();
  }
  void Initialize_Square(){
    Clear();
    aTopo.AddPolygon(4);
    aVtx.push_back(CVector2(-1.0, -1.0));
    aVtx.push_back(CVector2(+1.0, -1.0));
    aVtx.push_back(CVector2(+1.0, +1.0));
    aVtx.push_back(CVector2(-1.0, +1.0));
    int iedge0 = aEdge.size();
    aEdge.push_back(CEdgeGeo());
    aEdge.push_back(CEdgeGeo());
    aEdge.push_back(CEdgeGeo());
    aEdge.push_back(CEdgeGeo());
    int iface0 = aFace.size();
    aFace.push_back(CFaceGeo());
    ////
    aEdge[iedge0+0].GenMesh(iedge0+0,aTopo,aVtx);
    aEdge[iedge0+1].GenMesh(iedge0+1,aTopo,aVtx);
    aEdge[iedge0+2].GenMesh(iedge0+2,aTopo,aVtx);
    aEdge[iedge0+3].GenMesh(iedge0+3,aTopo,aVtx);
    aFace[iface0].GenMesh(aTopo,iface0,aEdge);
  }
  void Draw() const {
    ::glDisable(GL_LIGHTING);
    ::glColor3d(1,0,0);
    ::glPointSize(6);
    ::glBegin(GL_POINTS);
    for(int iv=0;iv<aVtx.size();++iv){
      ::glVertex3d( aVtx[iv].pos.x, aVtx[iv].pos.y, 0.0);
    }
    ::glEnd();
    /////
    ::glColor3d(0,0,0);
    ::glLineWidth(3);
    ::glBegin(GL_LINES);
    for(int ie=0;ie<aEdge.size();++ie){
      int iv0 = aTopo.aEdge[ie].iv0;
      int iv1 = aTopo.aEdge[ie].iv1;
      ::glVertex3d( aVtx[iv0].pos.x, aVtx[iv0].pos.y, -0.1);
      ::glVertex3d( aVtx[iv1].pos.x, aVtx[iv1].pos.y, -0.1);
    }
    ::glEnd();
    //////
    ::glColor3d(0.8,0.8,0.8);
    ::glLineWidth(1);
    glTranslated(0,0,-0.2);
    for(int iface=0;iface<aFace.size();++iface){
      const CFaceGeo& face = aFace[iface];
      DrawMeshTri2D_Face(face.aTri, face.aXY);
      DrawMeshTri2D_Edge(face.aTri, face.aXY);
    }
    glTranslated(0,0,+0.2);
  }
public:
public:
  class CVertexGeo{
  public:
    CVertexGeo(const CVector2& p) : pos(p){}
  public:
    CVector2 pos;
  };
  class CEdgeGeo{
  public:
    void GenMesh(int iedge, const CCadTopo& topo,
                 std::vector<CVertexGeo>& aVtxGeo)
    {
      assert( iedge>=0 && iedge<topo.aEdge.size() );
      const int iv0 = topo.aEdge[iedge].iv0;
      const int iv1 = topo.aEdge[iedge].iv1;
      this->p0 = aVtxGeo[iv0].pos;
      this->p1 = aVtxGeo[iv1].pos;
    }
  public:
    CVector2 p0,p1;
    std::vector<CVector2> aP;
  };
  class CFaceGeo{
  public:
//    std::vector<CVector2> aP;
    std::vector<int> aTri;
    std::vector<double> aXY;
  public:
    void GenMesh(const CCadTopo& topo,int iface0,
                 std::vector<CEdgeGeo>& aEdgeGeo){
      assert( iface0>=0 && iface0<topo.aFace.size() );
      const std::vector< std::pair<int,bool> >& aIE = topo.aFace[iface0].aIE;
      std::vector<double> aXY_corner;
      for(int iie=0;iie<aIE.size();++iie){
        int ie0 = aIE[iie].first;
        assert( ie0>=0 && ie0<topo.aEdge.size() );
        const bool dir0 = aIE[iie].second;
        int iv0 = (dir0) ? topo.aEdge[ie0].iv0 : topo.aEdge[ie0].iv1;
        {
          const CEdgeGeo& eg0 = aEdgeGeo[ie0];
          CVector2 p0 = (dir0) ? eg0.p0 : eg0.p1;
          aXY_corner.push_back(p0.x);
          aXY_corner.push_back(p0.y);
        }
      }
      for(int ixy=0;ixy<aXY_corner.size()/2;++ixy){
        std::cout << aXY_corner[ixy*2+0] << " " << aXY_corner[ixy*2+1] << std::endl;
      }
      {
        std::vector<int> aPtrVtxInd,aVtxInd;
        std::vector< std::vector<double> > aVecAry0;
        aVecAry0.push_back(aXY_corner);
        GenerateTesselation2(aTri, aXY,  aPtrVtxInd,aVtxInd,
                             -1, false, aVecAry0);
      }
      std::cout << aTri.size() << std::endl;
    }
  };
public:
  CCadTopo aTopo;
  /////
  std::vector<CVertexGeo> aVtx;
  std::vector<CEdgeGeo> aEdge;
  std::vector<CFaceGeo> aFace;
};

///////////////////////////////////////////////////////////////////////////////////////////////


CGlutWindowManager win;
const double view_height = 2.0;
bool is_animation = false;
int imode_draw = 0;

CCad2D cad;

////////////////////////////////////////////////////////////////////////////////////////////////

void myGlutDisplay(void)
{
//  	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
	::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
//	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );
  
  win.SetGL_Camera();
  
//  if(      imode_draw == 0 ){ cad.DrawFace_RightSelected(false); }
//  else if( imode_draw == 1 ){ cad.DrawFace_RightSelected(true); }
//  else if( imode_draw == 2 ){ cad.DrawFace_LeftRight(); }
//  cad.DrawVtxEdgeHandler(win.camera.view_height);
  
  cad.Draw();
  
  ::glColor3d(0,0,0);
  ShowFPS();
  
  ::glutSwapBuffers();
}

void myGlutIdle(){
  if( is_animation ){
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
  win.glutSpecial(Key,x,y);
}

void myGlutMotion( int x, int y ){
  win.glutMotion(x,y);
  if( win.imodifier != 0){ return; }
  float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
  float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
  CVector2 sp0(win.mouse_x-win.dx, win.mouse_y-win.dy);
  CVector2 sp1(win.mouse_x, win.mouse_y);
  const CVector3 src_pick = screenUnProjection(CVector3(sp0.x,sp0.y, 0.0), mMV,mPj);
  const CVector3 dir_pick = screenUnProjectionDirection(CVector3(0.0,  0, -1.0 ), mMV,mPj);
  /////
//  cad.MouseMotion(src_pick,dir_pick, sp0,sp1, mMV, mPj);
}

void myGlutMouse(int button, int state, int x, int y)
{
  win.glutMouse(button,state,x,y);
  if( win.imodifier == GLUT_ACTIVE_SHIFT || win.imodifier == GLUT_ACTIVE_ALT ) return;
  float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
  float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
  CVector2 sp0(win.mouse_x, win.mouse_y);
  const CVector3 src_pick = screenUnProjection(CVector3(sp0.x,sp0.y, 0.0), mMV,mPj);
  const CVector3 dir_pick = screenUnProjectionDirection(CVector3(0.0,  0, -1.0 ), mMV,mPj);
  if( state == GLUT_DOWN ){
//    cad.MouseDown(src_pick, dir_pick, sp0, mMV, mPj, view_height);
  }
  if( state == GLUT_UP ){
//    cad.MouseUp(mMV,mPj,view_height);
  }
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
	switch(Key)
	{
    case 'q':
    case 'Q':
    case '\033':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
      break;
    case 'a':
    {
      is_animation = !is_animation;
      break;
    }
    case 'd':
    {
      imode_draw = (imode_draw+1)%3;
      break;
    }
    case 'b':
    {
      break;
    }
    case 'p':
    {
      break;
    }
    case 'e':
    {
//      MakeItSmooth(cad.aVertex, cad.aEdge, cad.aFace);
      break;
    }
    case 'f':
    {
      break;
    }
    case 's':
    {
//      cad.Initialize_Sphere();
      break;
    }
    case 'c':
    {
//      cad.Initialize_Cube();
      break;
    }
    case 'n':
    {
      break;
    }
    case 't':
    {
      break;
    }
      /*
    case 'w':
    {
      std::ofstream fout;
      fout.open("hoge.txt",std::ios::out);
      cad.WriteFile(fout);
      break;
    }
    case 'r':
    {
      std::ifstream fin;
      fin.open("hoge.txt",std::ios::in);
      cad.ReadFile(fin);
      break;
    }
    case '1':
    {
      cad.imode_edit = CCad3D::EDIT_MOVE;
      break;
    }
    case '2':
    {
      cad.imode_edit = CCad3D::EDIT_ADD_CROSS_SECTION;
      break;
    }
    case '3':
    {
      cad.imode_edit = CCad3D::EDIT_ADD_POINT_EDGE;
      break;
    }
    case '4':
    {
      cad.imode_edit = CCad3D::EDIT_SKETCH;
      break;
    }
    case '5':
    {
      break;
    }
    case '+':
    {
      cad.ChangeEdgeLength(cad.elen*0.9);
      break;
    }
    case '-':
    {
      cad.ChangeEdgeLength(cad.elen/0.9);
      break;
    }
       */
  }
  ::glutPostRedisplay();
}

int main(int argc,char* argv[])
{
  glutInit(&argc, argv);
  
	// Initialize GLUT window 3D
  glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
// 	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_STENCIL);
  glutCreateWindow("3D View");
	glutDisplayFunc(myGlutDisplay);
	glutIdleFunc(myGlutIdle);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
  
  ////////////////////////
  win.camera.view_height = view_height;
//  win.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  win.camera.camera_rot_mode = CAMERA_ROT_YTOP;
//    win.camera.camera_rot_mode = CAMERA_ROT_ZTOP;
  
  setSomeLighting();
  
  cad.Initialize_Square();
  
  glutMainLoop();
	return 0;
}


