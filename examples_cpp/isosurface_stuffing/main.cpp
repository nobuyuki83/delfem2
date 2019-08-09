#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <fstream>
#include <time.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/glut_funcs.h"
#include "delfem2/gl_color.h"
#include "delfem2/gl_funcs.h"

#include "delfem2/isosurface_stuffing.h"
#include "delfem2/sdf.h"
#include "delfem2/msh.h"
#include "delfem2/mshio.h"
#include "delfem2/mshtopo.h"

//////////////////////////////////////////////////////////////


std::vector<double> aXYZ;
std::vector<unsigned int> aTet;
std::vector<unsigned int> aTetSurface;
std::vector<CColor> aTetColor;

std::vector<unsigned int> aTet1;
std::vector<CColor> aTetColor1;

CGlutWindowManager win;
double vis_cut_org[3] = {-0.0, 0.0, 0.0};
double vis_cut_nrm[3] = {0.0,-0.9, +0.2};
//double vis_cut_time = 0;
bool is_animation = false;
double cur_time = 0.5;
int imode_draw = 0;


void SetProblem()
{
  const unsigned int nprob = 4;
  static int iprob = 0;
  
  std::vector<int> aIsOnSurfXYZ;
  if( iprob == 0 ){
    class CSphere : public CInputIsosurfaceStuffing
    {
    public:
      CSphere(double rad){
        sp.radius_ = rad;
      }
      virtual double SignedDistance(double px, double py, double pz) const{
        double n[3];
        return sp.Projection(px,py,pz,n);
      }
      virtual void Level(int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
                         double px, double py, double pz) const
      {
        sdf = this->SignedDistance(px,py,pz);
        const double rad0 = sp.radius_;
        const double rad1 = sqrt(px*px+py*py+pz*pz);
        if( rad1 > rad0*0.5 ){ ilevel_vol = 0; }
        else{ ilevel_vol = 1; }
        ilevel_srf = 1;
        nlayer = 1;
        ilevel_vol = -1;
      }
    public:
      CSDF3_Sphere sp;
    };
    double rad = 1.5;
    CSphere sphere(rad);
    double cent[3] = {0,0,0};
    IsoSurfaceStuffing(aXYZ, aTet,aIsOnSurfXYZ,
                       sphere, 0.2, rad*4.0, cent);
  }
  else if( iprob ==  1 ){
    class CBox : public CInputIsosurfaceStuffing
    {
    public:
      CBox(double hwx, double hwy, double hwz){
        bx.hwx = hwx;
        bx.hwy = hwy;
        bx.hwz = hwz;
      }
      virtual double SignedDistance(double px, double py, double pz) const{
        double n[3];
        return bx.Projection(px,py,pz,n);
      }
      virtual void Level(int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
                         double px, double py, double pz) const
      {
        sdf = this->SignedDistance(px,py,pz);
        ilevel_vol = -1;
        ilevel_srf = 1;
        nlayer = 1;
      }
    public:
      CSDF3_Box bx;
    };
    const double hwx = 0.91;
    const double hwy = 0.61;
    const double hwz = 0.41;
    double cent[3] = {0,0,0};
    CBox box(hwx,hwy,hwz);
    IsoSurfaceStuffing(aXYZ, aTet,aIsOnSurfXYZ,
                       box, 0.1, 2.00, cent);
  }
  else if( iprob == 2 ){
    class CCavSphere : public CInputIsosurfaceStuffing
    {
    public:
      CCavSphere(){
        const double hwx = 0.91;
        const double hwy = 0.61;
        const double hwz = 0.61;
        box.hwx = hwx;
        box.hwy = hwy;
        box.hwz = hwz;
        ////
        const double rad = 0.2;
        sphere.radius_ = rad;
      }
      virtual double SignedDistance(double x, double y, double z ) const {
        double n[3];
        double dist0 = -sphere.Projection(x, y, z, n);
        double cx = 0.0;
        double cy = 0.0;
        double cz = 0.0;
        double dist1 = box.Projection(x-cx, y-cy, z-cz, n);
        return (dist0<dist1) ? dist0 : dist1;
      }
      virtual void Level(int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
                         double px, double py, double pz) const
      {
        sdf = this->SignedDistance(px,py,pz);
        double dist0 = sqrt(px*px+py*py+pz*pz);
        ilevel_vol = -1;
        if( dist0 < 0.5 && px > 0 ){ ilevel_srf = +3; }
        else{                        ilevel_srf = -1; }
        nlayer = 1;
      }
    public:
      CSDF3_Box box;
      CSDF3_Sphere sphere;
    } cav_sphere;
    double cent[3] = {0,0,0};
    IsoSurfaceStuffing(aXYZ, aTet, aIsOnSurfXYZ,
                       cav_sphere, 0.2, 3.00, cent);
  }
  if( iprob == 3 ){
    class CMesh : public CInputIsosurfaceStuffing
    {
    public:
      virtual double SignedDistance(double x, double y, double z) const {
        double n[3];
        return sdf_mesh.Projection(x, y, z,n);
      }
      virtual void Level(int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
                         double px, double py, double pz) const
      {
        sdf = this->SignedDistance(px,py,pz);
        ilevel_vol = 0;
        ilevel_srf = 2;
        nlayer = 2;
      }
    public:
      CSDF3_Mesh sdf_mesh;
    } mesh;
    {
      std::vector<unsigned int> aTri;
      std::vector<double> aXYZ_Tri;
      Read_Ply(std::string(PATH_INPUT_DIR)+"/bunny_1k.ply", aXYZ_Tri, aTri);
      Normalize(aXYZ_Tri,2.3);
      mesh.sdf_mesh.SetMesh(aTri, aXYZ_Tri);
      mesh.sdf_mesh.BuildBoxel();
    }
    double cent[3] = {0,0,0};
    IsoSurfaceStuffing(aXYZ, aTet, aIsOnSurfXYZ,
                       mesh, 0.18, 3.0, cent);
  }
  /*
  if( iprob == 4 ){
    class CCavMesh : public CInputIsosurfaceStuffing
    {
    public:
      CCavMesh(){
      }
      virtual double SignedDistance(double x, double y, double z ) const {
        double ntmp[3];
        double dist0 = -sdf.Projection(x, y, z, ntmp);
        double cx = 0.0;
        double cy = 1.0;
        double cz = 0.0;
        double dist1 = box.Projection(x-cx, y-cy, z-cz, ntmp);
        return (dist0<dist1) ? dist0 : dist1;
      }
      virtual void Level(int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
                         double px, double py, double pz) const {
        sdf = this->SignedDistance(px,py,pz);
        const double lx = box.hwx;
        const double ly = box.hwy;
        const double lz = box.hwz;
        if( fabs(px)<lx*0.6 && fabs(py)<ly*0.5 && fabs(pz)<lz*0.5 ){ ilevel_srf = 4; }
        else{ ilevel_srf = 0; }
        if( (px>0&&px<lx*0.7) && fabs(py)<ly*0.5 && fabs(pz)<lz*0.5 ){ ilevel_vol = 1; }
        else{ ilevel_vol = 0; }
        ilevel_vol = -1;
        nlayer = 2;
      }
    public:
      CSDF3_Box box;
      CSDF3_Mesh sdf;
    } cav_mesh;
    ////
    double lenx = 3;
    double leny = 3;
    double lenz = 6;
    std::vector<unsigned int> aTriIn;
    std::vector<double> aXYZIn;
    Read_Ply("models/car.ply", aXYZIn, aTriIn);
    {
      cav_mesh.sdf.SetMesh(aTriIn, aXYZIn);
      cav_mesh.sdf.BuildBoxel();
      cav_mesh.box.hwx = lenx;
      cav_mesh.box.hwy = leny;
      cav_mesh.box.hwz = lenz;
    }
    /////
    double resolution = 1.1;
    double centre[3] = {0,0,0};
    IsoSurfaceStuffing(aXYZ, aTet,aIsOnSurfXYZ,
                       cav_mesh, resolution, lenz*2.3, centre);
  }
   */
  
  aTetColor.resize(aTet.size()/4);
  for(int it=0;it<aTet.size()/4;++it){
    aTetColor[it].setRandomVividColor();
  }
  
  std::vector<int> aTetSurRel;
  makeSurroundingRelationship(aTetSurRel,
                              aTet.data(), aTet.size()/4,
                              MESHELEM_TET,
                              aXYZ.size()/3);
  aTetSurface.clear();
  for(int it=0;it<aTet.size()/4;++it){
    for(int ift=0;ift<4;++ift){
      if( aTetSurRel[it*8+ift*2+0] != -1 ){ continue; }
      aTetSurface.push_back(it);
      aTetSurface.push_back(ift);
    }
  }
  
  aTet1.clear();
  aTetColor1.clear();
  
  iprob++;
  if( iprob == nprob ){ iprob = 0; }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void myGlutResize(int w, int h)
{
  ::glViewport(0, 0, w, h);
	::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
  win.glutMotion(x,y);
}

void myGlutMouse(int button, int state, int x, int y){
  win.glutMouse(button,state,x,y);
}

void SetProblem();
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
		case 'b':
			cur_time = 0;
			break;
    case 'd':
      imode_draw = (imode_draw+1)%3;
      break;
		case ' ':	
			SetProblem();
			break;
	}
}


void myGlutSpecial(int Key, int x, int y)
{
  win.glutSpecial(Key,x,y);
}

void myGlutIdle(){
  if( is_animation ){
    cur_time += 0.005;
    if( cur_time > 1 ){ cur_time = 0.0; }
    vis_cut_org[1] = -1*(1-cur_time) + 1*cur_time;
  }
	::glutPostRedisplay();
}

void myGlutDisplay(void)
{
//	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
	::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);  
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);
	
	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );
  
  win.SetGL_Camera();

  ::glEnable(GL_LIGHTING);
  if( imode_draw == 0 ){
    DrawMeshTet3D_Cut(aXYZ,aTet,aTetColor,
                      vis_cut_org, vis_cut_nrm);
  }
  else if( imode_draw == 1 ){
    DrawMeshTet3DSurface_Edge(aXYZ, aTet, aTetSurface);
  }
  else if( imode_draw == 2 ){
    DrawMeshTet3D_Cut(aXYZ,aTet1,aTetColor1,
                      vis_cut_org, vis_cut_nrm);
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
	
//  Hoge();
  SetProblem();
  win.camera.view_height = 2;
  win.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  ////
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_LIGHT0);            
  float light0pos[4] = {0.5,0.5,+20,0};
  ::glLightfv(GL_LIGHT0, GL_POSITION, light0pos);      
  float white[3] = {1.0,1.0,1.0};
  ::glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
  
  ::glEnable(GL_LIGHT1);    
  float light1pos[4] = {0.5,0.5,-20,0}; 
  ::glLightfv(GL_LIGHT1, GL_POSITION, light1pos);          
  float white2[3] = {1.0,1.0,0.0};
  ::glLightfv(GL_LIGHT1, GL_DIFFUSE, white2);        
  ////
  
  
	glutMainLoop();
	return 0;
}

