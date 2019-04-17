#include <iostream>
#include <vector>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"

#include "delfem2/matrix_sparse.h"
#include "delfem2/ilu_sparse.h"
#include "delfem2/fem.h"
#include "delfem2/cloth_internal.h"

#include "delfem2/funcs_gl.h"
#include "delfem2/color_gl.h"
#include "delfem2/funcs_glut.h"

/* ------------------------------------------------------------------------ */

// Setting problem here
void SetClothShape_Square
(std::vector<double>& aXYZ0, // (out) undeformed vertex positions
 std::vector<int>& aBCFlag, // (out) boundary condition flag (0:free 1:fixed)
 std::vector<int>& aTri, // (out) index of triangles
 std::vector<int>& aQuad, // (out) index of 4 vertices required for bending
 double& total_area, // (out) total area of cloth
 ///
 int ndiv, // (in) number of division of the square cloth edge
 double cloth_size) // (in) size of square cloth
{
  // make vertex potision array
  const double elem_length = cloth_size/ndiv; // size of an element
  const int nxyz =(ndiv+1)*(ndiv+1); // number of points
  aXYZ0.reserve( nxyz*3 );
  for(int ix=0;ix<ndiv+1;ix++){
    for(int iy=0;iy<ndiv+1;iy++){
      aXYZ0.push_back( ix*elem_length );
      aXYZ0.push_back( iy*elem_length );
      aXYZ0.push_back( 0.0 );
    }
  }
  
  // make triangle index array
  const int ntri = ndiv*ndiv*2;
  aTri.reserve( ntri*3 );
  for(int ix=0;ix<ndiv;ix++){
    for(int iy=0;iy<ndiv;iy++){
      aTri.push_back(  ix   *(ndiv+1)+ iy    );
      aTri.push_back( (ix+1)*(ndiv+1)+ iy    );      
      aTri.push_back(  ix   *(ndiv+1)+(iy+1) );
      ////
      aTri.push_back( (ix+1)*(ndiv+1)+(iy+1) );
      aTri.push_back(  ix   *(ndiv+1)+(iy+1) );
      aTri.push_back( (ix+1)*(ndiv+1)+ iy    );
    }
  }
  
  // make quad index array
  const int nquad = ndiv*ndiv + ndiv*(ndiv-1)*2;
  aQuad.reserve( nquad*4 );
  for(int ix=0;ix<ndiv;ix++){
    for(int iy=0;iy<ndiv;iy++){
      aQuad.push_back( (ix+0)*(ndiv+1)+(iy+0) );
      aQuad.push_back( (ix+1)*(ndiv+1)+(iy+1) );
      aQuad.push_back( (ix+1)*(ndiv+1)+(iy+0) );
      aQuad.push_back( (ix+0)*(ndiv+1)+(iy+1) );
    }
  }
  for(int ix=0;ix<ndiv;ix++){
    for(int iy=0;iy<ndiv-1;iy++){
      aQuad.push_back( (ix+1)*(ndiv+1)+(iy+0) );
      aQuad.push_back( (ix+0)*(ndiv+1)+(iy+2) );
      aQuad.push_back( (ix+1)*(ndiv+1)+(iy+1) );
      aQuad.push_back( (ix+0)*(ndiv+1)+(iy+1) );
    }
  }
  for(int ix=0;ix<ndiv-1;ix++){
    for(int iy=0;iy<ndiv;iy++){
      aQuad.push_back( (ix+0)*(ndiv+1)+(iy+1) );
      aQuad.push_back( (ix+2)*(ndiv+1)+(iy+0) );
      aQuad.push_back( (ix+1)*(ndiv+1)+(iy+0) );
      aQuad.push_back( (ix+1)*(ndiv+1)+(iy+1) );
    }
  }
  
  total_area = cloth_size*cloth_size;
  
  aBCFlag = std::vector<int>(nxyz*3,0);
  for(int iy=0;iy<ndiv+1;iy++){
    aBCFlag[iy*3+0] = 1;
    aBCFlag[iy*3+1] = 1;
    aBCFlag[iy*3+2] = 1;
  }
}



class CInput_ContactNothing: public CInput_Contact
{
public:
  double penetrationNormal(double& nx, double &ny, double& nz,
                           double px, double py, double pz) const
  {
    return -100;
  }
};


class CInput_ContactPlane: public CInput_Contact
{
  double penetrationNormal(double& nx, double &ny, double& nz,
                                double px, double py, double pz) const
  {
    nx = 0.0;  ny = 0.0;  nz = 1.0; // normal of the plane
    return -0.5 - pz; // penetration depth
  }
};

class CInput_ContactSphere: public CInput_Contact
{
  double penetrationNormal(double& nx, double &ny, double& nz,
                           double px, double py, double pz) const
  {
    const double center[3] = { 0.1, 0.5, -0.8 };
    const double radius = 0.3;
    nx = px-center[0];
    ny = py-center[1];
    nz = pz-center[2];
    double len = sqrt(nx*nx+ny*ny+nz*nz);
    nx /= len;
    ny /= len;
    nz /= len;
    return radius-len; // penetration depth
  }
};


/* ------------------------------------------------------------------------ */
// input parameter for simulation
const int ndiv = 25;  // (in) number of division of the square cloth edge
const double cloth_size = 1; // square cloth 1m x 1m
std::vector<double> aXYZ0; // undeformed vertex positions
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aUVW; // deformed vertex velocity
std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)
std::vector<int> aTri;  // index of triangles
std::vector<int> aQuad; // index of 4 vertices required for bending
const double lambda = 1.0; // Lame's 1st parameter
const double myu    = 4.0; // Lame's 2nd parameter
const double stiff_bend = 1.0e-3; // bending stiffness
const double areal_density = 1.0; // areal density of a cloth
const double gravity[3] = {0,0,-10}; // gravitatinal accereration
double time_step_size = 0.02; // size of time step
const double stiff_contact = 1.0e+3;
const double contact_clearance = 0.02;
int imode_contact; // mode of contacting object
double mass_point; // mass for a point

// variables for sparse solver
CMatrixSquareSparse mat_A; // coefficient matrix
CPreconditionerILU  ilu_A; // ilu decomposition of the coefficient matrix

std::vector<double> aNormal; // deformed vertex noamals

// data for camera
bool is_animation;
CGlutWindowManager win;
int imode_draw = 0;
/* ------------------------------------------------------------------------ */


void StepTime()
{
  CInput_ContactPlane c0;
  CInput_ContactNothing c1;
  CInput_ContactSphere c2;
  CInput_Contact* ic = 0;
  if(      imode_contact == 1 ){ ic = &c0; }
  if(      imode_contact == 0 ){ ic = &c1; }
  if(      imode_contact == 2 ){ ic = &c2; }
  // solving lienar system using conjugate gradient method with ILU(0) preconditioner
  StepTime_InternalDynamicsILU(aXYZ, aUVW, mat_A, ilu_A,
                               aXYZ0, aBCFlag,
                               aTri, aQuad,
                               time_step_size,
                               lambda, myu, stiff_bend,
                               gravity, mass_point,
                               stiff_contact,contact_clearance,*ic);
  MakeNormal(aNormal, aXYZ, aTri);
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
  /*
   ::glMatrixMode(GL_MODELVIEW);
   ::glLoadIdentity();
   ::glTranslated(camera_trans[0],camera_trans[1],camera_trans[2]);
   {
   double R_view_3d[16];
   QuatRot(R_view_3d, camera_qrot);
   ::glMultMatrixd(R_view_3d);
   }
   */
  bool is_lighting = glIsEnabled(GL_LIGHTING);
  
  { // draw triangle edge
    ::glDisable(GL_LIGHTING);    
    ::glLineWidth(1);
    ::glColor3d(0,0,0);
    ::glBegin(GL_LINES);
    for(int itri=0;itri<aTri.size()/3;itri++){
      const int ip0 = aTri[itri*3+0];
      const int ip1 = aTri[itri*3+1];
      const int ip2 = aTri[itri*3+2];
      double c[3][3] = {
        { aXYZ[ip0*3+0],aXYZ[ip0*3+1],aXYZ[ip0*3+2] },
        { aXYZ[ip1*3+0],aXYZ[ip1*3+1],aXYZ[ip1*3+2] },
        { aXYZ[ip2*3+0],aXYZ[ip2*3+1],aXYZ[ip2*3+2] } };
      ::glVertex3dv(c[0]);  ::glVertex3dv(c[1]);
      ::glVertex3dv(c[1]);  ::glVertex3dv(c[2]);
      ::glVertex3dv(c[2]);  ::glVertex3dv(c[0]);
    }
    ::glEnd();
  }
  {
    ::glEnable(GL_LIGHTING);    
    float color[4] = {200.0/256.0, 200.0/256.0, 200.0/256.0,1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
    ::glShadeModel(GL_SMOOTH);
    ::glBegin(GL_TRIANGLES);
    for(int itri=0;itri<aTri.size()/3;itri++){
      const int ip0 = aTri[itri*3+0];
      const int ip1 = aTri[itri*3+1];
      const int ip2 = aTri[itri*3+2];
      double c[3][3] = {
        { aXYZ[ip0*3+0],aXYZ[ip0*3+1],aXYZ[ip0*3+2] },
        { aXYZ[ip1*3+0],aXYZ[ip1*3+1],aXYZ[ip1*3+2] },
        { aXYZ[ip2*3+0],aXYZ[ip2*3+1],aXYZ[ip2*3+2] } };
      ::glNormal3d(aNormal[ip0*3+0],aNormal[ip0*3+1],aNormal[ip0*3+2]);
      ::glVertex3dv(c[0]);
      ::glNormal3d(aNormal[ip1*3+0],aNormal[ip1*3+1],aNormal[ip1*3+2]);
      ::glVertex3dv(c[1]);
      ::glNormal3d(aNormal[ip2*3+0],aNormal[ip2*3+1],aNormal[ip2*3+2]);
      ::glVertex3dv(c[2]);
    }
    ::glEnd();    
  }
  
  { // fixed boundary condition
    ::glDisable(GL_LIGHTING);        
    ::glPointSize(5);
    ::glColor3d(0,0,1);
    ::glBegin(GL_POINTS);
    for(int ip=0;ip<aXYZ.size()/3;ip++){
      if( aBCFlag[ip*3+0] == 0 && aBCFlag[ip*3+1] == 0 && aBCFlag[ip*3+2] == 0 ) continue;
      ::glVertex3d(aXYZ0[ip*3+0],aXYZ0[ip*3+1],aXYZ0[ip*3+2]);
    }
    ::glEnd();
  }
  
  if(      imode_contact == 1 ){     // draw floor
    ::glDisable(GL_LIGHTING);        
    ::glLineWidth(1);
    ::glColor3d(1,0,0);
    ::glBegin(GL_LINES);
    double grid_x_min = -10;
    double grid_x_max = +10;
    double grid_y_min = -10;
    double grid_y_max = +10;
    int ndiv_grid = 30;
    for(unsigned int ix=0;ix<ndiv_grid+1;ix++){
      double x0 = (grid_x_max-grid_x_min) / ndiv_grid * ix + grid_x_min;
      ::glVertex3d(x0,grid_y_min,-0.5);
      ::glVertex3d(x0,grid_y_max,-0.5);
    }
    for(unsigned int iz=0;iz<ndiv_grid+1;iz++){
      double z0 = (grid_y_max-grid_y_min) / ndiv_grid * iz + grid_y_min;
      ::glVertex3d(grid_x_min,z0,-0.5);
      ::glVertex3d(grid_x_max,z0,-0.5);
    }
    ::glEnd();        
  }
  else if( imode_contact == 2 ){
    ::glDisable(GL_LIGHTING);        
    ::glLineWidth(1);
    ::glColor3d(1,0,0);    
    ::glPushMatrix();
    ::glTranslated(0.1, 0.5, -0.8);
    ::glutWireSphere(0.3, 16, 16);
    ::glPopMatrix();
  }
  
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }
  else{              ::glDisable(GL_LIGHTING); }  
  
  ::glutSwapBuffers();
}

void myGlutIdle(){
  
  if( is_animation ){
    StepTime();
  }
  
  ::glutPostRedisplay();
}


void myGlutResize(int w, int h)
{
  glViewport(0,0, w, h);
  ::glutPostRedisplay();
}

void myGlutSpecial(int Key, int x, int y)
{
  win.glutSpecial(Key, x, y);
  ::glutPostRedisplay();
}


void myGlutMotion( int x, int y ){
  win.glutMotion(x, y);
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y){
  win.glutMouse(button, state, x, y);
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
      StepTime();
      break;
    case ' ':
      imode_contact++;
      aXYZ = aXYZ0;
      aUVW.assign(aUVW.size(),0.0);
      if( imode_contact >= 3 ){
        imode_contact = 0;
      }
      break;
  }
  ::glutPostRedisplay();
}

int main(int argc,char* argv[])
{
  { // initialze data
    double total_area;
    SetClothShape_Square(aXYZ0,aBCFlag,aTri,aQuad,total_area,
                         ndiv,cloth_size);
    const int np = aXYZ0.size()/3.0;
    mass_point = total_area*areal_density / (double)np;
    // initialize deformation
    aXYZ = aXYZ0;
    aUVW.assign(np*3,0.0);
    MakeNormal(aNormal, aXYZ, aTri);
    mat_A.Initialize(np,3,true);
    std::vector<int> psup_ind,psup;
    JArray_MeshOneRingNeighborhood(psup_ind, psup,
                                   aQuad.data(),aQuad.size()/4, 4, np);
    JArray_Sort(psup_ind, psup);
    mat_A.SetPattern(psup_ind.data(),psup_ind.size(), psup.data(),psup.size());
    ilu_A.Initialize_ILU0(mat_A);
  }
  
  
  glutInit(&argc, argv);
  
  // Initialize GLUT window 3D
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
  
  ////////////////////////
  
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);
  setSomeLighting();
  win.camera.camera_rot_mode = CAMERA_ROT_ZTOP;
  win.camera.camera_rot_mode = CAMERA_ROT_ZTOP;
  win.camera.psi = 3.1415*0.2;
  win.camera.theta = 3.1415*0.1;
  win.camera.view_height = 2;
  
  
  glutMainLoop();
  return 0;
}


