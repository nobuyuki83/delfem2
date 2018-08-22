#if defined(USE_GL)
  #if defined(__APPLE__) && defined(__MACH__)
    #include <GLUT/glut.h>
  #else
    #include <GL/glut.h>
  #endif
#endif

#include "delfem2/vec3.h"
#include "delfem2/voxel.h"

const unsigned int noelElemFace_Vox[6][4] = {
  { 0, 4, 6, 2 }, // -x
  { 1, 3, 7, 5 }, // +x
  { 0, 1, 5, 4 }, // -y
  { 2, 6, 7, 3 }, // +y
  { 0, 2, 3, 1 }, // -z
  { 4, 5, 7, 6 } }; // +z

const CVector3 normalHexFace[6] = {
  CVector3(-1, 0, 0),
  CVector3(+1, 0, 0),
  CVector3( 0,-1, 0),
  CVector3( 0,+1, 0),
  CVector3( 0, 0,-1),
  CVector3( 0, 0,+1)
};

//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USE_GL

static void myGlNormal(const CVector3& n){ ::glNormal3d(n.x,n.y,n.z); }
static void myGlVertex(const CVector3& v){ ::glVertex3d(v.x,v.y,v.z); }
static void myGlColorDiffuse(float r, float g, float b, float a){
  ::glColor4d(r, g, b, a );
  float c[4] = {r, g, b, a};
  ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);
}

void Draw_CubeGrid
(bool is_picked, int iface_picked,
 double elen, const CVector3& org,
 const CCubeGrid& cube)
{
  if( !cube.is_active ) return;
  int ih = cube.ivx;
  int jh = cube.ivy;
  int kh = cube.ivz;
  CVector3 aP[8] = {
    org + elen*CVector3(ih+0,jh+0,kh+0),
    org + elen*CVector3(ih+1,jh+0,kh+0),
    org + elen*CVector3(ih+0,jh+1,kh+0),
    org + elen*CVector3(ih+1,jh+1,kh+0),
    org + elen*CVector3(ih+0,jh+0,kh+1),
    org + elen*CVector3(ih+1,jh+0,kh+1),
    org + elen*CVector3(ih+0,jh+1,kh+1),
    org + elen*CVector3(ih+1,jh+1,kh+1) };
  ::glEnable(GL_LIGHTING);
  ::glBegin(GL_QUADS);
  for(int iface=0;iface<6;++iface){
    if( is_picked && iface_picked == iface ){
      ::myGlColorDiffuse(1,1,0,1);
    }
    else{
      ::myGlColorDiffuse(0.8,0.8,0.8,1);
    }
    myGlNormal(normalHexFace[iface]);
    myGlVertex(aP[noelElemFace_Vox[iface][0]]);
    myGlVertex(aP[noelElemFace_Vox[iface][1]]);
    myGlVertex(aP[noelElemFace_Vox[iface][2]]);
    myGlVertex(aP[noelElemFace_Vox[iface][3]]);
  }
  ::glEnd();
  //////////////////////////////
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINE_STRIP);
  myGlVertex(aP[0]);
  myGlVertex(aP[1]);
  myGlVertex(aP[3]);
  myGlVertex(aP[2]);
  myGlVertex(aP[0]);
  myGlVertex(aP[4]);
  myGlVertex(aP[5]);
  myGlVertex(aP[7]);
  myGlVertex(aP[6]);
  myGlVertex(aP[4]);
  glEnd();
  ::glBegin(GL_LINES);
  myGlVertex(aP[1]); myGlVertex(aP[5]);
  myGlVertex(aP[2]); myGlVertex(aP[6]);
  myGlVertex(aP[3]); myGlVertex(aP[7]);
  ::glEnd();
}
#elif
void Draw_CubeGrid
(bool is_picked, int iface_picked,
 double elen, const CVector3& org,
 const CCubeGrid& cube)
{
  std::cout << "Error!->define #USE_GL to use this funciton" << std::endl;
}
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////

int Adj_Grid
(int igridvox, int iface,
 int ndivx, int ndivy, int ndivz)
{
  int ivx0 = igridvox/(ndivy*ndivz);
  int ivy0 = (igridvox-ivx0*(ndivy*ndivz))/ndivz;
  int ivz0 = igridvox-ivx0*(ndivy*ndivz)-ivy0*ndivz;
  if( iface == 0 ){ ivx0 -= 1; }
  if( iface == 1 ){ ivx0 += 1; }
  if( iface == 2 ){ ivy0 -= 1; }
  if( iface == 3 ){ ivy0 += 1; }
  if( iface == 4 ){ ivz0 -= 1; }
  if( iface == 5 ){ ivz0 += 1; }
  if( ivx0 < 0 || ivx0 >= ndivx ){ return -1; }
  if( ivy0 < 0 || ivy0 >= ndivy ){ return -1; }
  if( ivz0 < 0 || ivz0 >= ndivz ){ return -1; }
  int igv1 = ivx0*(ndivy*ndivz)+ivy0*ndivz+ivz0;
  assert( igv1 >= 0 && igv1 < ndivx*ndivy*ndivz );
  return igv1;
}

void Pick_CubeGrid
(int& icube_pic, int& iface_pic,
 const CVector3& src_pic, const CVector3& dir_pic,
 double elen,
 const CVector3& org,
 const std::vector<CCubeGrid>& aCube)
{
  icube_pic = -1;
  double depth_min = 0;
  for(int ivox=0;ivox<aCube.size();++ivox){
    if( !aCube[ivox].is_active ) continue;
    int ih = aCube[ivox].ivx;
    int jh = aCube[ivox].ivy;
    int kh = aCube[ivox].ivz;
    CVector3 cnt =  org + elen*CVector3(ih+0.5,jh+0.5,kh+0.5);
    {
      CVector3 q = nearest_Line_Point(cnt, src_pic, dir_pic);
      if( (q-cnt).Length() > elen  ) continue;
    }
    CVector3 aP[8] = {
      org + elen*CVector3(ih+0,jh+0,kh+0),
      org + elen*CVector3(ih+1,jh+0,kh+0),
      org + elen*CVector3(ih+0,jh+1,kh+0),
      org + elen*CVector3(ih+1,jh+1,kh+0),
      org + elen*CVector3(ih+0,jh+0,kh+1),
      org + elen*CVector3(ih+1,jh+0,kh+1),
      org + elen*CVector3(ih+0,jh+1,kh+1),
      org + elen*CVector3(ih+1,jh+1,kh+1) };
    for(int iface=0;iface<6;++iface){
      const CVector3& n = normalHexFace[iface];
      const CVector3& p0 = aP[noelElemFace_Vox[iface][0]];
      const CVector3& p1 = aP[noelElemFace_Vox[iface][1]];
      //      const CVector3& p2 = aP[noelHexFace[iface][2]];
      const CVector3& p3 = aP[noelElemFace_Vox[iface][3]];
      const CVector3 pi = intersection_Plane_Line(p0,n, src_pic,dir_pic);
      const double r0 = (pi-p0)*(p1-p0)/(elen*elen);
      const double r1 = (pi-p0)*(p3-p0)/(elen*elen);
      if( r0>0 && r0<1 && r1>0 && r1< 1 ){
        double depth = -(pi-src_pic)*dir_pic/dir_pic.DLength();
        if( icube_pic == -1 || depth < depth_min ){
          depth_min = depth;
          icube_pic = ivox;
          iface_pic = iface;
        }
      }
    }
  }
}

void Adj_CubeGrid
(int& ivx, int& ivy, int& ivz,
 int icube, int iface,
 std::vector<CCubeGrid>& aCube)
{
  ivx = aCube[icube].ivx;
  ivy = aCube[icube].ivy;
  ivz = aCube[icube].ivz;
  if( iface == 0 ){ ivx -= 1; }
  if( iface == 1 ){ ivx += 1; }
  if( iface == 2 ){ ivy -= 1; }
  if( iface == 3 ){ ivy += 1; }
  if( iface == 4 ){ ivz -= 1; }
  if( iface == 5 ){ ivz += 1; }
}

void Add_CubeGrid
(std::vector<CCubeGrid>& aCube,
 int ivx1, int ivy1, int ivz1)
{
  for(int ic=0;ic<aCube.size();++ic){
    if( aCube[ic].ivx == ivx1 && aCube[ic].ivy == ivy1 && aCube[ic].ivz == ivz1 ){
      if( aCube[ic].is_active ){
        return;
      }
      else{
        aCube[ic].is_active = true;
        return;
      }
    }
  }
  {
    CCubeGrid v(ivx1,ivy1,ivz1);
    aCube.push_back(v);
  }
}

void Del_CubeGrid
(std::vector<CCubeGrid>& aCube,
 int i1, int j1, int k1)
{
  for(int ic=0;ic<aCube.size();++ic){
    if( aCube[ic].ivx == i1 && aCube[ic].ivy == j1 && aCube[ic].ivz == k1 ){
      if( aCube[ic].is_active ){
        aCube[ic].is_active = false;
        return;
      }
      else{
        return;
      }
    }
  }
}


void AABB_CubeGrid
(int aabb[6],
 const std::vector<CCubeGrid>& aCube)
{
  if( aCube.size() == 0 ){
    aabb[0] = +1; aabb[1] = -1;
    aabb[2] = +1; aabb[3] = -1;
    aabb[4] = +1; aabb[5] = -1;
    return;
  }
  aabb[0] = aCube[0].ivx;  aabb[1] = aabb[0]+1;
  aabb[2] = aCube[0].ivy;  aabb[3] = aabb[0]+1;
  aabb[4] = aCube[0].ivz;  aabb[5] = aabb[0]+1;

  for(unsigned int ic=1;ic<aCube.size();++ic){
    if( aCube[ic].ivx+0 < aabb[0] ){ aabb[0] = aCube[ic].ivx+0; }
    if( aCube[ic].ivx+1 > aabb[1] ){ aabb[1] = aCube[ic].ivx+1; }
    if( aCube[ic].ivy+0 < aabb[2] ){ aabb[2] = aCube[ic].ivy+0; }
    if( aCube[ic].ivy+1 > aabb[3] ){ aabb[3] = aCube[ic].ivy+1; }
    if( aCube[ic].ivz+0 < aabb[4] ){ aabb[4] = aCube[ic].ivz+0; }
    if( aCube[ic].ivz+1 > aabb[5] ){ aabb[5] = aCube[ic].ivz+1; }
  }
  
}
