#ifndef FUNCS_GL_H
#define FUNCS_GL_H

#include <string>
#include <vector>

// x = ax*[x] + bx
// y = ay*[y] + by
void DrawCharacter(int* pChr,
                   double ax, double bx,
                   double ay, double by);

// x = ax*[x] + bx
// y = ay*[y] + by
void DrawCharacter(char ic,
                   double ax, double bx,
                   double ay, double by);


void setSomeLighting();
void setSomeLighting2();
void setSomeLighting3();
void drawFloorShadow(void (*DrawObject)(), float yfloor, float wfloor);


void DrawRectangle_FullCanvas();
void showdepth();

//////////////////////////////////////////////////////////////////////////////
// texture related funcitons

void LoadImage_PPM(const std::string& filename,
                   std::vector<unsigned char>& image,
                   int& width, int& height);

void SaveImage(const std::string& path);

class SFile_TGA
{
public:
  unsigned char imageTypeCode;
  short int imageWidth;
  short int imageHeight;
  unsigned char bitCount;
  unsigned char *imageData;
};

bool LoadTGAFile(const char *filename, SFile_TGA *tgaFile);
int ReadPPM_SetTexture(const std::string& fname);

unsigned int LoadTexture(const unsigned char* image,
                   const int width, const int height, const int bpp);
void DrawTextureBackground(const unsigned int tex,
                           const int imgWidth,
                           const int imgHeight,
                           const int winWidth,
                           const int winHeight);

class CTextureInfo
{
public:
  std::string full_path;
  int width, height, bpp; // byte par pixel
  int id_tex_gl;
};

class CTexManager
{
public:
  void Clear();
  void AddTexture(const unsigned char* pixels,
                  const std::string& path,
                  int width, int height, int bpp)
  {
    const int id_tex_gl = LoadTexture(pixels, width,height,bpp);
    CTextureInfo texinfo;
    texinfo.full_path = path;
    texinfo.height = height;
    texinfo.width = width;
    texinfo.bpp = bpp;
    texinfo.id_tex_gl = id_tex_gl;
    /////
    bool is_new = true;
    for(int itex=0;itex<(int)aTexInfo.size();++itex){
      if( aTexInfo[itex].full_path != path ) continue;
      aTexInfo[itex] = texinfo;
      is_new = false;
    }
    if( is_new ){
      aTexInfo.push_back(texinfo);
    }
  }
  void AddPath(const std::string& path){
    CTextureInfo tex;
    tex.width = -1;
    tex.height = -1;
    tex.bpp = -1;
    tex.id_tex_gl = -1;
    tex.full_path = path;
    aTexInfo.push_back(tex);
  }
  void BindTexturePath(const std::string& path) const;
public:
  std::vector<CTextureInfo> aTexInfo;
};

////////////////////////////////////////////////////////////////////////
// draw functions

void DrawAxis(double s);
void DrawSphere(int nla, int nlo);
void DrawSphereAt(int nla, int nlo, double rad, double x, double y, double z);
void DrawBox_MinMaxXYZ(double x_min, double x_max,
                       double y_min, double y_max,
                       double z_min, double z_max);
void DrawBox_MinMaxXYZ(double aabbMinMaxXYZ[6]);
void DrawAABB3D_Edge(double cx, double cy, double cz,
                     double wx, double wy, double wz);
void DrawAABB3D_Edge(const double cw[6]);
void Draw_AABB3D_MinMaxXYZ_Edge(double x_min, double x_max,
                                double y_min, double y_max,
                                double z_min, double z_max);

///////////////
// Draw Point

void DrawPoints2D_Vectors(std::vector<double>& aXY,
                          std::vector<double>& aVal,
                          int nstride,
                          int noffset,
                          double mag);
void DrawPoints2D_Points(std::vector<double>& aXY);
void DrawPoints3D_Points(std::vector<double>& aXYZ);

///////////////
// Draw Tri


void DrawMeshTri2D_Face(const std::vector<int>& aTri,
                        const std::vector<double>& aXY);
void DrawMeshTri2D_FaceDisp2D(const std::vector<int>& aTri,
                              const std::vector<double>& aXY,
                              const std::vector<double>& aDisp);
void DrawMeshTri2D_Edge(const std::vector<int>& aTri,
                        const std::vector<double>& aXY);

void DrawMeshTri3D_FaceEdge(const std::vector<double>& aXYZ,
                        const std::vector<int>& aTri);
void DrawMeshTri3D_FaceNorm(const std::vector<double>& aXYZ,
                        const std::vector<int>& aTri,
                        const std::vector<double>& aNorm);
void DrawMeshTri3D_FaceNorm(const std::vector<double>& aXYZ,
                        const std::vector<int>& aTriVtx,
                        const std::vector<double>& aNorm,
                        const std::vector<int>& aTriNrm);
void DrawMeshTri3D_FaceNorm(const std::vector<double>& aXYZ,
                        const std::vector<int>& aTri);
void DrawMeshTri3DPart_FaceNorm(const std::vector<double>& aXYZ,
                            const std::vector<int>& aTri,
                            const std::vector<int>& aIndTri);
void DrawMeshTri3D_FaceNorm_Flg(const std::vector<double>& aXYZ,
                            const std::vector<int>& aTri,
                            int iflg,
                            const std::vector<int>& aFlgTri);
void DrawMeshTri3D_FaceNorm_XYsym(const std::vector<double>& aXYZ,
                              const std::vector<int>& aTri);
void DrawMeshTri3D_FaceNormEdge(const std::vector<double>& aXYZ,
                            const std::vector<int>& aTri);
void DrawMeshTri3D_FaceNorm_TexFace(const std::vector<double>& aXYZ,
                                    const std::vector<int>& aTri,
                                    const std::vector<double>& aTex);
void DrawMeshTri3D_FaceNorm_TexVtx(const std::vector<double>& aXYZ,
                                   const std::vector<int>& aTri,
                                   const std::vector<double>& aTex);
void DrawMeshTri3D_Edge(const std::vector<double>& aXYZ,
                    const std::vector<int>& aTri);
void DrawMeshTriMap3D_Edge(const std::vector<double>& aXYZ,
                       const std::vector<int>& aTri,
                       const std::vector<int>& map);

///////////////
// Draw Quad

void DrawMeshQuad3D_Edge(const std::vector<double>& aXYZ,
                     const std::vector<int>& aQuad);
void DrawMeshQuad3D_FaceNorm(const std::vector<double>& aXYZ,
                             const std::vector<int>& aQuad);
void DrawMeshQuad2D_Edge(const std::vector<double>& aXY,
                         const std::vector<int>& aQuad);

///////////////
// Draw Tet

void DrawMeshTet3D_Edge(const std::vector<double>& aXYZ,
                    const std::vector<int>& aTet);
void DrawMeshTet3D_EdgeDisp(const std::vector<double>& aXYZ,
                        const std::vector<int>& aTet,
                        const std::vector<double>& aDisp);
void DrawMeshTet3D_FaceNormal(const std::vector<double>& aXYZ,
                          const std::vector<int>& aTet);
void DrawMeshTet3DSurface_FaceNorm(const std::vector<double>& aXYZ,
                               const std::vector<int>& aTet,
                               const std::vector<int>& aTetFace);
void DrawMeshTet3DSurface_Edge(const std::vector<double>& aXYZ,
                           const std::vector<int>& aTet,
                           const std::vector<int>& aTetFace);
void DrawMeshTet3D_FaceNormDisp(const std::vector<double>& aXYZ,
                            const std::vector<int>& aTet,
                            const std::vector<double>& aDisp);

///////////////
// Draw Hex

void DrawHex3D_Edge(const std::vector<double>& aXYZ,
                    const std::vector<int>& aHex);
void DrawHex3D_FaceNorm(const std::vector<double>& aXYZ,
                        const std::vector<int>& aHex);
void Draw_HexMeshFaceDisp(const std::vector<double>& aXYZ,
                          const std::vector<int>& aHex,
                          const std::vector<double>& aDisp);


////////////////
// Draw Mix

void DrawMeshElem3D_FaceNorm(const std::vector<double>& aXYZ,
                             const std::vector<int>& aElemInd,
                             const std::vector<int>& aElem);

void DrawMeshElem3D_FaceNorm(const std::vector<double>& aXYZ,
                             const std::vector<int>& aElemInd,
                             const std::vector<int>& aElem,
                             const std::vector<double>& aUV);
void DrawMeshElemPart3D_FaceNorm_TexPoEl(const std::vector<double>& aXYZ,
                                         const std::vector<int>& aElemInd,
                                         const std::vector<int>& aElem,
                                         const std::vector<int>& aIndElem,
                                         const std::vector<double>& aUV);

#endif /* utility_gl_h */
