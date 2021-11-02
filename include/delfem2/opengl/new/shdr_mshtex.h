#ifndef DFM2_OPENGL_NEW_SHDR_MSHTEX_H
#define DFM2_OPENGL_NEW_SHDR_MSHTEX_H

#include <cstdio>
#include <vector>
#include <set>
#include <climits>

#include "delfem2/opengl/new/funcs.h" // CGL4_VAO_Mesh
#include "delfem2/mshuni.h"
#include "delfem2/dfm2_inline.h"

// -------------------------------------

namespace delfem2::opengl {

class Drawer_MeshTex {
public:
  void SetElement(
      std::vector<unsigned int> &elem_vtx,
      int gl_primitive_type);

  template <typename REAL>
  void SetCoords(
      std::vector<REAL> &vtx_coords,
      unsigned int ndim);

  template <typename REAL>
  void SetTexUV(
      std::vector<REAL> &aTex);

  virtual void InitGL();

  void Draw(const float mat4_projection[16],
            const float mat4_modelview[16]) const;

public:
  VertexArrayObject vao; // gl4
  int shaderProgram = -1;
  int Loc_MatrixProjection = -1;
  int Loc_MatrixModelView = -1;
  int Loc_Texture = -1;
};

// ----------------------------------

class Drawer_RectangleTex : public Drawer_MeshTex {
public:
  Drawer_RectangleTex(float xmin, float xmax, float ymin, float ymax, float z=0.f)
  : xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), z(z){}
  
  Drawer_RectangleTex(float half_size)
  : xmin(-half_size), xmax(+half_size), ymin(-half_size), ymax(+half_size){}
  
  Drawer_RectangleTex() = default;
  
  // --------------------
  
  void InitGL() override {
    Drawer_MeshTex::InitGL();
    std::vector<float> aPos3d = {
      xmin, ymin, z,
      xmax, ymin, z,
      xmax, ymax, z,
      xmin, ymax, z,
    };
    std::vector<unsigned int> aTri = {
      0,1,2,
      0,2,3,
    };
    std::vector<float> aTex2d = {
      0.0, 0.0,
      1.0, 0.0,
      1.0, 1.0,
      0.0, 1.0
    };
    Drawer_MeshTex::SetCoords(aPos3d,3);
    Drawer_MeshTex::SetTexUV(aTex2d);
    Drawer_MeshTex::SetElement( aTri, GL_TRIANGLES);
  }
private:
  float xmin = -1;
  float xmax = +1;
  float ymin = -1;
  float ymax = +1;
  float z = 0.0;
};

class Drawer_MeshTexSeparateIndexing : public Drawer_MeshTex {
 public:
  void Clear(){
    tri_uni.clear();
    uni_xyz.clear();
    uni_tex.clear();
  }

  template <typename REAL>
  void SetMesh(
      const std::vector<REAL> &vtx_xyz,
      const std::vector<REAL> &vtx_tex,
      const std::vector<unsigned int> &elem_vtx_xyz,
      const std::vector<unsigned int> &elem_vtx_tex)
  {
    this->Clear();
    namespace dfm2 = delfem2;
    const unsigned int nxyz = vtx_xyz.size()/3;
    // const unsigned int ntex = vtx_tex.size()/2;
    assert( elem_vtx_xyz.size() / 3 == elem_vtx_tex.size() / 3 );
    std::cout << vtx_xyz.size() / 3 << " " << vtx_tex.size() /2 << std::endl;
    std::vector<unsigned int> elsup0_ind, elsup0;
    dfm2::JArray_ElSuP_MeshTri(
        elsup0_ind, elsup0,
        elem_vtx_xyz, nxyz);
    std::vector<unsigned int> elsup1_ind, elsup1;
    dfm2::JArray_ElSuP_MeshTri(
        elsup1_ind, elsup1, elem_vtx_tex, vtx_tex.size()/2);
    tri_uni.resize(elem_vtx_xyz.size(),UINT_MAX);
    for(unsigned int itri=0;itri<elem_vtx_xyz.size()/3;++itri){
      for(unsigned int ino=0;ino<3;++ino){
        if( tri_uni[itri*3+ino] != UINT_MAX ) { continue; }
        unsigned int ixyz = elem_vtx_xyz[itri*3+ino];
        unsigned int itex = elem_vtx_tex[itri*3+ino];
        const std::set s0(elsup0.data()+elsup0_ind[ixyz], elsup0.data()+elsup0_ind[ixyz+1]);
        const std::set s1(elsup1.data()+elsup1_ind[itex], elsup1.data()+elsup1_ind[itex+1]);
        std::vector<unsigned int> tris;
        set_intersection(
            s0.begin(),s0.end(),
            s1.begin(),s1.end(),
            std::back_inserter(tris));
        if( tris.empty() ){ continue; }
        unsigned int iuni = uni_tex.size();
        uni_xyz.push_back(ixyz);
        uni_tex.push_back(itex);
        for(unsigned jtri : tris){
          const unsigned int jno = this->FindIndexTri(elem_vtx_xyz.data()+jtri*3, ixyz);
          assert( elem_vtx_xyz[jtri*3+jno] == ixyz );
          assert( elem_vtx_tex[jtri*3+jno] == itex );
          assert( tri_uni[jtri*3+jno] == UINT_MAX );
          tri_uni[jtri*3+jno] = iuni;
        }
      }
    }
    this->SetElement(tri_uni, GL_TRIANGLES);
    const unsigned int nuni = uni_xyz.size();
    {
      std::vector<double> uni_xyz0(nuni*3);
      for(unsigned int iuni=0;iuni<nuni;++iuni){
        unsigned int ixyz = uni_xyz[iuni];
        uni_xyz0[iuni*3+0] = vtx_xyz[ixyz*3+0];
        uni_xyz0[iuni*3+1] = vtx_xyz[ixyz*3+1];
        uni_xyz0[iuni*3+2] = vtx_xyz[ixyz*3+2];
      }
      this->SetCoords(uni_xyz0, 3);
    }
    {
      std::vector<double> uni_tex0(nuni*2);
      for(unsigned int iuni=0;iuni<nuni;++iuni){
        unsigned int itex = uni_tex[iuni];
        uni_tex0[iuni*2+0] = vtx_tex[itex*2+0];
        uni_tex0[iuni*2+1] = vtx_tex[itex*2+1];
      }
      this->SetTexUV(uni_tex0);
    }
  }
 private:
  unsigned int FindIndexTri(
      const unsigned int *tri_vtx,
      unsigned int ixyz) {
    if( tri_vtx[0] == ixyz ){ return 0; }
    if( tri_vtx[1] == ixyz ){ return 1; }
    if( tri_vtx[2] == ixyz ){ return 2; }
    assert(0);
    return UINT_MAX;
  }
 public:
  std::vector<unsigned int> tri_uni;
  std::vector<unsigned int> uni_xyz, uni_tex;
};

class Drawer_MeshTexSeparateIndexingMixTwoTexture
    : public Drawer_MeshTexSeparateIndexing{
 public:
  void CompileShader();
  void Draw(
      const float mat4_projection[16],
      const float mat4_modelview[16],
      float ratio) const;
 public:
  unsigned int loc_ratio;
  unsigned int loc_tex0;
  unsigned int loc_tex1;
};


}  // delfem2::opengl


#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/new/shdr_mshtex.cpp"
#endif



#endif  // DFM2_OPENGL_NEW_SHDR_MSHTEX_H
