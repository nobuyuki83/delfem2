#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

#include "delfem2/funcs_gl.h"

#include "delfem2/mshio.h"
#include "delfem2/mat3.h"
#include "delfem2/rigmesh.h"
#include "delfem2/bv.h"    // include gl

#include "delfem2/../../external/io_fbx.h"

#define STB_IMAGE_IMPLEMENTATION
#include "delfem2/../../external/stb/stb_image.h"

CTexManager GetTexManager(const std::vector<std::string>& aPath){
  CTexManager tm;
  for(unsigned int ipath=0;ipath<aPath.size();++ipath){
    const std::string path = aPath[ipath];
    int width, height, bpp;
    unsigned char* pixels; pixels = stbi_load(path.c_str(), &width, &height, &bpp, 0);
    if( pixels == 0){ continue; }
    stbi__vertical_flip(pixels, width, height, bpp);
    tm.AddTexture(pixels, path, width, height, bpp);
    stbi_image_free(pixels);
  }
  return tm;
}

class CRigMshTex{
public:
  CRigMshTex(){
    is_draw_skeleton = true;
    color_bone_weight.assign(4, 0.0);
    color_bone_weight[0] = 1;
    color_bone_weight[3] = 1;
    color_bone_weight_back.assign(4, 1.0);
  }
  CRigMshTex(const std::string& fpath){
    is_draw_skeleton = true;
    color_bone_weight.assign(4, 0.0);
    color_bone_weight[0] = 1;
    color_bone_weight[3] = 1;
    color_bone_weight_back.assign(4, 1.0);
    this->Read(fpath);
  }
  void Draw(){
    rigmsh.is_draw_bone = is_draw_skeleton;
    rigmsh.color_bone_weight_back = color_bone_weight_back;
    rigmsh.color_bone_weight = color_bone_weight;
    rigmsh.Draw(tm);
  }
  std::vector<double> MinMaxXYZ(){
    return rigmsh.MinMaxXYZ();
  }
  void Scale(double scale){
    double A[16]; SetAffine_Scale(A,scale);
    rigmsh.Affine(A);
  }
  void Translate(double dx, double dy, double dz){
    double A[16]; SetAffine_Trans(A,
                                  dx,dy,dz);
    rigmsh.Affine(A);
  }
  void Rotate(double dx, double dy, double dz){
    double A[16]; SetAffine_Rotate_Rodriguez(A,
                                             dx,dy,dz);
    rigmsh.Affine(A);
  }
  void Read(const std::string& fpath){
    Read_FBX(fpath,rigmsh);
    rigmsh.FixTexPath(fpath);
    ////
    tm.Clear();
    std::vector< std::string > aTexPath = rigmsh.GetArrayTexPath();
    for(unsigned int ipath=0;ipath<aTexPath.size();++ipath){
      const std::string path = aTexPath[ipath];
      tm.AddPath(path);
    }
  }
  void LoadTex(){
    for(unsigned int it=0;it<tm.aTexInfo.size();++it){
      const std::string path = tm.aTexInfo[it].full_path;
      int width, height, bpp;
      unsigned char* pixels; pixels = stbi_load(path.c_str(), &width, &height, &bpp, 0);
      if( pixels == 0){ continue; }
      stbi__vertical_flip(pixels, width, height, bpp);
      tm.AddTexture(pixels, path, width, height, bpp);
      stbi_image_free(pixels);
    }
  }
  void PrintInfo(){
    rigmsh.PrintInfo();
  }
  CBV3D_AABB AABB3D(){
    std::vector<double> minmax_xyz = this->MinMaxXYZ();
    return CBV3D_AABB(minmax_xyz);
  }
  void SetBone_DrawBoneWeightOnMesh(int ibone){
    rigmsh.DisplayBoneWeightOnMesh(ibone);
  }
public:
  CRigMsh rigmsh;
  CTexManager tm;
  bool is_draw_skeleton;
  std::vector<double> color_bone_weight_back;
  std::vector<double> color_bone_weight;
};

namespace py = pybind11;

void init_fbx(py::module &m){
  py::class_<CRigMshTex>(m, "RigMshTex")
  .def(py::init<>())
  .def(py::init<const std::string&>())
  .def("draw",          &CRigMshTex::Draw)
  .def("minmax_xyz",    &CRigMshTex::MinMaxXYZ)
  .def("init_gl",       &CRigMshTex::LoadTex)
  .def("open",          &CRigMshTex::Read)
  .def("scale",         &CRigMshTex::Scale)
  .def("translate",     &CRigMshTex::Translate)
  .def("rotate",        &CRigMshTex::Rotate)
  .def("aabb",          &CRigMshTex::AABB3D)
  .def("info",          &CRigMshTex::PrintInfo)
  .def("set_bone_draw_bone_weight", &CRigMshTex::SetBone_DrawBoneWeightOnMesh)
  .def_readwrite("is_draw_skeleton", &CRigMshTex::is_draw_skeleton)
  .def_readwrite("color_bone_weight", &CRigMshTex::color_bone_weight)
  .def_readwrite("color_bone_weight_back", &CRigMshTex::color_bone_weight_back);

  
  m.def("get_texture_manager",GetTexManager);
  py::class_<CTexManager>(m,"TexManager")
  .def(py::init<>());
}

