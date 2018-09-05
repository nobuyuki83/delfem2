#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

#include "delfem2/funcs_gl.h"

#include "delfem2/mshio.h"

#include "delfem2/rigmesh.h"

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
  void Draw(){
    rigmsh.Draw(tm);
  }
  std::vector<double> MinMaxXYZ(){
    return rigmsh.MinMaxXYZ();
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
public:
  CRigMsh rigmsh;
  CTexManager tm;
};

namespace py = pybind11;

void init_fbx(py::module &m){
  py::class_<CRigMshTex>(m, "RigMshTex")
  .def(py::init<>())
  .def("draw",          &CRigMshTex::Draw)
  .def("minmax_xyz",    &CRigMshTex::MinMaxXYZ)
  .def("init_gl",       &CRigMshTex::LoadTex)
  .def("open",          &CRigMshTex::Read)
  .def("info",          &CRigMshTex::PrintInfo);

  m.def("get_texture_manager",GetTexManager);
  py::class_<CTexManager>(m,"TexManager")
  .def(py::init<>());
}

