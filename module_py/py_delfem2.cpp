#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vector>
#include <map>
#include <deque>

#include "delfem2/funcs_glew.h" // have to be included in the beginning
#include "delfem2/mshtopoio_gl.h"

#include "delfem2/voxel.h"
#include "delfem2/bv.h"    // include gl
#include "delfem2/cad2d.h"

namespace py = pybind11;

// TODO:Make a wrapper class of the VoxelGrid?
CMeshElem MeshQuad3D_VoxelGrid(const CVoxelGrid& vg){
  CMeshElem me;
  vg.GetQuad(me.aPos, me.aElem);
  me.elem_type = MESHELEM_QUAD;
  me.ndim = 3;
  return me;
}

class CAxisXYZ {
public:
  CAxisXYZ(): len(1.0){
    line_width = 1.0;
  }
  CAxisXYZ(double len): len(len){
    line_width=1.0;
  }
  void Draw() const{
    glLineWidth(line_width);
    DrawAxis(len);
  }
  std::vector<double> MinMaxXYZ() const{
    std::vector<double> mm(6,0);
    mm[1] = len;  mm[3] = len;  mm[5] = len;
    return mm;
  }
public:
  double len;
  double line_width;
};

int add_edge
(int ip0, int ip1,
 double v0, double v1,
 std::map<std::pair<int,int>, std::pair<int,double> >& aNode)
{
  assert(v0*v1<0);
  auto key = std::make_pair(ip0,ip1);
  if( ip0 > ip1 ){
    key = std::make_pair(ip1,ip0);
  }
  auto res = aNode.find(key);
  if( res == aNode.end() ){
    int ino0 = aNode.size();
    auto val = std::make_pair(ino0,-v1/(v0-v1));
    aNode.insert(std::make_pair(key,val));
    return ino0;
  }
  else{
    return res->second.first;
  }
  return -1;
}

class CLevSetNode{
public:
  int ip0;
  int ip1;
  double ratio;
  int in0, in1;
public:
  CLevSetNode(){
    in0 = -1;
    in1 = -1;
  }
  void setNeighbour(int ino_in){
    if( in0 == - 1 ){ in0 = ino_in; }
    else{
      assert( in1 == -1 );
      in1 = ino_in;
    }
  }
  void pos(double& x0, double& y0, int nw){
    assert(ip0<ip1);
    int ih0 = ip0/nw;
    int iw0 = ip0-ih0*nw;
    int ih1 = ip1/nw;
    int iw1 = ip1-ih1*nw;
    x0 = ratio * iw0 + (1-ratio)*iw1;
    y0 = ratio * ih0 + (1-ratio)*ih1;
  }
};


std::vector<std::vector<int> > getLevelSet(const py::array_t<float>& a, double mag)
{
  assert(a.ndim()==2);
  const int mh = a.shape()[0];
  const int mw = a.shape()[1];
  const int nh = mh-1;
  const int nw = mw-1;
  const float* pVal = a.data();
  double thre = 0.5;
  std::map<std::pair<int,int>, std::pair<int,double> > mapNode;
  std::vector< std::pair<int,int> > aEdge;
  for(int ih=0;ih<nh;++ih){ // x
  for(int iw=0;iw<nw;++iw){ // y
    int ip0 = (ih+0)*mw+(iw+0);
    int ip1 = (ih+0)*mw+(iw+1);
    int ip2 = (ih+1)*mw+(iw+0);
    int ip3 = (ih+1)*mw+(iw+1);
    float v0 = pVal[ip0]-thre;
    float v1 = pVal[ip1]-thre;
    float v2 = pVal[ip2]-thre;
    float v3 = pVal[ip3]-thre;
    if( (v0>0 && v1<0 && v2<0 && v3<0) ||  (v0<0 && v1>0 && v2>0 && v3>0) ){
      int in0 = add_edge(ip0,ip1,v0,v1,mapNode);
      int in1 = add_edge(ip0,ip2,v0,v2,mapNode);
      aEdge.push_back(std::make_pair(in0,in1));
    }
    if( (v0<0 && v1>0 && v2<0 && v3<0) ||  (v0>0 && v1<0 && v2>0 && v3>0) ){
      int in0 = add_edge(ip0,ip1,v0,v1,mapNode);
      int in1 = add_edge(ip1,ip3,v1,v3,mapNode);
      aEdge.push_back(std::make_pair(in0,in1));
    }
    if( (v0<0 && v1<0 && v2>0 && v3<0) ||  (v0>0 && v1>0 && v2<0 && v3>0) ){
      int in0 = add_edge(ip0,ip2,v0,v2,mapNode);
      int in1 = add_edge(ip2,ip3,v2,v3,mapNode);
      aEdge.push_back(std::make_pair(in0,in1));
    }
    if( (v0<0 && v1<0 && v2<0 && v3>0) ||  (v0>0 && v1>0 && v2>0 && v3<0) ){
      int in0 = add_edge(ip1,ip3,v1,v3,mapNode);
      int in1 = add_edge(ip2,ip3,v2,v3,mapNode);
      aEdge.push_back(std::make_pair(in0,in1));
    }
    if( (v0>0 && v1>0 && v2<0 && v3<0) ||  (v0<0 && v1<0 && v2>0 && v3>0) ){
      int in0 = add_edge(ip0,ip2,v0,v2,mapNode);
      int in1 = add_edge(ip1,ip3,v1,v3,mapNode);
      aEdge.push_back(std::make_pair(in0,in1));
    }
    if( (v0>0 && v1<0 && v2>0 && v3<0) ||  (v0<0 && v1>0 && v2<0 && v3>0) ){
      int in0 = add_edge(ip0,ip1,v0,v1,mapNode);
      int in1 = add_edge(ip2,ip3,v2,v3,mapNode);
      aEdge.push_back(std::make_pair(in0,in1));
    }
    if( (v0>0 && v1<0 && v2<0 && v3>0) ||  (v0<0 && v1>0 && v2>0 && v3<0) ){
      int in0 = add_edge(ip0,ip1,v0,v1,mapNode);
      int in1 = add_edge(ip0,ip2,v0,v2,mapNode);
      aEdge.push_back(std::make_pair(in0,in1));
      int in2 = add_edge(ip1,ip3,v1,v3,mapNode);
      int in3 = add_edge(ip2,ip3,v2,v3,mapNode);
      aEdge.push_back(std::make_pair(in2,in3));
    }
  }
  }
  const int nNode = mapNode.size();
  std::vector<CLevSetNode> aNode(nNode);
  for(auto itr=mapNode.begin();itr!=mapNode.end();++itr){
    CLevSetNode node;
    node.ip0 = itr->first.first;
    node.ip1 = itr->first.second;
    node.ratio = itr->second.second;
    int in0 = itr->second.first;
    aNode[in0] = node;
  }
  for(unsigned int ie=0;ie<aEdge.size();++ie){
    int in0 = aEdge[ie].first;
    int in1 = aEdge[ie].second;
    aNode[in0].setNeighbour(in1);
    aNode[in1].setNeighbour(in0);
  }
  ///////////////////////////////////
  std::vector< std::deque<int> > aLoop;
  unsigned int ino_ker = 0;
  for(;;){
    for(;ino_ker<aNode.size();++ino_ker){
      bool is_included = false;
      for(unsigned int iloop=0;iloop<aLoop.size();++iloop){
        for(unsigned int jino=0;jino<aLoop[iloop].size();++jino){
          const int jno0 = aLoop[iloop][jino];
          if( jno0 == (int)ino_ker ){ is_included = true; break; }
        }
        if( is_included ){ break; }
      }
      if( !is_included ){ break; }
    }
    if( ino_ker == aNode.size() ){ break; }
//    std::cout << "inoker: " << ino_ker << " " << aNode.size() << std::endl;
    const int iloop0 = aLoop.size();
    aLoop.resize(aLoop.size()+1);
    aLoop[iloop0].push_back(ino_ker);
    if( aNode[ino_ker].in0 != -1 ){
      aLoop[iloop0].push_front(aNode[ino_ker].in0);
      for(;;){
        int in1 = aLoop[iloop0][0];
        int in2 = aLoop[iloop0][1];
        int jn0a = aNode[in1].in0;
        int jn0b = aNode[in1].in1;
        int in0 = -1;
        if( jn0a == in2 ){ in0 = jn0b; }
        if( jn0b == in2 ){ in0 = jn0a; }
        if( in0 == -1 ){ break; }
        if( in0 == aLoop[iloop0].back() ){ break; }
        aLoop[iloop0].push_front(in0);
      }
    }
    if( aLoop[iloop0].front() != aNode[ino_ker].in1 && aNode[ino_ker].in1 != -1 ){
      aLoop[iloop0].push_back(aNode[ino_ker].in1);
      for(;;){
        int in1 = aLoop[iloop0][aLoop[iloop0].size()-1];
        int in0 = aLoop[iloop0][aLoop[iloop0].size()-2];
        int jn2a = aNode[in1].in0;
        int jn2b = aNode[in1].in1;
        int in2 = -1;
        if( jn2a == in0 ){ in2 = jn2b; }
        if( jn2b == in0 ){ in2 = jn2a; }
        if( in2 == -1 ){ break; }
        if( in2 == aLoop[iloop0].front() ){ break; }
        aLoop[iloop0].push_back(in2);
      }
    }
//    break;
  }
  /////////
  std::vector< std::vector<int> > aRes;
  aRes.resize(aLoop.size());
  for(unsigned int iloop=0;iloop<aLoop.size();++iloop){
    for(unsigned int iin=0;iin<aLoop[iloop].size();++iin){
      int in0 = aLoop[iloop][iin];
      assert( in0>=0 && in0<(int)aNode.size() );
      double x0, y0; aNode[in0].pos(x0,y0,mw);
      aRes[iloop].push_back(int(x0*mag));
      aRes[iloop].push_back(int(y0*mag));
    }
  }
  return aRes;
}



//////////////////////////////////////////////////////////////////////////////////////////



void init_mshtopoio_gl(py::module &m);
void init_sampler(py::module &m);
void init_fbx(py::module &m);


PYBIND11_MODULE(dfm2, m) {
  m.doc() = "pybind11 delfem2 binding";
  ///////////////////////////////////
  
#ifdef USE_FBX
  init_fbx(m);
#endif

  ///////////////////////////////////
  // mesh
  init_mshtopoio_gl(m);
  
  ///////////////////////////////////
  init_sampler(m);
  
  ///////////////////////////////////
  // axis arrigned boudning box
  py::class_<CBV3D_AABB>(m,"AABB3", "3D axis aligned bounding box class")
  .def(py::init<>())
  .def(py::init<const std::vector<double>&>())
  .def("__str__",            &CBV3D_AABB::str, "print x_min,x_max,y_min,y_max,z_min,z_max")
  .def("minmax_xyz",         &CBV3D_AABB::MinMaxXYZ)
  .def("draw",               &CBV3D_AABB::Draw, "draw edge of the bounding box to opengl")
  .def("set_minmax_xyz",     &CBV3D_AABB::SetMinMaxXYZ)
  .def("add_minmax_xyz",     &CBV3D_AABB::Add_AABBMinMax)
  .def("list_xyz",           &CBV3D_AABB::Point3D_Vox, "corner xyz coords in voxel point order")
  .def("diagonal_length",    &CBV3D_AABB::DiagonalLength, "diagonal length of the bounding box")
  .def("max_length",         &CBV3D_AABB::MaxLength, "diagonal length of the bounding box")
  .def("center",             &CBV3D_AABB::Center, "center position")
  .def_readwrite("isActive", &CBV3D_AABB::is_active);
  
  py::class_<CAxisXYZ>(m,"AxisXYZ","3D axis class")
  .def(py::init<>())
  .def(py::init<double>(), py::arg("len"))
  .def("draw",                 &CAxisXYZ::Draw)
  .def("minmax_xyz",           &CAxisXYZ::MinMaxXYZ)
  .def_readwrite("len",        &CAxisXYZ::len)
  .def_readwrite("line_width", &CAxisXYZ::line_width);
  
  ///////////////////////////////////
  // voxel
  py::class_<CVoxelGrid>(m, "VoxelGrid", "voxel grid class")
  .def(py::init<>())
  .def("add",&CVoxelGrid::Add,"add voxel at the integer coordinate");
  
  m.def("meshQuad3d_voxelGrid",  &MeshQuad3D_VoxelGrid, "get quad mesh from voxel grid");
  
  ///////////////////////////////////
  // cad
  py::class_<CCad2D>(m, "Cad2D", "2D CAD class")
  .def(py::init<>())
  .def("add_square", &CCad2D::Add_Square)
  .def("draw",       &CCad2D::Draw)
  .def("minmax_xyz", &CCad2D::MinMaxXYZ);

  ///////////////////////////////////
  // gl misc
  m.def("setSomeLighting",  &setSomeLighting, "set some lighting that looks good for me");
  m.def("setup_glsl", setUpGLSL, "compile shader program");
  m.def("glew_init", glewInit);
  
  m.def("get_level_set", getLevelSet);
}
