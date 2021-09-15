/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/cad2_io_svg.h"

#include <cstdio>
#include <deque>

#include "delfem2/geo_polyline2.h"
#include "delfem2/str.h"
#include "delfem2/file.h"

// ------------------------------------------

namespace delfem2 {
namespace cad2_io_svg {

DFM2_INLINE std::vector<std::string> SVG_Split_Path_d
    (std::string &s0) {
  unsigned int imark = 0;
  std::vector<std::string> aS;
  for (unsigned int i = 0; i < s0.size(); ++i) {
    if (isNumeric(s0[i])) continue;
    if (s0[i] == ',') {
      std::string s1(s0.begin() + imark, s0.begin() + i);
      aS.push_back(s1);
      imark = i + 1;  // mark shoud be the begining position of the string so move next
    }
    if (s0[i] == ' ') { // sometimes the space act as delimiter in the SVG (inkscape version)
      if (i > imark) {
        std::string s1(s0.begin() + imark, s0.begin() + i);
        aS.push_back(s1);
      }
      imark = i+1; // mark shoud be the begining position of the string so move next
    }
    if (s0[i] == '-') {
      if (i > imark) {
        std::string s1(s0.begin() + imark, s0.begin() + i);
        aS.push_back(s1);
      }
      imark = i;
    }
    if (isAlphabet(s0[i])) {
      if (i > imark) {
        std::string s1(s0.begin() + imark, s0.begin() + i);
        aS.push_back(s1);
      }
      const char s2[2] = {s0[i], '\0'};
      aS.emplace_back(s2);
      imark = i + 1;
    }
  }
  return aS;
}

DFM2_INLINE void LoopEdgeCad2D_SVGPathD(
    std::vector<CCad2D_EdgeGeo> &aEdge,
    std::vector<std::string> &aStr1)
{
  assert(aStr1[0] == "M" || aStr1[0] == "m" );
  assert(aStr1[aStr1.size() - 1] == "Z" || aStr1[aStr1.size() - 1] == "z");
  CVec2d pos_cur;
  for (int is = 0;;) {
    if (aStr1[is] == "M") {
      pos_cur.p[0] = myStod(aStr1[is + 1]);
      pos_cur.p[1] = myStod(aStr1[is + 2]);
      is += 3;
      for (;;) {
        if (isAlphabet(aStr1[is][0])) { break; }
        CCad2D_EdgeGeo e;
        e.p0 = pos_cur;
        e.p1.p[0] = myStod(aStr1[is + 0]);
        e.p1.p[1] = myStod(aStr1[is + 1]);
        pos_cur = e.p1;
        aEdge.push_back(e);
        is += 2;
      }
    } else if (aStr1[is] == "m") {
      pos_cur.p[0] = myStod(aStr1[is + 1]);
      pos_cur.p[1] = myStod(aStr1[is + 2]);
      is += 3;
      for (;;) {
        if (isAlphabet(aStr1[is][0])) { break; }
        CCad2D_EdgeGeo e;
        e.p0 = pos_cur;
        e.p1.p[0] = pos_cur.x + myStod(aStr1[is + 0]);
        e.p1.p[1] = pos_cur.y + myStod(aStr1[is + 1]);
        pos_cur = e.p1;
        aEdge.push_back(e);
        is += 2;
      }
    } else if (aStr1[is] == "C") { // cubic Bezier absolute coordinates
      ++is;
      for (;;) { // loop for poly Bezier
        CCad2D_EdgeGeo e;
        e.p0 = pos_cur;
        e.p1 = CVec2d(myStod(aStr1[is + 4]), myStod(aStr1[is + 5]));
        double len01 = (e.p1 - e.p0).norm();
        const CVec2d lx = (e.p1 - e.p0) / (len01 * len01);
        const CVec2d ly = CVec2d(lx.y, -lx.x);
        e.type_edge = CCad2D_EdgeGeo::BEZIER_CUBIC;
        e.param.resize(4, 0.0);
        CVec2d p2(myStod(aStr1[is + 0]), myStod(aStr1[is + 1]));
        CVec2d p3(myStod(aStr1[is + 2]), myStod(aStr1[is + 3]));
        e.param[0] = (p2 - e.p0).dot(lx);
        e.param[1] = (p2 - e.p0).dot(ly);
        e.param[2] = (p3 - e.p0).dot(lx);
        e.param[3] = (p3 - e.p0).dot(ly);
        aEdge.push_back(e);
        pos_cur = e.p1;
        is += 6;
        if (isAlphabet(aStr1[is][0])) { break; }
      }
    } else if (aStr1[is] == "c") {
      ++is; // 'c'
      for (;;) { // loop for poly-Bezeir curve
        CCad2D_EdgeGeo e;
        e.p0 = pos_cur;
        e.p1 = CVec2d(pos_cur.x + myStod(aStr1[is + 4]),
                      pos_cur.y + myStod(aStr1[is + 5]));
        const CVec2d p2(pos_cur.x+ myStod(aStr1[is + 0]),
                        pos_cur.y+ myStod(aStr1[is + 1]));
        const CVec2d p3(pos_cur.x + myStod(aStr1[is + 2]),
                        pos_cur.y + myStod(aStr1[is + 3]));
        const double len01 = (e.p1 - e.p0).norm();
        const CVec2d lx = (e.p1 - e.p0) / (len01 * len01);
        const CVec2d ly = CVec2d(lx.y, -lx.x);
        e.type_edge = CCad2D_EdgeGeo::BEZIER_CUBIC;
        e.param.resize(4, 0.0);
        e.param[0] = (p2 - e.p0).dot(lx);
        e.param[1] = (p2 - e.p0).dot(ly);
        e.param[2] = (p3 - e.p0).dot(lx);
        e.param[3] = (p3 - e.p0).dot(ly);
        aEdge.push_back(e);
        pos_cur = e.p1;
        is += 6;
        if (isAlphabet(aStr1[is][0])) { break; }
      }
    } else if (aStr1[is] == "l") {
      CCad2D_EdgeGeo e;
      e.p0 = pos_cur;
      e.p1.p[0] = pos_cur.x + myStod(aStr1[is + 1]);
      e.p1.p[1] = pos_cur.y + myStod(aStr1[is + 2]);
      pos_cur = e.p1;
      aEdge.push_back(e);
      is += 3;
      for (;;) {
        if (isAlphabet(aStr1[is][0])) { break; }
        CCad2D_EdgeGeo e0;
        e0.p0 = pos_cur;
        e0.p1.p[0] = pos_cur.x + myStod(aStr1[is + 0]);
        e0.p1.p[1] = pos_cur.y + myStod(aStr1[is + 1]);
        pos_cur = e0.p1;
        aEdge.push_back(e0);
        is += 2;
      }
    } else if (aStr1[is] == "L") {
      CCad2D_EdgeGeo e;
      e.p0 = pos_cur;
      e.p1.p[0] = myStod(aStr1[is + 1]);
      e.p1.p[1] = myStod(aStr1[is + 2]);
      pos_cur = e.p1;
      aEdge.push_back(e);
      is += 3;
      for (;;) {
        if (isAlphabet(aStr1[is][0])) { break; }
        CCad2D_EdgeGeo e0;
        e0.p0 = pos_cur;
        e0.p1.p[0] = myStod(aStr1[is + 0]);
        e0.p1.p[1] = myStod(aStr1[is + 1]);
        pos_cur = e0.p1;
        aEdge.push_back(e0);
        is += 2;
      }
    } else if (aStr1[is] == "S") {
      CCad2D_EdgeGeo e;
      e.p0 = pos_cur;
      e.p1 = CVec2d(myStod(aStr1[is + 3]), myStod(aStr1[is + 4]));
      double len01 = (e.p1 - e.p0).norm();
      const CVec2d lx = (e.p1 - e.p0) / (len01 * len01);
      const CVec2d ly = CVec2d(lx.y, -lx.x);
      e.type_edge = CCad2D_EdgeGeo::BEZIER_CUBIC;
      CVec2d p2(myStod(aStr1[is + 1]), myStod(aStr1[is + 2]));
      CVec2d p3 = e.p1;
      e.param.resize(4, 0.0);
      e.param[0] = (p2 - e.p0).dot(lx);
      e.param[1] = (p2 - e.p0).dot(ly);
      e.param[2] = (p3 - e.p0).dot(lx);
      e.param[3] = (p3 - e.p0).dot(ly);
      aEdge.push_back(e);
      pos_cur = e.p1;
      is += 5;
    } else if (aStr1[is] == "v") {
      CCad2D_EdgeGeo e;
      {
        e.p0 = pos_cur;
        e.p1 = CVec2d(pos_cur.x, pos_cur.y + myStod(aStr1[is + 1]));
        e.type_edge = CCad2D_EdgeGeo::LINE;
      }
      pos_cur = e.p1;
      aEdge.push_back(e);
      is += 2;
    } else if (aStr1[is] == "V") {
      CCad2D_EdgeGeo e;
      {
        e.p0 = pos_cur;
        e.p1 = CVec2d(pos_cur.x, myStod(aStr1[is + 1]));
        e.type_edge = CCad2D_EdgeGeo::LINE;
      }
      aEdge.push_back(e);
      pos_cur = e.p1;
      is += 2;
    } else if (aStr1[is] == "H") {
      CCad2D_EdgeGeo e;
      {
        e.p0 = pos_cur;
        e.p1 = CVec2d(myStod(aStr1[is + 1]), pos_cur.y);
        e.type_edge = CCad2D_EdgeGeo::LINE;
      }
      aEdge.push_back(e);
      pos_cur = e.p1;
      is += 2;
    } else if (aStr1[is] == "h") {
      CCad2D_EdgeGeo e;
      {
        e.p0 = pos_cur;
        e.p1 = CVec2d(e.p0.x + myStod(aStr1[is + 1]), pos_cur.y);
        e.type_edge = CCad2D_EdgeGeo::LINE;
      }
      aEdge.push_back(e);
      pos_cur = e.p1;
      is += 2;
    } else if (aStr1[is] == "q") {
      const CVec2d p0 = pos_cur;
      const CVec2d p1(p0.x + myStod(aStr1[is + 1]), p0.y + myStod(aStr1[is + 2]));
      const CVec2d p2(p0.x + myStod(aStr1[is + 3]), p0.y + myStod(aStr1[is + 4]));
      const CVec2d lx = (p2 - p0) / (p2 - p0).squaredNorm();
      const CVec2d ly = -CVec2d(lx.y, -lx.x);
      CCad2D_EdgeGeo e;
      {
        e.p0 = pos_cur;
        e.p1 = p2;
        e.type_edge = CCad2D_EdgeGeo::BEZIER_QUADRATIC;
        e.param.resize(2);
        e.param[0] = (p1 - p0).dot(lx);
        e.param[1] = (p1 - p0).dot(ly);
      }
      aEdge.push_back(e);
      pos_cur = e.p1;
      is += 5;
    } else if (aStr1[is] == "Q") {
      const CVec2d p0 = pos_cur;
      const CVec2d p1(myStod(aStr1[is + 1]), myStod(aStr1[is + 2]));
      const CVec2d p2(myStod(aStr1[is + 3]), myStod(aStr1[is + 4]));
      const CVec2d lx = (p2 - p0) / (p2 - p0).squaredNorm();
      const CVec2d ly = -CVec2d(lx.y, -lx.x);
      CCad2D_EdgeGeo e;
      {
        e.p0 = pos_cur;
        e.p1 = p2;
        e.type_edge = CCad2D_EdgeGeo::BEZIER_QUADRATIC;
        e.param.resize(2);
        e.param[0] = (p1 - p0).dot(lx);
        e.param[1] = (p1 - p0).dot(ly);
      }
      aEdge.push_back(e);
      pos_cur = e.p1;
      is += 5;
    } else if (aStr1[is] == "z" || aStr1[is] == "Z") {
      const CVec2d p1 = aEdge[0].p0;
      const CVec2d p0 = aEdge[aEdge.size() - 1].p1;
      double dist0 = (p0 - p1).norm();
      if (dist0 > 1.0e-9) {
        CCad2D_EdgeGeo e;
        e.p0 = p0;
        e.p1 = p1;
        e.type_edge = CCad2D_EdgeGeo::LINE;
        aEdge.push_back(e);
      }
      break;
    } else {
      std::cout << "error!--> " << aStr1[is] << std::endl;
      break;
    }
  }
}

DFM2_INLINE void LoopEdgeCad2D_SVGPolygonPoints
    (std::vector<CCad2D_EdgeGeo> &aEdge,
     std::vector<std::string> &aS) {
  const size_t np = aS.size() / 2;
  std::vector<CVec2d> aP;
  for (size_t ip = 0; ip < np; ++ip) {
    aP.emplace_back(myStod(aS[ip * 2 + 0]), myStod(aS[ip * 2 + 1]));
  }
  for (size_t ie = 0; ie < np; ++ie) {
    CCad2D_EdgeGeo e;
    e.p0 = aP[(ie + 0) % np];
    e.p1 = aP[(ie + 1) % np];
    aEdge.push_back(e);
  }
}

} // cad2
} // delfem2

// =========================================================

DFM2_INLINE void delfem2::ReadSVG_LoopEdgeCCad2D(
    std::vector< std::vector<CCad2D_EdgeGeo>> & aaEdge,
    const std::string& fname)
{
  namespace lcl = ::delfem2::cad2_io_svg;
  aaEdge.clear();
  std::vector<char> aC;
  if( !GetFileContents(aC, fname) ){ return; }
  
  // ----
  /*
  std::cout << "svg file content: ";
  for(unsigned int ic=0;ic<aC.size();++ic){ std::cout << aC[ic]; }
  std::cout << std::endl;
  */
  // ----
  
  std::vector< std::string > aStrTagContent;
  XML_SeparateTagContent(aStrTagContent,
                         aC);
  
  { // get path
    for(auto & sTagContent : aStrTagContent){
//      std::cout << "tagcontent: " << sTagContent << std::endl;
      std::string str_path;
      if( sTagContent.compare(0,4,"path") == 0 || // adobe illustrator
          sTagContent.compare(0,5,"path\r") == 0){ // inkscape
        str_path = std::string(sTagContent.begin()+5,sTagContent.end());
      }
      if( str_path == "" ){ continue; }
      // remove new line codes as inkscape svg has new line in side path tag
      // don't remove spaces here as inkscape svg uses space for the delimiter
      str_path = Remove(str_path, "\n\r");
//      std::cout << "str_path: " << str_path << std::endl;
      std::map< std::string, std::string > mapAttr;
      ParseAttributes(mapAttr,
                      str_path);
      std::string str_path_d = mapAttr["d"];
//      std::cout << "str_path_d: " << str_path_d << std::endl;
      if( str_path_d.empty() ){ continue; }
//      std::cout << "str_path_d: " << str_path_d << std::endl;
      std::vector<std::string> aStr1 = lcl::SVG_Split_Path_d(str_path_d);
      /*
      for(unsigned int is=0;is<aStr1.size();++is){
        std::cout << is << " " << aStr1[is] << std::endl;
      }
       */
      std::vector<CCad2D_EdgeGeo> aEdge;
      lcl::LoopEdgeCad2D_SVGPathD(aEdge,
          aStr1);
      /*
      for(int ie=0;ie<aEdge.size();++ie){
        std::cout << ie << " " << aEdge.size() << " " << aEdge[ie].param.size() << std::endl;
      }
       */
      aaEdge.push_back(aEdge);
    }
  }
  
  { // get polygon
    for(auto & sTagContent : aStrTagContent){
      std::string str_polygon;
      if( sTagContent.compare(0,8,"polygon ") == 0 ){
        str_polygon = std::string(sTagContent.begin()+8,sTagContent.end());
      }
//    std::cout << "str_polygon: " << str_polygon << std::endl;
      if( str_polygon == "" ){ continue; }
      std::map< std::string, std::string > mapAttr;
      ParseAttributes(mapAttr,
                      str_polygon);
      std::string str_polygon_points = mapAttr["points"];
      std::vector<std::string> aS = Split(str_polygon_points, "  ,");
      /*
      for(unsigned int is=0;is<aS.size();++is){
        std::cout << is << " " << aS[is] << std::endl;
      }
       */
      std::vector<CCad2D_EdgeGeo> aEdge;
      lcl::LoopEdgeCad2D_SVGPolygonPoints(
          aEdge,
          aS);
      aaEdge.push_back(aEdge);
    }
  }
}


DFM2_INLINE void delfem2::ReadSVG_Cad2D(
    delfem2::CCad2D& cad,
    const std::string& fpath,
    double scale)
{
  std::vector< std::vector<delfem2::CCad2D_EdgeGeo> > aaEdge;
  ReadSVG_LoopEdgeCCad2D(aaEdge,
                         fpath);
  cad.Clear();
  for(unsigned int iae=0;iae<aaEdge.size();++iae){
    std::vector<delfem2::CCad2D_EdgeGeo> aEdge = aaEdge[iae];
    Transform_LoopEdgeCad2D(aEdge,false,true,scale,scale);
    if( AreaLoop(aEdge) < 0 ){ aEdge = InvertLoop(aEdge); }
    aEdge = RemoveEdgeWithZeroLength(aEdge);
    for(auto & ie : aEdge){ ie.GenMeshLength(-1); }
    cad.AddFace(aEdge);
  }
}