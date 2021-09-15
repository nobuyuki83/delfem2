/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/cad2_io_dxf.h"

#include <cstdio>
#include <deque>
#include <climits>

#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/geo_polyline2.h"
#include "delfem2/str.h"
#include "delfem2/file.h"

// ------------------------------------------

DFM2_INLINE bool delfem2::WriteCAD_DXF(
    const std::string& file_name,
    const CCad2D& cad,
    double scale)
{
  FILE *fp;
  if( (fp = ::fopen(file_name.c_str(),"w"))== nullptr ){
    fclose(fp);
    assert(0);
    return false;
  }
  CBoundingBox2<double> bb;
  {  // Get Bounding Box of this Object
    for(const auto & edge : cad.aEdge){
      bb += edge.BB();
    }
  }
  // header section
  fprintf(fp, "  0\nSECTION\n");
  fprintf(fp, "  2\nHEADER\n");
  fprintf(fp, "  9\n$ACADVER\n  1\nAC1009\n");
  fprintf(fp, "  9\n$EXTMIN\n  10\n%lf\n  20\n%lf\n",bb.x_min*scale,bb.y_min*scale);
  fprintf(fp, "  9\n$EXTMAX\n  10\n%lf\n  20\n%lf\n",bb.x_max*scale,bb.y_max*scale);
  fprintf(fp, "  0\nENDSEC\n");
  // table section
  fprintf(fp, "  0\nSECTION\n");
  fprintf(fp, "  2\nTABLES\n");
  fprintf(fp, "  0\nENDSEC\n");
  // block section
  fprintf(fp, "  0\nSECTION\n");
  fprintf(fp, "  2\nBLOCKS\n");
  fprintf(fp, "  0\nENDSEC\n");
  // entity section
  fprintf(fp,"  0\nSECTION\n");
  fprintf(fp,"  2\nENTITIES\n");
  for(size_t ifc=0;ifc<cad.aFace.size();++ifc){
    const std::vector<int>& aIL = cad.topo.aFace[ifc].aIL;
    for(int il0 : aIL){
      const std::vector< std::pair<int,bool> >& aIE = cad.topo.aLoop[il0].aIE;
      for(const auto & iie : aIE){
        unsigned int ie0 = iie.first;
//        bool dir0 = aIE[iie].second;
        unsigned int id_vs = cad.topo.aEdge[ie0].iv0;
        unsigned int id_ve = cad.topo.aEdge[ie0].iv1;
        const CVec2d& ps = cad.aVtx[id_vs].pos;
        const CVec2d& pe = cad.aVtx[id_ve].pos;
        if( cad.aEdge[ie0].type_edge == 0 ){
          fprintf(fp,"  0\nLINE\n  8\n%d\n  6\nCONTINUOUS\n  62\n7\n",il0);
          fprintf(fp,"  10\n%lf\n",ps.x*scale);
          fprintf(fp,"  20\n%lf\n",ps.y*scale);
          fprintf(fp,"  11\n%lf\n",pe.x*scale);
          fprintf(fp,"  21\n%lf\n",pe.y*scale);
        }
        /*
        else if( this->GetEdgeCurveType(id_e) == 1 ){ // Arc
          const CEdge2D& edge = this->GetEdge(id_e);
          CVector2D pc;  double r;
          edge.GetCenterRadius(pc,r);
          double d1, d2;
          {
            CVector2D vs = ps - pc;
            CVector2D ve = pe - pc;
            double ds = atan2(vs.y,vs.x); ds = ds * 180.0 / 3.14159265; if( ds < 0.0 ) ds += 360;
            double de = atan2(ve.y,ve.x); de = de * 180.0 / 3.14159265; if( de < 0.0 ) de += 360;
            if( edge.is_left_side ){ d1 = de; d2 = ds; }
            else{                    d1 = ds; d2 = de; }
          }
          fprintf(fp,"  0\nARC\n  8\n%d\n  6\nCONTINUOUS\n  62\n7\n  100\nAcDbCircle\n",id_l);
          fprintf(fp,"  10\n%lf\n",pc.x*scale);  // x coord
          fprintf(fp,"  20\n%lf\n",pc.y*scale);  // y coord
          fprintf(fp,"  40\n%lf\n",r*scale);  // radius
          fprintf(fp,"  100\nAcDbArc\n");
          fprintf(fp,"  50\n%lf\n",d1);
          fprintf(fp,"  51\n%lf\n",d2);
        }
        else if( this->GetEdgeCurveType(id_e) == 2 ){ // polyline
          const CEdge2D& edge = this->GetEdge(id_e);
          fprintf(fp,"  0\nPOLYLINE\n  8\n%d\n  6\nCONTINUOUS\n",id_l);
          fprintf(fp,"  10\n0.0\n");
          fprintf(fp,"  20\n0.0\n");
          fprintf(fp,"  30\n0.0\n");
          fprintf(fp,"  70\n8\n");
          fprintf(fp,"  66\n1\n");
          ////
          const std::vector<double>& axys = edge.aRelCoMesh;
          assert( axys.size() % 2 == 0 );
          const unsigned int nno = axys.size()/2;
          const Com::CVector2D& po_s = this->GetVertexCoord( edge.id_v_s );
          const Com::CVector2D& po_e = this->GetVertexCoord( edge.id_v_e );
          Com::CVector2D v0 = po_e-po_s;
          Com::CVector2D v1(-v0.y,v0.x);
          fprintf(fp,"  0\nVERTEX\n 8\n0\n 10\n%lf\n 20\n%lf\n 30\n%lf\n", po_s.x*scale, po_s.y*scale, 0.0);
          for(unsigned int ino=0;ino<nno;ino++){
            const Com::CVector2D& p = po_s + v0*axys[ino*2+0] + v1*axys[ino*2+1];
            fprintf(fp,"  0\nVERTEX\n 8\n0\n 10\n%lf\n 20\n%lf\n 30\n%lf\n", p.x*scale, p.y*scale, 0.0);
          }
          fprintf(fp,"  0\nVERTEX\n 8\n0\n 10\n%lf\n 20\n%lf\n 30\n%lf\n", po_e.x*scale, po_e.y*scale, 0.0);
          fprintf(fp,"  0\nSEQEND\n");
        }
         */
      }
    }
  }
  fprintf(fp, "  0\nENDSEC\n  0\nEOF\n");
  fclose(fp);
  return true;
}
