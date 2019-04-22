#ifndef DYNTRI_V2_H
#define DYNTRI_V2_H

#include <map>
#include <algorithm>
#include <stack>

#include "delfem2/vec3.h"
#include "delfem2/dyntri.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

class CMeshDensity
{
public:
  virtual double edgeLengthRatio(double px, double py) const = 0;
};

int delaunay_triangulation2(std::vector<int>& aTri_out,		// out
                            std::vector<double>& aXY_out, // out
                            std::vector<int>& aPtrVtxInd, // out
                            std::vector<int>& aVtxInd, // out
                            ////
                            bool is_add_point_boundary,
                            const std::vector<int>& aIndXYs, // in
                            const std::vector<double>& aXY_in, // ind
                            const double max_edge_length, // ind
                            const CMeshDensity& mesh_density); // ind

// TODO: there should be three optoions for adding point on edge, 0:none, 1:only between input points, 2: resample everything
bool GenerateTesselation2(std::vector<int>& aTri_out, // out
                          std::vector<double>& aXY_out, // out
                          std::vector<int>& aPtrVtxInd,
                          std::vector<int>& aVtxInd,
                          ////
                          double elen,
                          const CMeshDensity& mesh_density,
                          bool is_uniform_resample_loop, // good for polyline curve in
                          const std::vector< std::vector<double> >& aVecAry0); // in



bool GenerateTesselation2(std::vector<int>& aTri_out, // out
                          std::vector<double>& aXY_out, // out
                          std::vector<int>& aPtrVtxInd,
                          std::vector<int>& aVtxInd,
                          ////
                          double elen,
                          bool is_uniform_resample_loop, // good for polyline curve
                          const std::vector< std::vector<double> >& aVecAry0); // in

#endif // #endif SURFACE_MESH_H
