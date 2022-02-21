#include "openvdb/openvdb.h"
#include <openvdb/tools/LevelSetSphere.h>

// https://www.openvdb.org/documentation/doxygen/codeExamples.html#sAllocatingGrids

int main()
{
  openvdb::initialize();
  // Create a FloatGrid and populate it with a narrow-band
  // signed distance field of a sphere.
  openvdb::FloatGrid::Ptr grid =
      openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
          50.0, openvdb::Vec3f(1.5, 2, 3),
          0.5, 4.0);
  // Associate some metadata with the grid.
  grid->insertMeta("radius", openvdb::FloatMetadata(50.0));
  // Name the grid "density".
  grid->setName("density");
  // Create a VDB file object and write out the grid.
  openvdb::io::File("mygrids.vdb").write({grid});
}