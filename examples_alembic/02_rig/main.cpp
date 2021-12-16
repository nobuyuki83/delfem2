#include <stack>
#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreOgawa/All.h>

#include "delfem2/msh_io_obj.h"
#include "delfem2/msh_rig.h"
#include "delfem2/rig_bvh.h"

int main() {
  namespace dfm2 = delfem2;
  std::vector<dfm2::CRigBone> aBone;
  std::vector<dfm2::CChannel_BioVisionHierarchy> aChannelRotTransBone;
  size_t nframe = 0;
  std::vector<double> aValRotTransBone;

  std::string path_bvh = std::string(PATH_INPUT_DIR) + "/jump.bvh";
  {
    double frame_time;
    std::string header_bvh;
    Read_BioVisionHierarchy(
      aBone, aChannelRotTransBone, nframe, frame_time, aValRotTransBone, header_bvh,
      path_bvh);
  }
  UpdateBoneRotTrans(aBone);

  const float dt = 1.0f / 60.0f;
  // -------------------------
  Alembic::Abc::OArchive archive(Alembic::AbcCoreOgawa::WriteArchive(), "simple.abc");
  Alembic::AbcGeom::OPolyMesh mesh_obj(Alembic::Abc::OObject(archive, Alembic::Abc::kTop), "mesh");
  { // set time sampling
    const Alembic::Abc::TimeSampling time_sampling(dt, 0);
    const uint32_t time_sampling_index = archive.addTimeSampling(time_sampling);
    mesh_obj.getSchema().setTimeSampling(time_sampling_index);
  }
  for(unsigned int iframe=0;iframe<nframe;++iframe){ // set geometry
    const size_t nch = aChannelRotTransBone.size();
    SetPose_BioVisionHierarchy(
      aBone, aChannelRotTransBone,
      aValRotTransBone.data() + iframe * nch);
    std::vector<int> aTri;
    std::vector<float> aXYZ;
    dfm2::Mesh_RigBones_Octahedron<float>(
      aXYZ, aTri,
      aBone);
    /*
    dfm2::Write_Obj<float>(
      "hoge"+std::to_string(iframe)+".obj",
      aXYZ.data(), aXYZ.size() / 3,
      aTri.data(), aTri.size() / 3);
      */
    if( iframe == 0 ) {
      std::vector<int> aElmSize(aTri.size() / 3, 3);
      const Alembic::AbcGeom::OPolyMeshSchema::Sample mesh_samp(
      Alembic::AbcGeom::V3fArraySample((const Alembic::Abc::V3f *) aXYZ.data(), aXYZ.size() / 3),
        Alembic::AbcGeom::Int32ArraySample(aTri.data(), aTri.size()),
        Alembic::AbcGeom::Int32ArraySample(aElmSize.data(), aElmSize.size()));
      mesh_obj.getSchema().set(mesh_samp);
    }
    else{
      const Alembic::AbcGeom::OPolyMeshSchema::Sample mesh_samp(
      Alembic::AbcGeom::V3fArraySample((const Alembic::Abc::V3f *) aXYZ.data(), aXYZ.size() / 3));
      mesh_obj.getSchema().set(mesh_samp);
    }
  }
  return 0;
}
