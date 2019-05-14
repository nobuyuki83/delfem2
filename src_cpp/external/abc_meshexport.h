#ifndef abc_meshexport_h
#define abc_meshexport_h

#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreOgawa/All.h>

class CAlembic_ExportMesh
{
public:
  CAlembic_ExportMesh(const std::string& fname, double dt){
    Alembic::AbcGeom::OArchive archive(Alembic::AbcCoreOgawa::WriteArchive(), fname);
    const Alembic::AbcGeom::TimeSampling time_sampling(dt, 0);
    const uint32_t time_sampling_index = archive.addTimeSampling(time_sampling);
    mesh_obj = Alembic::AbcGeom::OPolyMesh(Alembic::AbcGeom::OObject(archive, Alembic::AbcGeom::kTop), "mesh");
    Alembic::AbcGeom::OPolyMeshSchema& mesh = mesh_obj.getSchema();
    mesh.setTimeSampling(time_sampling_index);
  }
  void Add_Vtx3D_TopoElm
  (const double* paXYZd, size_t nP,
   const int32_t* paElm, size_t nE,
   int nNoEl)
  {
    Alembic::AbcGeom::OPolyMeshSchema& mesh = mesh_obj.getSchema();
    std::vector<float> aXYZf(paXYZd,paXYZd+nP*3);
    auto vrt = Alembic::AbcGeom::V3fArraySample((const Alembic::AbcGeom::V3f*) aXYZf.data(), nP);
    auto ind = Alembic::AbcGeom::Int32ArraySample(paElm, nE*nNoEl);
    std::vector<int> counts(nE,nNoEl);
    auto cnt = Alembic::AbcGeom::Int32ArraySample(counts.data(), nE);
    const Alembic::AbcGeom::OPolyMeshSchema::Sample mesh_samp_1(vrt,ind,cnt);
    mesh.set(mesh_samp_1);
  }
  void Add_Vtx3D
  (const double* paXYZd, size_t nP)
  {
    Alembic::AbcGeom::OPolyMeshSchema& mesh = mesh_obj.getSchema();
    std::vector<float> aXYZf(paXYZd,paXYZd+nP*3);
    auto vrt = Alembic::AbcGeom::V3fArraySample((const Alembic::AbcGeom::V3f*) aXYZf.data(), nP);
    const Alembic::AbcGeom::OPolyMeshSchema::Sample mesh_samp_1(vrt);
    mesh.set(mesh_samp_1);
  }
public:
  Alembic::AbcGeom::OPolyMesh mesh_obj;
};


#endif /* abc_meshexport_h */
