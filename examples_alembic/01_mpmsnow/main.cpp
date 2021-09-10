

#include <iostream>
#include <unordered_map>
#include <array>
#include <random>
#include "delfem2/mat3.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/points.h"
#include "delfem2/grid3hash.h"

// ---------------

struct GridNodeData {
  double mass;
  double velo[3];
  double force[3];
  double velo_old[3];
};

void AddMass(
    std::unordered_map<delfem2::GridCoordinate3, GridNodeData>& grid,
    const delfem2::GridCoordinate3& gc,
    double mass)
{
  const auto it = grid.find(gc);
  if (it == grid.end()) {
    grid[gc] = GridNodeData{
      mass,
      {0,0,0},
      {0,0,0} };
  }
  else {
    grid[gc].mass += mass;
  }
}

void AddVelocity(
    std::unordered_map<delfem2::GridCoordinate3, GridNodeData>& grid,
    const delfem2::GridCoordinate3& gc,
    const double velo[3])
{
  const auto it = grid.find(gc);
  if (it == grid.end()) {
    grid[gc] = GridNodeData{
      0,
      {velo[0],velo[1],velo[2]},
      {0,0,0} };
  }
  else {
    grid[gc].velo[0] += velo[0];
    grid[gc].velo[1] += velo[1];
    grid[gc].velo[2] += velo[2];
  }
}

void AddForce(
    std::unordered_map<delfem2::GridCoordinate3, GridNodeData>& grid,
    const delfem2::GridCoordinate3& gc,
    const double force[3])
{
  const auto it = grid.find(gc);
  if (it == grid.end()) {
    grid[gc] = GridNodeData{
        0,
        {0,0,0},
        {force[0], force[1], force[2]} };
  }
  else {
    grid[gc].force[0] += force[0];
    grid[gc].force[1] += force[1];
    grid[gc].force[2] += force[2];
  }
}

double GetMass(
    const std::unordered_map<delfem2::GridCoordinate3, GridNodeData>& grid,
    const delfem2::GridCoordinate3& gc)
{
  auto it = grid.find(gc);
  if (it != grid.end()) { return it->second.mass; }
  return 0.0;
}

void GetVelocity(
    double velo[3],
    const std::unordered_map<delfem2::GridCoordinate3, GridNodeData>& grid,
    const delfem2::GridCoordinate3& gc)
{
  auto it = grid.find(gc);
  if (it == grid.end()) {
    velo[0] = 0.0;
    velo[1] = 0.0;
    velo[2] = 0.0;
    return;
  }
  velo[0] = it->second.velo[0];
  velo[1] = it->second.velo[1];
  velo[2] = it->second.velo[2];
}

void GetVelocityOld(
    double velo_old[3],
    const std::unordered_map<delfem2::GridCoordinate3, GridNodeData>& grid,
    const delfem2::GridCoordinate3& gc)
{
  auto it = grid.find(gc);
  if (it == grid.end()) {
    velo_old[0] = 0.0;
    velo_old[1] = 0.0;
    velo_old[2] = 0.0;
    return;
  }
  velo_old[0] = it->second.velo_old[0];
  velo_old[1] = it->second.velo_old[1];
  velo_old[2] = it->second.velo_old[2];
}

// ---------------------

void RasterizeParticlesMassWeightOnGrid(
    std::vector< std::unordered_map<delfem2::GridCoordinate3,std::array<double,4>> >& aPtclWdWGrid,
    std::unordered_map<delfem2::GridCoordinate3, GridNodeData>& aGrid,
    const std::vector<double> &aPtclPosition,
    const std::vector<double> &aPtclMass,
    double grid_cell_size,
    const int aGridIndxOffset[64][3])
{
  unsigned int np = aPtclPosition.size() / 3;
  assert(aPtclMass.size() == np);
  assert(aPtclWdWGrid.size() == np);
  const double h_inverse = 1.0 / grid_cell_size;
  for (unsigned int ip = 0; ip < np; ++ip) {
    const double particle_mass = aPtclMass[ip];
    const int ix = static_cast<int>(std::floor(aPtclPosition[ip * 3 + 0] * h_inverse));
    const int iy = static_cast<int>(std::floor(aPtclPosition[ip * 3 + 1] * h_inverse));
    const int iz = static_cast<int>(std::floor(aPtclPosition[ip * 3 + 2] * h_inverse));
    for(int ioff=0;ioff<64;++ioff){
      const delfem2::GridCoordinate3 gc(
          ix + aGridIndxOffset[ioff][0],
          iy + aGridIndxOffset[ioff][1],
          iz + aGridIndxOffset[ioff][2]);
      const double ofx = aPtclPosition[ip * 3 + 0] - static_cast<double>(gc.i) * grid_cell_size;
      const double ofy = aPtclPosition[ip * 3 + 1] - static_cast<double>(gc.j) * grid_cell_size;
      const double ofz = aPtclPosition[ip * 3 + 2] - static_cast<double>(gc.k) * grid_cell_size;
      const double weight = delfem2::WGrid3Cubic(h_inverse, ofx, ofy, ofz);
      AddMass(aGrid, gc, particle_mass * weight);
      aPtclWdWGrid[ip][gc][0] = weight;
      aPtclWdWGrid[ip][gc][1] = delfem2::dWGrid3Cubic(h_inverse, ofx, ofy, ofz);
      aPtclWdWGrid[ip][gc][2] = delfem2::dWGrid3Cubic(h_inverse, ofy, ofz, ofx);
      aPtclWdWGrid[ip][gc][3] = delfem2::dWGrid3Cubic(h_inverse, ofz, ofx, ofy);
    }
  }
}

void RasterizeVelocityOnGrid(
    std::unordered_map<delfem2::GridCoordinate3, GridNodeData>& aGrid,
    const std::vector<double> &aPtclPosition,
    const std::vector<double> &aPtclVelocity,
    const std::vector<double> &aPtclMass,
    const std::vector< std::unordered_map<delfem2::GridCoordinate3,std::array<double,4>> >& aPtclWdWGrid,
    double grid_cell_size,
    const int aGridIndxOffset[64][3])
{
  const unsigned int np = aPtclPosition.size()/3;
  assert(aPtclMass.size()==np);
  assert(aPtclVelocity.size()==np*3);
  const double h_inverse = 1.0 / grid_cell_size;
  for (unsigned int ip=0;ip<np;++ip) {
    const double* velo_ptcl = aPtclVelocity.data()+ip*3;
    const double mass_ptcl = aPtclMass[ip];
    const int ix0 = static_cast<int>(std::floor(aPtclPosition[ip * 3 + 0]*h_inverse));
    const int iy0 = static_cast<int>(std::floor(aPtclPosition[ip * 3 + 1]*h_inverse));
    const int iz0 = static_cast<int>(std::floor(aPtclPosition[ip * 3 + 2]*h_inverse));
    for(int ioff=0;ioff<64;++ioff){
      const  delfem2::GridCoordinate3 gc(
          ix0 + aGridIndxOffset[ioff][0],
          iy0 + aGridIndxOffset[ioff][1],
          iz0 + aGridIndxOffset[ioff][2]);
      const double mass_node = GetMass(aGrid,gc);
      assert(aPtclWdWGrid[ip].find(gc) != aPtclWdWGrid[ip].end() );
      const double weight_node = aPtclWdWGrid[ip].find(gc)->second[0];
      if (mass_node == 0.0 || weight_node == 0.0) { continue; }
      const double velo[3] = {
          velo_ptcl[0] * mass_ptcl * weight_node / mass_node,
          velo_ptcl[1] * mass_ptcl * weight_node / mass_node,
          velo_ptcl[2] * mass_ptcl * weight_node / mass_node };
      AddVelocity(aGrid,gc,velo);
    }
  }
}


void VolumeParticle(
    std::vector<double> &aPtclVolume0,
    const std::unordered_map<delfem2::GridCoordinate3, GridNodeData>& aGrid,
    const std::vector< std::unordered_map<delfem2::GridCoordinate3,std::array<double,4>> >& aPtclWdWGrid,
    double grid_cell_size,
    const std::vector<double> &aPtclPos,
    const std::vector<double> &aPtclMass,
    const int aGridIndxOffset[64][3])
{
  const unsigned int np = aPtclPos.size()/3;
  assert(aPtclMass.size()==np);
  aPtclVolume0.resize(np);
  const double h_inverse = 1.0 / grid_cell_size;
  for (unsigned int ip=0;ip<np;++ip) {
    double density0 = 0.0;
    const int ix0 = static_cast<int>(std::floor( aPtclPos[ip*3+0]*h_inverse ) );
    const int iy0 = static_cast<int>(std::floor( aPtclPos[ip*3+1]*h_inverse ) );
    const int iz0 = static_cast<int>(std::floor( aPtclPos[ip*3+2]*h_inverse ) );
    for(int ioff=0;ioff<64;++ioff){
      const delfem2::GridCoordinate3 gc(
          ix0 + aGridIndxOffset[ioff][0],
          iy0 + aGridIndxOffset[ioff][1],
          iz0 + aGridIndxOffset[ioff][2]);
      const double node_mass = GetMass(aGrid, gc);
      const double cached_weight = aPtclWdWGrid[ip].find(gc)->second[0];
      if (node_mass == 0.0 || cached_weight == 0.0) { continue; }
      density0 += node_mass * cached_weight;
    }
    density0 /= std::pow(grid_cell_size, 3);
    aPtclVolume0[ip] = aPtclMass[ip] / density0;
  }
}

void UpdateGridVelocity(
    std::unordered_map<delfem2::GridCoordinate3, GridNodeData>& aGrid,
    double dt)
{
  for (auto& node : aGrid){
    GridNodeData& data = node.second;
    data.velo_old[0] = data.velo[0];
    data.velo_old[1] = data.velo[1];
    data.velo_old[2] = data.velo[2];
    if (data.mass == 0.0) {
      data.velo[0] = 0.0;
      data.velo[1] = 0.0;
      data.velo[2] = 0.0;
    }
    else {
      data.velo[0] = data.velo[0] + dt * data.force[0] / data.mass;
      data.velo[1] = data.velo[1] + dt * data.force[1] / data.mass;
      data.velo[2] = data.velo[2] + dt * data.force[2] / data.mass;
    }
  }
}

void UpdateParticleVelocity(
    std::vector<double> &aPtclVelo,
    const std::vector<double> &aPtclPos,
    std::unordered_map<delfem2::GridCoordinate3, GridNodeData>& aGrid,
    const std::vector< std::unordered_map<delfem2::GridCoordinate3,std::array<double,4>> >& aPtclWdWGrid,
    double grid_cell_size,
    double flip_alpha,
    const int aGridIndxOffset[64][3])
{
  const unsigned int np = aPtclPos.size()/3;
  assert( aPtclPos.size() == np*3 );
  double vtmp[3], votmp[3];
  for (unsigned int ip=0;ip<np;++ip) {
    double v_PIC[3] = {0.0, 0.0, 0.0};
    double v_FLIP[3] = {aPtclVelo[ip*3+0], aPtclVelo[ip*3+1], aPtclVelo[ip*3+2]};
    const int ix0 = static_cast<int>(std::floor( aPtclPos[ip*3+0]/ grid_cell_size ));
    const int iy0 = static_cast<int>(std::floor( aPtclPos[ip*3+1]/ grid_cell_size ));
    const int iz0 = static_cast<int>(std::floor( aPtclPos[ip*3+2]/ grid_cell_size ));
    for(int ioff=0;ioff<64;++ioff){
      const delfem2::GridCoordinate3 gc(
          ix0 + aGridIndxOffset[ioff][0],
          iy0 + aGridIndxOffset[ioff][1],
          iz0 + aGridIndxOffset[ioff][2]);
      const double weight_ip = aPtclWdWGrid[ip].find(gc)->second[0];
      GetVelocity(vtmp, aGrid, gc);
      GetVelocityOld(votmp, aGrid, gc);
      v_PIC[0] += vtmp[0] * weight_ip;
      v_PIC[1] += vtmp[1] * weight_ip;
      v_PIC[2] += vtmp[2] * weight_ip;
      v_FLIP[0] += (vtmp[0] - votmp[0]) * weight_ip;
      v_FLIP[1] += (vtmp[1] - votmp[1]) * weight_ip;
      v_FLIP[2] += (vtmp[2] - votmp[2]) * weight_ip;
    }
    aPtclVelo[ip*3+0] = (1 - flip_alpha) * v_PIC[0] + flip_alpha * v_FLIP[0];
    aPtclVelo[ip*3+1] = (1 - flip_alpha) * v_PIC[1] + flip_alpha * v_FLIP[1];
    aPtclVelo[ip*3+2] = (1 - flip_alpha) * v_PIC[2] + flip_alpha * v_FLIP[2];
  }
}

void AddSphere(
    std::vector<double>& aPosPtcl,
    double radius,
    unsigned int ndiv,
    const std::array<double,3>& position,
    std::mt19937& gen)
{
  const double dDia = radius * 2.0 / static_cast<double>(ndiv);
  std::uniform_real_distribution<double> dist(-dDia * 0.2, +dDia * 0.2);
  for (unsigned int ix = 0; ix < ndiv; ++ix) {
    for (unsigned int iy = 0; iy < ndiv; ++iy) {
      for (unsigned int iz = 0; iz < ndiv; ++iz) {
        const double px = dDia * ix - radius + dist(gen);
        const double py = dDia * iy - radius + dist(gen);
        const double pz = dDia * iz - radius + dist(gen);
        if ( px*px+py*py+pz*pz > radius * radius) { continue; }
        aPosPtcl.push_back(position[0]+px);
        aPosPtcl.push_back(position[1]+py);
        aPosPtcl.push_back(position[2]+pz);
      }
    }
  }
}

template <typename T>
T MyClamp(
    T vin,
    T vmin,
    T vmax )
{
  if( vin < vmin ){ return vmin; }
  if( vin > vmax ){ return vmax; }
  return vin;
}

// -----------------------------------

class CParticleData {
public:
  double def_elastic[9];
  double def_plastic[9];
};

void ForceOnGrid(
    std::unordered_map<delfem2::GridCoordinate3, GridNodeData>& aGrid,
    const std::vector<double> &aPtclPos,
    const std::vector<double> &aPtclVolume0,
    const std::vector<CParticleData>& aPtclData,
    const std::vector< std::unordered_map<delfem2::GridCoordinate3,std::array<double,4>> >& aPtclWdWGrid,
    double dt,
    double grid_cell_size,
    double mu,
    double lambda,
    double hardening_coefficient,
    const double gravity[3],
    const int aGridIndxOffset[64][3])
{
  const unsigned int np = aPtclPos.size()/3;
  const double h_inverse = 1.0 / grid_cell_size;
  for (unsigned int ip=0;ip<np;++ip) {
    const int ix0 = static_cast<int>(std::floor( aPtclPos[ip*3+0]*h_inverse ));
    const int iy0 = static_cast<int>(std::floor( aPtclPos[ip*3+1]*h_inverse ));
    const int iz0 = static_cast<int>(std::floor( aPtclPos[ip*3+2]*h_inverse ));
    const delfem2::CMat3d Fp(aPtclData[ip].def_plastic);
    const delfem2::CMat3d Fe(aPtclData[ip].def_elastic);
    const double detFp = Fp.determinant();
    const double volume = detFp * aPtclVolume0[ip];
    double sigma[9];
    {
      const double detFe = Fe.determinant();
      const double detF = detFe * detFp;
      delfem2::CMat3d Rot;
      {
        double U0[9], G0[3], V0[9];
        delfem2::svd3(U0,G0,V0,Fe.data(),30);
        delfem2::MatMatT3(Rot.mat, U0,V0);
      }
      const double mu0 = mu * std::exp(hardening_coefficient * (1.0 - detFp));
      const double lambda0 = lambda * std::exp(hardening_coefficient * (1.0 - detFp));
      const delfem2::CMat3d EigSigma =
          2.0 * mu0 / detF * (Fe - Rot) * Fe.transpose() +
          lambda0 / detF * (detFe - 1.0) * detFe * delfem2::CMat3d::Identity();
      EigSigma.CopyTo(sigma);
    }
    for(int ioff=0;ioff<64;++ioff){
      const delfem2::GridCoordinate3 gc(
          ix0 + aGridIndxOffset[ioff][0],
          iy0 + aGridIndxOffset[ioff][1],
          iz0 + aGridIndxOffset[ioff][2]);
      const double dw[3] = {
          aPtclWdWGrid[ip].find(gc)->second[1],
          aPtclWdWGrid[ip].find(gc)->second[2],
          aPtclWdWGrid[ip].find(gc)->second[3] };
      double force[3];
      {
        delfem2::MatVec3(force, sigma, dw);
        force[0] *= -volume;
        force[1] *= -volume;
        force[2] *= -volume;
      }
      AddForce(aGrid, gc, force);
    }
  }
  // Add gravity forces
  for (auto &node : aGrid) {
    GridNodeData &data = node.second;
    data.force[0] += data.mass * gravity[0];
    data.force[1] += data.mass * gravity[1];
    data.force[2] += data.mass * gravity[2];
  }
}

void UpdateDeformationGradient(
    CParticleData& pData,
    const std::array<double,3>& pos,
    std::unordered_map<delfem2::GridCoordinate3, GridNodeData>& aGrid,
    std::unordered_map<delfem2::GridCoordinate3,std::array<double,4>>& aWdWGrid,
    double dt,
    double grid_cell_size,
    double critical_compression,
    double critical_stretch,
    const int aGridIndxOffset[64][3])
{
  const double h_inverse = 1.0 / grid_cell_size;
  delfem2::CMat3d gradV;
  gradV.setZero();
  const int ix0 = static_cast<int>(std::floor( pos[0]*h_inverse ));
  const int iy0 = static_cast<int>(std::floor( pos[1]*h_inverse ));
  const int iz0 = static_cast<int>(std::floor( pos[2]*h_inverse ));
  for(int ioff=0;ioff<64;++ioff){
    const delfem2::GridCoordinate3 gc(
        ix0+aGridIndxOffset[ioff][0],
        iy0+aGridIndxOffset[ioff][1],
        iz0+aGridIndxOffset[ioff][2]);
    double dw[3] = {0,0,0};
    if(aWdWGrid.find(gc) != aWdWGrid.end() ) {
      dw[0] = aWdWGrid.find(gc)->second[1];
      dw[1] = aWdWGrid.find(gc)->second[2];
      dw[2] = aWdWGrid.find(gc)->second[3];
    }
    double vtmp[3]; GetVelocity(vtmp, aGrid, gc);
    for(int i=0;i<3;++i){
      for(int j=0;j<3;++j){
        gradV.mat[i*3+j] += vtmp[i]*dw[j];
      }
    }
  }
  const delfem2::CMat3d Fe0(pData.def_elastic);
  const delfem2::CMat3d Fp0(pData.def_plastic);
  const delfem2::CMat3d Fe1 = (delfem2::CMat3d::Identity()+dt*delfem2::CMat3d(gradV))*Fe0;
  delfem2::CMat3d U, V;
  double G[3];
  delfem2::svd3(U.mat,G,V.mat,Fe1.mat,30);
  const double a = 1 - critical_compression;
  const double b = 1 + critical_stretch;
  G[0] = MyClamp(G[0], a, b);
  G[1] = MyClamp(G[1], a, b);
  G[2] = MyClamp(G[2], a, b);
  delfem2::CMat3d Fe2 = U * delfem2::CMat3d(G[0],G[1],G[2]) * V.transpose();
  delfem2::CMat3d Fp2 = V * delfem2::CMat3d(1./G[0],1./G[1],1./G[2]) * U.transpose() * Fe1 * Fp0;
  Fe2.CopyTo(pData.def_elastic);
  Fp2.CopyTo(pData.def_plastic);
}

template <typename T>
T Dot3(
    const T* v0,
    const T* v1)
{
  return v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
}

void ResolveCollision(
    double velocity[3],
    const double position[3],
    double dt,
    double friction_coeff,
    const double normal[3])
{
  const double p2[3] = {
      position[0] + dt * velocity[0],
      position[1] + dt * velocity[1],
      position[2] + dt * velocity[2] };

  if (p2[1]>0) { return; }

  const double v_n = Dot3(normal, velocity );
  if (v_n >= 0) { return; }

  const double v_t[3] = {
      velocity[0] - v_n * normal[0],
      velocity[1] - v_n * normal[1],
      velocity[2] - v_n * normal[2] };

  const double v_t_norm = sqrt( v_t[0]*v_t[0] + v_t[1]*v_t[1] + v_t[2]*v_t[2] );

  if (v_t_norm <= -friction_coeff * v_n) {
    velocity[0] = 0.0;
    velocity[1] = 0.0;
    velocity[2] = 0.0;
    return;
  }

  velocity[0] = v_t[0] + friction_coeff * v_n / v_t_norm * v_t[0];
  velocity[1] = v_t[1] + friction_coeff * v_n / v_t_norm * v_t[1];
  velocity[2] = v_t[2] + friction_coeff * v_n / v_t_norm * v_t[2];
}

// -------------------------------------------------

int main(){
  std::vector<double> aPtclPosition;
  {
//    std::mt19937 gen(static_cast<std::mt19937::result_type>(std::random_device()()));
    std::mt19937 gen(0);
    AddSphere(
        aPtclPosition,
        1, 15,
        std::array<double, 3>{-4, 4, 0},
        gen);
  }
  // ----------
  std::vector<double> aPtclMass, aPtclVelocity;
  {
    const unsigned int np = aPtclPosition.size()/3;
    assert(aPtclPosition.size()==np*3);
    for(int ip=0;ip<np;++ip) {
      aPtclVelocity.push_back(5.0);
      aPtclVelocity.push_back(1.0);
      aPtclVelocity.push_back(0.0);
      aPtclMass.push_back(0.2);
    }
  }
  // ------
  int aGridIndxOffset[64][3];
  {
    int icnt = 0;
    for (int a = -1; a < 3; ++a) {
      for (int b = -1; b < 3; ++b) {
        for (int c = -1; c < 3; ++c) {
          aGridIndxOffset[icnt][0] = a;
          aGridIndxOffset[icnt][1] = b;
          aGridIndxOffset[icnt][2] = c;
          ++icnt;
        }
      }
    }
    assert(icnt==64);
  }
  std::vector<double> aPtclVolume0;
  std::vector<CParticleData> aPtclData;
  aPtclData.resize(aPtclPosition.size()/3);
  for (auto &data : aPtclData){
    delfem2::Mat3_Identity(data.def_elastic, 1.0);
    delfem2::Mat3_Identity(data.def_plastic, 1.0);
  }
  std::vector< std::unordered_map<delfem2::GridCoordinate3,std::array<double,4>> > aPtclWdWGrid;
  aPtclWdWGrid.resize(aPtclPosition.size()/3);
  std::unordered_map<delfem2::GridCoordinate3, GridNodeData> aGrid;
  const unsigned int NSTEP_PER_FRAME = 20;
  const double grid_cell_size = 0.1;
  const double dt = 1.0/1200;
  const double youngs_modulus = 1.4e5;
  const double poissons_ratio = 0.2;
  const double mu = (youngs_modulus / (2 * (1 + poissons_ratio)));
  const double lambda = (youngs_modulus * poissons_ratio / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio)));
  const double hardening_coefficient = 10;
  const double critical_compression = 2.5e-2;
  const double critical_stretch = 7.5e-3;
  const double gravity[3] = {0.0, -9.81, 0.0};
  const double flip_alpha = 0.95;
  const double friction_coeff = 0.35;
  const double normal[3] = {0.0, 1.0, 0.0};
  // ------------
  for(int itr=0;itr<120*NSTEP_PER_FRAME;++itr) {
    std::cout << itr/NSTEP_PER_FRAME << " " << itr-(itr/NSTEP_PER_FRAME)*NSTEP_PER_FRAME << std::endl;
    // ---
    aGrid.clear();
    for (auto &wdw : aPtclWdWGrid) { wdw.clear(); }
    // ---
    RasterizeParticlesMassWeightOnGrid(
        aPtclWdWGrid, aGrid,
        aPtclPosition, aPtclMass,
        grid_cell_size,
        aGridIndxOffset);
    RasterizeVelocityOnGrid(
        aGrid,
        aPtclPosition, aPtclVelocity, aPtclMass,
        aPtclWdWGrid,
        grid_cell_size,
        aGridIndxOffset);
    if( itr == 0 ) {
      VolumeParticle(
          aPtclVolume0,
          aGrid, aPtclWdWGrid, grid_cell_size,
          aPtclPosition, aPtclMass, aGridIndxOffset);
    }
    ForceOnGrid(
        aGrid,
        aPtclPosition, aPtclVolume0, aPtclData, aPtclWdWGrid,
        dt, grid_cell_size,
        mu, lambda, hardening_coefficient, gravity,
        aGridIndxOffset);
    UpdateGridVelocity(
        aGrid, dt);
    for(auto& node : aGrid){
      const double position[3] = {
          node.first.i * grid_cell_size,
          node.first.j * grid_cell_size,
          node.first.k * grid_cell_size };
      ResolveCollision(
          node.second.velo,
          position, dt, friction_coeff, normal);
    }
    //
    for (unsigned int ip = 0; ip < aPtclPosition.size() / 3; ++ip) {
      UpdateDeformationGradient(
          aPtclData[ip],
          {aPtclPosition[ip * 3 + 0], aPtclPosition[ip * 3 + 1], aPtclPosition[ip * 3 + 2]},
          aGrid, aPtclWdWGrid[ip], dt, grid_cell_size,
          critical_compression,
          critical_stretch,
          aGridIndxOffset);
    }
    //
    UpdateParticleVelocity(
        aPtclVelocity,
        aPtclPosition, aGrid, aPtclWdWGrid,
        grid_cell_size, flip_alpha, aGridIndxOffset);
    //
    for(unsigned int ip=0;ip<aPtclPosition.size()/3;++ip){
      ResolveCollision(
          aPtclVelocity.data()+ip*3,
          aPtclPosition.data()+ip*3, dt, friction_coeff, normal);
    }
    //
    for (unsigned int ip = 0; ip < aPtclPosition.size() / 3; ++ip) {
      aPtclPosition[ip * 3 + 0] += dt * aPtclVelocity[ip * 3 + 0];
      aPtclPosition[ip * 3 + 1] += dt * aPtclVelocity[ip * 3 + 1];
      aPtclPosition[ip * 3 + 2] += dt * aPtclVelocity[ip * 3 + 2];
    }
    if( itr % NSTEP_PER_FRAME == 0 ){
      std::vector< std::pair< std::vector<double>, std::vector<unsigned int> > > aMesh;
      for(unsigned ip=0;ip<aPtclPosition.size()/3;++ip) {
        std::vector<double> aXYZ;
        std::vector<unsigned int> aTri;
        delfem2::MeshTri3D_Sphere(aXYZ,aTri,0.1,2,4);
        delfem2::Translate_Points3(aXYZ,
            aPtclPosition[ip*3+0],aPtclPosition[ip*3+1],aPtclPosition[ip*3+2]);
        aMesh.emplace_back(aXYZ,aTri );
      }
      delfem2::Write_Obj(std::to_string(itr/NSTEP_PER_FRAME)+".obj",aMesh);
    }
  }
}
