

/**
 * @file demo inspired by taichi_mpm (https://github.com/yuanming-hu/taichi_mpm)
 */

#include <cstdlib>
#include <cstdio>
#include <random>
#include <algorithm>  // std::max
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/color.h"
#include "delfem2/vec2.h"
#include "delfem2/vec3.h"
#include "delfem2/mat2.h"
#include "delfem2/glfw/viewer2.h"
#include "delfem2/glfw/util.h"

class CParticle {
 public:
  CParticle(const delfem2::CVec2f &x,
            int c)
      : pos(x), velo(0.f), F(1), C(0), Jp(1), icolor(c) {}
 public:
  delfem2::CVec2f pos;   //! position
  delfem2::CVec2f velo;  //! velocity
  delfem2::CMat2f F;  //! deformation gradient
  delfem2::CMat2f C;  //! Cauchy stress
  float Jp;  //! volume change ratio
  int icolor;  //! integer RGB color
};

template<typename T>
T myclamp(T v, T vmin, T vmax) {
  if (v < vmin) { return vmin; }
  if (v > vmax) { return vmax; }
  return v;
}

void MPM_Particle2Grid2(
    std::vector<delfem2::CVec3f> &aVeloGrid,
    const std::vector<CParticle> &aParticles,
    float dt,
    int ngrid,
    const float hgrid,
    float vol,
    float particle_mass,
    float mu_0,
    float lambda_0,
    float hardening) {
  const unsigned int mgrid = ngrid + 1;
  const float inv_dx = 1.f / hgrid;
  for (auto &g : aVeloGrid) { g.setZero(); }
  for (auto &p : aParticles) {
    delfem2::CVec2i bc = (p.pos * inv_dx - delfem2::CVec2f(0.5)).cast<int>();
    delfem2::CVec2f fx0 = p.pos * inv_dx - bc.cast<float>();
    const float dfx = fx0.p[0];
    const float dfy = fx0.p[1];
    const delfem2::CVec2f w0[3]{
        {0.5f * (1.5f - dfx) * (1.5f - dfx), 0.5f * (1.5f - dfy) * (1.5f - dfy)},
        {0.75f - (dfx - 1.f) * (dfx - 1.f), 0.75f - (dfy - 1.f) * (dfy - 1.f)},
        {0.5f * (dfx - 0.5f) * (dfx - 0.5f), 0.5f * (dfy - 0.5f) * (dfy - 0.5f)}};
    delfem2::CMat2f affine;
    {
      auto e = std::exp(hardening * (1.0f - p.Jp));
      auto mu = mu_0 * e;
      auto lambda = lambda_0 * e;
      const float J = p.F.determinant();
      delfem2::CMat2f r, s;
      polar_decomposition(
          r, s,
          p.F);
      const delfem2::CMat2f cauchy = 2 * mu * (p.F - r) * p.F.transpose() +
          lambda * (J - 1) *
              J;
      const delfem2::CMat2f stress = -4 * inv_dx * inv_dx * dt * vol * cauchy;
      affine = stress + particle_mass * p.C;
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {  // Scatter to grid
        delfem2::CVec2f dpos0((float(i) - fx0.x) * hgrid, (float(j) - fx0.y) * hgrid);
        const delfem2::CVec3f mv(
            p.velo.x * particle_mass,
            p.velo.y * particle_mass,
            particle_mass);  //translational momentum of particle, and mass
        const unsigned int ig = bc.x + i;
        assert(ig < mgrid);
        const unsigned int jg = bc.y + j;
        assert(jg < mgrid);
        delfem2::CVec2f t1;
        affine.multiply_vec2(t1.p, dpos0.p);
        const auto t0 = mv + delfem2::CVec3f(t1.x, t1.y, 0);
        aVeloGrid[ig * mgrid + jg] += (w0[i].x * w0[j].y) * t0;
      }
    }
  }
}

void MPM_Particle2Grid_Snow2(
    std::vector<CParticle> &aParticles,
    const std::vector<delfem2::CVec3f> &aVeloGrid,
    float dt,
    int ngrid,
    const float hgrid,
    bool is_plastic) {
  const unsigned int mgrid = ngrid + 1;
  const float inv_dx = 1.f / hgrid;
  for (auto &particle : aParticles) {
    delfem2::CVec2i bc = (particle.pos * inv_dx - delfem2::CVec2f(0.5)).cast<int>();
    delfem2::CVec2f fx0 = particle.pos * inv_dx - bc.cast<float>();
    const float dfx = fx0.p[0];
    const float dfy = fx0.p[1];
    const delfem2::CVec2f w0[3]{
        {0.5f * (1.5f - dfx) * (1.5f - dfx), 0.5f * (1.5f - dfy) * (1.5f - dfy)},
        {0.75f - (dfx - 1.f) * (dfx - 1.f), 0.75f - (dfy - 1.f) * (dfy - 1.f)},
        {0.5f * (dfx - 0.5f) * (dfx - 0.5f), 0.5f * (dfy - 0.5f) * (dfy - 0.5f)}};
    particle.C.setZero();
    particle.velo.setZero();
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        const delfem2::CVec2f dpos(float(i) - fx0.x, float(j) - fx0.y);
        const unsigned int ig = bc.x + i;
        assert(ig < mgrid);
        const unsigned int jg = bc.y + j;
        assert(jg < mgrid);
        const delfem2::CVec2f grid_v(
            aVeloGrid[ig * mgrid + jg].x,
            aVeloGrid[ig * mgrid + jg].y);
        const float weight = w0[i].x * w0[j].y;
        particle.velo += weight * delfem2::CVec2f(grid_v.x, grid_v.y);
        particle.C += 4 * inv_dx * weight * delfem2::CMat2f::outer_product(grid_v.p, dpos.p);
      }
    }
    particle.pos += dt * particle.velo;
    const delfem2::CMat2f tmp0 = delfem2::CMat2f(1.) + dt * particle.C;
    delfem2::CMat2f F0 = tmp0 * particle.F;
    const float oldJ0 = F0.determinant();
    if (is_plastic) {  // updating deformation gradient tensor by clamping the eignvalues
      delfem2::CMat2f svd_u0, sig0, svd_v0;
      delfem2::svd<float>(
          svd_u0, sig0, svd_v0,
          F0);
      for (int i = 0; i < 2; i++) {  // Snow Plasticity
        sig0(i, i) = myclamp(sig0(i, i), 1.0f - 2.5e-2f, 1.0f + 7.5e-3f);
      }
      F0 = svd_u0 * sig0 * svd_v0.transpose();
    }
    float Jp_new0 = myclamp(particle.Jp * oldJ0 / F0.determinant(), 0.6f, 20.0f);
    particle.Jp = Jp_new0;
    particle.F = F0;

  }
}

void StepTime_Mpm2Snow(
    std::vector<CParticle> &aParticles,
    std::vector<delfem2::CVec3f> &aVeloGrid,
    float dt,
    int ngrid,
    const float hgrid,
    //
    float vol,
    float particle_mass,
    float mu_0,
    float lambda_0,
    float hardening,
    bool plastic) {
  MPM_Particle2Grid2(
      aVeloGrid,
      aParticles, dt, ngrid, hgrid,
      vol, particle_mass,
      mu_0, lambda_0, hardening);
  const unsigned int mgrid = ngrid + 1;
  for (unsigned int igrid = 0; igrid < mgrid; igrid++) {
    for (unsigned int jgrid = 0; jgrid < mgrid; jgrid++) {  // For all grid nodes
      auto &grid = aVeloGrid[igrid * mgrid + jgrid];
      if (grid[2] <= 0) { continue; } // grid is empty
      grid /= grid[2];                                   // Normalize by mass
      grid += dt * delfem2::CVec3f(0, -200, 0); // Gravity
      const float boundary = 0.05f;
      const float x = float(igrid) / float(ngrid);
      const float y = float(jgrid) / float(ngrid); // boundary thick.,node coord
      if (x < boundary || x > 1 - boundary || y > 1 - boundary) {
        grid = delfem2::CVec3f(0.f, 0.f, 0.f); // Sticky
      }
      if (y < boundary) {
        grid[1] = std::max(0.0f, grid[1]);  // "Separate"
      }
    }
  }
  MPM_Particle2Grid_Snow2(
      aParticles,
      aVeloGrid, dt, ngrid, hgrid, plastic);
}

void AddParticlesCircle(
    std::vector<CParticle> &particles,
    const delfem2::CVec2f &center,
    int icolor) {   // Seed particles with position and color
  std::mt19937 rnd_eng(0);
  std::uniform_real_distribution<float> dist01(0, 1);
  for (int i = 0; i < 500; i++) {
    float x0 = dist01(rnd_eng) * 2.0f - 1;
    float y0 = dist01(rnd_eng) * 2.0f - 1;
    if (x0 * x0 + y0 * y0 > 1) continue;
    particles.emplace_back(
        delfem2::CVec2f(x0, y0) * 0.08f + center,
        icolor);
  }
}

// ---------------------------
// OpenGL dependency from here

static void error_callback(
    [[maybe_unused]] int error, const char *description) {
  fputs(description, stderr);
}

int main() {
  std::vector<CParticle> aParticles;
  AddParticlesCircle(aParticles, delfem2::CVec2f(0.3f, 0.8f), 0xFF0000);
  AddParticlesCircle(aParticles, delfem2::CVec2f(0.4f, 0.6f), 0x00FF00);
  AddParticlesCircle(aParticles, delfem2::CVec2f(0.5f, 0.8f), 0x0000FF);
  const int n = 80;
  const float dx = 1.0f / n;
  std::vector<delfem2::CVec3f> aVeloGrid((n + 1) * (n + 1));
  const float dt = 1e-4f;
  float particle_mass = 1.0f;
  float vol = 1.0f;
  float hardening = 10.0f;
  float E = 1e4f;
  float nu = 0.2f;
  float mu_0 = E / (2 * (1 + nu));
  float lambda_0 = E * nu / ((1 + nu) * (1 - 2 * nu));
  bool is_plastic = true;
  //
  delfem2::glfw::CViewer2 viewer;
  {
    viewer.view_height = 1.0;
    viewer.trans[0] -= 0.5;
    viewer.trans[1] -= 0.5;
    viewer.title = "633_MpmSnow2";
  }
  glfwSetErrorCallback(error_callback);
  if (!glfwInit()) { exit(EXIT_FAILURE); }
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  // ------
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();

  while (!glfwWindowShouldClose(viewer.window)) {
    StepTime_Mpm2Snow(aParticles, aVeloGrid,
                      dt, n, dx,
                      vol, particle_mass, mu_0, lambda_0, hardening, is_plastic);  //  Advance simulation
    //
    viewer.DrawBegin_oldGL();
    ::glDisable(GL_LIGHTING);
    ::glPointSize(5);
    ::glBegin(GL_POINTS);
    for (const auto &p : aParticles) {
      double rgb[3];
      delfem2::ColorRGB_Int(rgb, p.icolor);
      ::glColor3dv(rgb);
      ::glVertex2d(p.pos[0], p.pos[1]);
    }
    ::glVertex2d(0.0, 0.0);
    ::glEnd();
    ::glColor3d(0.0, 0.0, 0.0);
    ::glBegin(GL_LINES);
    ::glVertex2d(0, 0);
    ::glVertex2d(1, 0);
    ::glEnd();
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
