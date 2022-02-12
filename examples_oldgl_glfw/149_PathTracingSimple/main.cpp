
/**
 * @file demo based on smallpt, a Path Tracer by Kevin Beason, 2008
 * (https://www.kevinbeason.com/smallpt/)
 */

#include <cmath>
#include <cstdlib>
#include <vector>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/sampling.h"
#include "delfem2/opengl/tex.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/vec3.h"
#include "delfem2/thread.h"

struct Ray {
  delfem2::CVec3d o, d;
  Ray(const delfem2::CVec3d &o_, const delfem2::CVec3d &d_) : o(o_), d(d_) {}
};

enum Refl_t {
  DIFF, SPEC, REFR
};

class SphereForRayTracing {
 public:
  SphereForRayTracing(
      double rad_,
      const delfem2::CVec3d &p_,
      const delfem2::CVec3d &e_,
      const delfem2::CVec3d &c_,
      Refl_t refl_)
      : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
  [[nodiscard]] double intersect(const Ray &r) const {
    static constexpr double eps = 1e-4;
    const delfem2::CVec3d op = p - r.o;
    double b = op.dot(r.d);
    double det = b * b - op.dot(op) + rad * rad;
    if (det < 0) { return 0; } // no intersection
    det = sqrt(det);
    //
    if (b - det > eps) { return b - det; }
    else if (b + det > eps) { return b + det; }
    else { return 0; }
  }
 public:
  double rad;
  delfem2::CVec3d p, e, c;
  Refl_t refl;
};

template<typename VEC>
std::tuple<bool, VEC, VEC, VEC, VEC, Refl_t> RayIntersection(
    const Ray &r,
    const std::vector<SphereForRayTracing> &spheres) {
  double dist_nearest;
  int isphere_nearest;
  double inf = dist_nearest = 1e20;
  for (int isphere = 0; isphere < spheres.size(); ++isphere) {
    const double d = spheres[isphere].intersect(r);
    if ((d > 0) && d < dist_nearest) {
      dist_nearest = d;
      isphere_nearest = isphere;
    }
  }
  if (dist_nearest == inf) {
    return {false, VEC(), VEC(), VEC(), VEC(), DIFF};
  }
  //
  assert(isphere_nearest >= 0 && isphere_nearest < spheres.size());
  const SphereForRayTracing &obj = spheres[isphere_nearest];
  const delfem2::CVec3d x = r.o + r.d * dist_nearest; // hitting point
  const delfem2::CVec3d n = (x - obj.p).normalized();  // normal at hitting point
  return {true, x, n, obj.c, obj.e, obj.refl};
}

template<class RENDER_TARGET>
delfem2::CVec3d Radiance(
    const Ray &r,
    int depth,
    std::array<unsigned short, 3> &Xi,
    const RENDER_TARGET &render_target,
    int max_depth) {
  using VEC = delfem2::CVec3d;
  auto[is_intersect, hit_pos, hit_nrm, mtrl_refl, mtrl_emis, mtrl_type]
  = RayIntersection<delfem2::CVec3d>(r, render_target);

  if (!is_intersect) { return {0, 0, 0}; } // no intersection

  VEC nl = hit_nrm.dot(r.d) < 0 ? hit_nrm : hit_nrm * -1.;  // normal should have component of ray direction
  depth++;
  if (depth > max_depth) { return mtrl_emis; }
  if (depth > 5) { // Russian roulette
    double p = (mtrl_refl.x > mtrl_refl.y && mtrl_refl.x) >
        mtrl_refl.z ? mtrl_refl.x : mtrl_refl.y > mtrl_refl.z ? mtrl_refl.y : mtrl_refl.z;  // maximum reflectance
    if (delfem2::MyERand48<double>(Xi) < p) {
      mtrl_refl = mtrl_refl * (1 / p);
    } else {
      return mtrl_emis;
    }
  }
  if (mtrl_type == DIFF) {
    const VEC d = SampleHemisphereNormalCos(nl, delfem2::RandomVec2<double>(Xi));
    return mtrl_emis + mtrl_refl.mult(Radiance(Ray(hit_pos, d), depth, Xi, render_target, max_depth));
  } else if (mtrl_type == SPEC) {
    const VEC r0 = Radiance(Ray(hit_pos, r.d - hit_nrm * 2. * hit_nrm.dot(r.d)), depth, Xi, render_target, max_depth);
    return mtrl_emis + mtrl_refl.mult(r0);
  }
  assert(mtrl_type == REFR);
  Ray reflRay(hit_pos, r.d - hit_nrm * 2. * hit_nrm.dot(r.d)); // reflection
  bool into = hit_nrm.dot(nl) > 0;
  double nc = 1;
  double nt = 1.5;
  double nnt = into ? nc / nt : nt / nc;
  double ddn = r.d.dot(nl);
  double cos2t;
  if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) {
    return mtrl_emis + mtrl_refl.mult(Radiance(reflRay, depth, Xi, render_target, max_depth));
  }
  VEC tdir = (r.d * nnt - hit_nrm * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalized();
  double a = nt - nc;
  double b = nt + nc;
  double R0 = a * a / (b * b);
  double c = 1 - (into ? -ddn : tdir.dot(hit_nrm));
  double Re = R0 + (1 - R0) * c * c * c * c * c;
  double Tr = 1 - Re;
  double P = .25 + .5 * Re;
  double RP = Re / P;
  double TP = Tr / (1 - P);
  if (depth > 2) {
    if (delfem2::MyERand48<double>(Xi) < P) {
      const VEC r0 = Radiance(reflRay, depth, Xi, render_target, max_depth) * RP;
      return mtrl_emis + mtrl_refl.mult(r0);
    } else {
      Ray ray0{hit_pos, tdir};
      VEC r0 = Radiance(ray0, depth, Xi, render_target, max_depth) * TP;
      return mtrl_emis + mtrl_refl.mult(r0);
    }
  } else {
    VEC r1 = Radiance(reflRay, depth, Xi, render_target, max_depth) * Re;
    Ray ray0{hit_pos, tdir};
    VEC r2 = Radiance(ray0, depth, Xi, render_target, max_depth) * Tr;
    return mtrl_emis + mtrl_refl.mult(r1 + r2);
  }
}

int main() {
  using namespace delfem2;
  std::vector<SphereForRayTracing> aSphere = {//Scene: radius, position, emission, color, material
      SphereForRayTracing(1e5, CVec3d(1e5 + 1, 40.8, 81.6), CVec3d(), CVec3d(.75, .25, .25), DIFF),//Left
      SphereForRayTracing(1e5, CVec3d(-1e5 + 99, 40.8, 81.6), CVec3d(), CVec3d(.25, .25, .75), DIFF),//Rght
      SphereForRayTracing(1e5, CVec3d(50, 40.8, 1e5), CVec3d(), CVec3d(.75, .75, .75), DIFF),//Back
      SphereForRayTracing(1e5, CVec3d(50, 40.8, -1e5 + 170), CVec3d(), CVec3d(), DIFF),//Frnt
      SphereForRayTracing(1e5, CVec3d(50, 1e5, 81.6), CVec3d(), CVec3d(.75, .75, .75), DIFF),//Botm
      SphereForRayTracing(1e5, CVec3d(50, -1e5 + 81.6, 81.6), CVec3d(), CVec3d(.75, .75, .75), DIFF),//Top
      SphereForRayTracing(16.5, CVec3d(27, 16.5, 47), CVec3d(), CVec3d(1, 1, 1) * .999, SPEC),//Mirr
      SphereForRayTracing(16.5, CVec3d(73, 16.5, 78), CVec3d(), CVec3d(1, 1, 1) * .999, REFR),//Glas
      SphereForRayTracing(600, CVec3d(50, 681.6 - .27, 81.6), CVec3d(12, 12, 12), CVec3d(), DIFF) //Lite
  };

  delfem2::opengl::CTexRGB_Rect2D tex;
  {
    tex.width = 400;
    tex.height = 400;
    tex.channels = 3;
    tex.pixel_color.resize(tex.width * tex.height * tex.channels);
  }
  delfem2::glfw::CViewer3 viewer(2);
  viewer.width = 300;
  viewer.height = 300;
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  tex.InitGL();

  const unsigned int nh = tex.height;
  const unsigned int nw = tex.width;
  Ray cam(
      CVec3d(50, 52, 295.6),
      CVec3d(0, -0.042612, -1).normalized()); // cam pos, dir
  CVec3d cx = CVec3d(nw * .5135 / nh, 0, 0);
  CVec3d cy = (cx.cross(cam.d)).normalized() * .5135;
  std::vector<float> afRGB(tex.height * tex.width * 3, 0.f);

  unsigned int isample = 0;
  auto render = [&](int iw, int ih) {
    std::array<unsigned short, 3> Xi = {
        (unsigned short) ih,
        (unsigned short) iw,
        (unsigned short) (isample * isample)}; // random seed
    CVec3f r_ave(0,0,0);
    for (int sy = 0; sy < 2; sy++) {
      for (int sx = 0; sx < 2; sx++) {
        const auto dx = SampleTent<double>(Xi);
        const auto dy = SampleTent<double>(Xi);
        CVec3d d = cx * (((sx + .5 + dx) / 2 + iw) / nw - .5) +
            cy * (((sy + .5 + dy) / 2 + ih) / nh - .5) + cam.d;
        CVec3d r = Radiance(Ray(cam.o + d * 140., d.normalized()), 0, Xi, aSphere, 10);
        r_ave += r.cast<float>() * 0.25;
      }
    }
    float *ptr = afRGB.data() + (ih * nw + iw) * 3;
    const float isamplef = static_cast<float>(isample);
    ptr[0] = (isamplef * ptr[0] + r_ave[0]) / (isamplef + 1.f);
    ptr[1] = (isamplef * ptr[1] + r_ave[1]) / (isamplef + 1.f);
    ptr[2] = (isamplef * ptr[2] + r_ave[2]) / (isamplef + 1.f);
  };

  while (!glfwWindowShouldClose(viewer.window)) {
    parallel_for(nw, nh, render);
    isample++;
    for (unsigned int ih = 0; ih < tex.height; ++ih) {
      for (unsigned int iw = 0; iw < tex.width; ++iw) {
        for (int ic = 0; ic < 3; ++ic) {
          float fc = afRGB[(ih * tex.width + iw) * 3 + ic];
          fc = (fc > 1.f) ? 1.f : fc;
          fc = (fc < 0.f) ? 0.f : fc;
          int ifc = int(pow(fc, 1 / 2.2) * 255 + .5);
          tex.pixel_color[(ih * tex.width + iw) * 3 + ic] = ifc;
        }
      }
    }
    std::cout << " " << isample << std::endl;
    tex.InitGL();
    //
    viewer.DrawBegin_oldGL();
    ::glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    ::glClear(GL_COLOR_BUFFER_BIT);
    ::glDisable(GL_LIGHTING);
    ::glMatrixMode(GL_PROJECTION);
    ::glLoadIdentity();
    ::glMatrixMode(GL_MODELVIEW);
    ::glLoadIdentity();
    tex.Draw_oldGL();
    viewer.SwapBuffers();
    glfwPollEvents();
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();

  exit(EXIT_SUCCESS);
}
