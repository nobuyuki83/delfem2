
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
    delfem2::CVec3d op = p - r.o;
    double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
    if (det < 0) {
      return 0;
    } else {
      det = sqrt(det);
    }
    return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
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
  double t;
  int id;
  auto nsphere = static_cast<double>(spheres.size());
  double inf = t = 1e20;
  for (int i = int(nsphere); i--;) {
    const double d = spheres[i].intersect(r);
    if ((d > 0) && d < t) {
      t = d;
      id = i;
    }
  }
  if (t > inf) {
    return {false, VEC(), VEC(), VEC(), VEC(), DIFF};
  }
  //
  const SphereForRayTracing &obj = spheres[id];
  const delfem2::CVec3d x = r.o + r.d * t; // hitting point
  const delfem2::CVec3d n = (x - obj.p).normalized();  // normal at hitting point
  return {true, x, n, obj.c, obj.e, obj.refl};
}

template<class RENDER_TARGET>
delfem2::CVec3d Radiance(
    const Ray &r,
    int depth,
    std::array<unsigned short, 3> &Xi,
    const RENDER_TARGET &render_target) {
  using VEC = delfem2::CVec3d;
  auto[is_intersect,
       his_pos, hit_nrm,
       mtrl_refl, mtrl_emis, mtrl_type] = RayIntersection<delfem2::CVec3d>(r, render_target);
  if (!is_intersect) {
    return {0, 0, 0};
  }
  VEC nl = hit_nrm.dot(r.d) < 0 ? hit_nrm : hit_nrm * -1.;  // normal should have component of ray direction
  if (++depth > 5) { // Russian roulette
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
    return mtrl_emis + mtrl_refl.mult(Radiance(Ray(his_pos, d), depth, Xi, render_target));
  } else if (mtrl_type == SPEC) {
    const VEC r0 = Radiance(Ray(his_pos, r.d - hit_nrm * 2. * hit_nrm.dot(r.d)), depth, Xi, render_target);
    return mtrl_emis + mtrl_refl.mult(r0);
  }
  Ray reflRay(his_pos, r.d - hit_nrm * 2. * hit_nrm.dot(r.d));
  bool into = hit_nrm.dot(nl) > 0;
  double nc = 1;
  double nt = 1.5;
  double nnt = into ? nc / nt : nt / nc;
  double ddn = r.d.dot(nl);
  double cos2t;
  if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) {
    return mtrl_emis + mtrl_refl.mult(Radiance(reflRay, depth, Xi, render_target));
  }
  VEC tdir = (r.d * nnt - hit_nrm * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalized();
  double a = nt - nc;
  double b = nt + nc;
  double R0 = a * a / (b * b);
  double c = 1 - (into ? -ddn : tdir.dot(hit_nrm));
  double Re = R0 + (1 - R0) * c * c * c * c * c;
  double Tr = 1 - Re;
  double P = .25 + .5 * Re, RP = Re / P;
  double TP = Tr / (1 - P);
  if (depth > 2) {
    if (delfem2::MyERand48<double>(Xi) < P) {
      return mtrl_emis + mtrl_refl.mult(
          Radiance(reflRay, depth, Xi, render_target) * RP);
    } else {
      return mtrl_emis + mtrl_refl.mult(Radiance(Ray(his_pos, tdir), depth, Xi, render_target) * TP);
    }
  } else {
    return mtrl_emis + mtrl_refl.mult(
        Radiance(reflRay, depth, Xi, render_target) * Re + Radiance(Ray(his_pos, tdir), depth, Xi, render_target) * Tr);
  }
}

inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toIntGammaCorrection(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

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
  viewer.width = 400;
  viewer.height = 400;
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
  while (!glfwWindowShouldClose(viewer.window)) {
    for (unsigned int ih = 0; ih < nh; ih++) {
      std::array<unsigned short, 3> Xi = {0, 0, (unsigned short) (ih * ih * ih + isample * isample)}; // random seed
      for (unsigned int iw = 0; iw < nw; iw++) {
        for (int sy = 0; sy < 2; sy++) {
          for (int sx = 0; sx < 2; sx++) {
            const auto dx = SampleTent<double>(Xi);
            const auto dy = SampleTent<double>(Xi);
            CVec3d d = cx * (((sx + .5 + dx) / 2 + iw) / nw - .5) +
                cy * (((sy + .5 + dy) / 2 + ih) / nh - .5) + cam.d;
            CVec3d r = Radiance(Ray(cam.o + d * 140., d.normalized()), 0, Xi, aSphere);
            afRGB[(ih * nw + iw) * 3 + 0] += static_cast<float>(r.x);
            afRGB[(ih * nw + iw) * 3 + 1] += static_cast<float>(r.y);
            afRGB[(ih * nw + iw) * 3 + 2] += static_cast<float>(r.z);
          }
        }
      }
    }
    isample++;
    for (unsigned int ih = 0; ih < tex.height; ++ih) {
      for (unsigned int iw = 0; iw < tex.width; ++iw) {
        for (int ic = 0; ic < 3; ++ic) {
          float fc = afRGB[(ih * tex.width + iw) * 3 + ic] * 0.25f / float(isample);
          fc = (fc > 1.f) ? 1.f : fc;
          tex.pixel_color[(ih * tex.width + iw) * 3 + ic] = toIntGammaCorrection(fc);
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
