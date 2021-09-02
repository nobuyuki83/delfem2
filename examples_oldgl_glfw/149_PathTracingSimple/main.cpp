
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

#include "delfem2/opengl/tex.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/vec3.h"

using namespace delfem2;

// based on https://github.com/postgres/postgres/blob/master/src/port/erand48.c
double my_erand48(unsigned short xseed[3])
{
  constexpr unsigned short my_rand48_mult[3] = {
      0xe66d,
      0xdeec,
      0x0005
  };
  constexpr unsigned short my_rand48_add = 0x000b;

  unsigned long accu;
  unsigned short temp[2];

  accu = (unsigned long) my_rand48_mult[0] * (unsigned long) xseed[0] +
         (unsigned long) my_rand48_add;
  temp[0] = (unsigned short) accu;	/* lower 16 bits */
  accu >>= sizeof(unsigned short) * 8;
  accu += (unsigned long) my_rand48_mult[0] * (unsigned long) xseed[1] +
          (unsigned long) my_rand48_mult[1] * (unsigned long) xseed[0];
  temp[1] = (unsigned short) accu;	/* middle 16 bits */
  accu >>= sizeof(unsigned short) * 8;
  accu += my_rand48_mult[0] * xseed[2] + my_rand48_mult[1] * xseed[1] + my_rand48_mult[2] * xseed[0];
  xseed[0] = temp[0];
  xseed[1] = temp[1];
  xseed[2] = (unsigned short) accu;
  // --------
  return ldexp((double) xseed[0], -48) +
         ldexp((double) xseed[1], -32) +
         ldexp((double) xseed[2], -16);
}

struct Ray {
  CVec3d o, d;
  Ray(const CVec3d& o_, const CVec3d& d_) : o(o_), d(d_) {}
};

enum Refl_t {
  DIFF, SPEC, REFR
};

class CSphere {
public:
  CSphere(
      double rad_,
      const CVec3d& p_,
      const CVec3d& e_,
      const CVec3d& c_,
      Refl_t refl_)
      : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
  double intersect(const Ray &r) const {
    CVec3d op = p - r.o;
    double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
    if (det < 0) {
      return 0;
    } else{
      det = sqrt(det);
    }
    return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
  }
public:
  double rad;
  CVec3d p, e, c;
  Refl_t refl;
};

inline bool intersect(
    const Ray &r,
    double &t,
    int &id,
    const std::vector<CSphere>& spheres)
{
  double n = spheres.size();
  double inf = t = 1e20;
  for (int i = int(n); i--;) {
    const double d = spheres[i].intersect(r);
    if ((d>0) && d < t) {
      t = d;
      id = i;
    }
  }
  return t < inf;
}

CVec3d radiance(
    const Ray &r,
    int depth,
    unsigned short *Xi,
    const std::vector<CSphere>& aSphere)
{
  double t;
  int id = 0;
  if (!intersect(r, t, id, aSphere)){
    return CVec3d();
  }
  const CSphere &obj = aSphere[id];
  CVec3d x = r.o + r.d * t;
  CVec3d n = (x - obj.p).normalized();
  CVec3d nl = n.dot(r.d) < 0 ? n : n * -1., f = obj.c;
  double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;
  if (++depth > 5){
    if (my_erand48(Xi) < p){
      f = f * (1 / p);
    } else{
      return obj.e;
    }
  }
  if (obj.refl == DIFF) {
    double r1 = 2 * M_PI * my_erand48(Xi);
    double r2 = my_erand48(Xi);
    double r2s = sqrt(r2);
    CVec3d w = nl;
    CVec3d u = ((fabs(w.x) > .1 ? CVec3d(0, 1, 0) : CVec3d(1,0,0)) ^ w).normalized(), v = w ^ u;
    CVec3d d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalized();
    return obj.e + f.mult(radiance(Ray(x, d), depth, Xi, aSphere));
  }
  else if (obj.refl == SPEC) {
    return obj.e + f.mult(radiance(Ray(x, r.d - n * 2. * n.dot(r.d)), depth, Xi, aSphere));
  }
  Ray reflRay(x, r.d - n * 2. * n.dot(r.d));
  bool into = n.dot(nl) > 0;
  double nc = 1;
  double nt = 1.5;
  double nnt = into ? nc / nt : nt / nc;
  double ddn = r.d.dot(nl);
  double cos2t;
  if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) {
    return obj.e + f.mult(radiance(reflRay, depth, Xi, aSphere));
  }
  CVec3d tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalized();
  double a = nt - nc;
  double b = nt + nc;
  double R0 = a * a / (b * b);
  double c = 1 - (into ? -ddn : tdir.dot(n));
  double Re = R0 + (1 - R0) * c * c * c * c * c;
  double Tr = 1 - Re;
  double P = .25 + .5 * Re, RP = Re / P;
  double TP = Tr / (1 - P);
  if( depth > 2 ){
    if( my_erand48(Xi) < P ){
      return obj.e + f.mult(
          radiance(reflRay, depth, Xi, aSphere) * RP );
    }
    else {
      return obj.e + f.mult(radiance(Ray(x, tdir), depth, Xi, aSphere) * TP);
    }
  } else {
    return obj.e + f.mult(radiance(reflRay, depth, Xi, aSphere) * Re + radiance(Ray(x, tdir), depth, Xi, aSphere) * Tr);
  }
}

int main(int argc, char *argv[]) {
  std::vector<CSphere> aSphere = {//Scene: radius, position, emission, color, material
      CSphere(1e5, CVec3d(1e5 + 1, 40.8, 81.6), CVec3d(), CVec3d(.75, .25, .25), DIFF),//Left
      CSphere(1e5, CVec3d(-1e5 + 99, 40.8, 81.6), CVec3d(), CVec3d(.25, .25, .75), DIFF),//Rght
      CSphere(1e5, CVec3d(50, 40.8, 1e5), CVec3d(), CVec3d(.75, .75, .75), DIFF),//Back
      CSphere(1e5, CVec3d(50, 40.8, -1e5 + 170), CVec3d(), CVec3d(), DIFF),//Frnt
      CSphere(1e5, CVec3d(50, 1e5, 81.6), CVec3d(), CVec3d(.75, .75, .75), DIFF),//Botm
      CSphere(1e5, CVec3d(50, -1e5 + 81.6, 81.6), CVec3d(), CVec3d(.75, .75, .75), DIFF),//Top
      CSphere(16.5, CVec3d(27, 16.5, 47), CVec3d(), CVec3d(1, 1, 1) * .999, SPEC),//Mirr
      CSphere(16.5, CVec3d(73, 16.5, 78), CVec3d(), CVec3d(1, 1, 1) * .999, REFR),//Glas
      CSphere(600, CVec3d(50, 681.6 - .27, 81.6), CVec3d(12, 12, 12), CVec3d(), DIFF) //Lite
  };

  delfem2::opengl::CTexRGB_Rect2D tex;
  {
    tex.width = 256;
    tex.height = 256;
    tex.aRGB.resize(tex.width*tex.height*3);
  }
  delfem2::glfw::CViewer3 viewer;
  viewer.width = 400;
  viewer.height = 400;
  viewer.camera.view_height = 2;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;

  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  tex.InitGL();

  const unsigned int nh = tex.height;
  const unsigned int nw = tex.width;
  Ray cam(
      CVec3d(50, 52, 295.6),
      CVec3d(0, -0.042612, -1).normalized()); // cam pos, dir
  CVec3d cx = CVec3d(nw * .5135 / nh, 0, 0);
  CVec3d cy = (cx ^ cam.d).normalized() * .5135;
  std::vector<float> afRGB(tex.height*tex.width*3, 0.f);
  unsigned int isample = 0.0;
  while (!glfwWindowShouldClose(viewer.window))
  {
    for (unsigned int ih = 0; ih < nh; ih++) {
      unsigned short Xi[3] = {0, 0, (unsigned short)(ih * ih * ih + isample * isample)}; // random seed
      for (unsigned int iw = 0; iw < nw; iw++){
        for (int sy = 0; sy < 2; sy++) {
          for (int sx = 0; sx < 2; sx++) {
            const double r1 = 2 * my_erand48(Xi);
            const double r2 = 2 * my_erand48(Xi);
            const double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
            const double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
            CVec3d d = cx * (((sx + .5 + dx) / 2 + iw) / nw - .5) +
                       cy * (((sy + .5 + dy) / 2 + ih) / nh - .5) + cam.d;
            CVec3d r = radiance(Ray(cam.o + d * 140., d.normalized()), 0, Xi, aSphere);
            afRGB[(ih*nw+iw)*3+0] += r.x;
            afRGB[(ih*nw+iw)*3+1] += r.y;
            afRGB[(ih*nw+iw)*3+2] += r.z;
          }
        }
      }
    }
    isample++;
    for(unsigned int ih=0;ih<tex.height;++ih){
      for(unsigned int iw=0;iw<tex.width;++iw) {
        for(int ic=0;ic<3;++ic) {
          float fc = afRGB[(ih * tex.width + iw) * 3 + ic]*0.25f/float(isample);
          fc = (fc>1.f) ? 1.f:fc;
          tex.aRGB[(ih * tex.width + iw) * 3 + ic] = 255*fc;
        }
      }
    }
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
