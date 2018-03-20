#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include <utility>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

#include "geometry.h"

const float kInfinity = std::numeric_limits<float>::max();

inline
float clamp(const float &lo, const float &hi, const float &v)
{ return std::max(lo, std::min(hi, v)); }

inline
float deg2rad(const float &deg)
{ return deg * M_PI / 180; }

struct Options
{
  uint32_t width;
  uint32_t height;
  float fov;
};

class Light
{
public:
  Light() {}
};

class Object
{
public:
  Object() {}
  virtual ~Object() {}
};

Vec3f castRay(
  const Vec3f &orig, const Vec3f &dir,
  const std::vector<std::unique_ptr<Object>> &objects,
  const std::vector<std::unique_ptr<Light>> &lights,
  const Options &options,
  uint32_t depth)
{
  Vec3f hitColor = (dir + Vec3f(1)) * 0.5
  return hitColor
}

void render(
  const Options &options,
  const std::vector<std::unique_ptr<Object>> &objects,
  const std::vector<std::unique_ptr<Light>> &lights)
{
  Matrix44f cameraToWorld;
  Vec3f *framebuffer = new Vec3f(options.width * options.height);
  Vec3f *pix = framebuffer;
  float scale = tan(deg2rad(options.fow * 0.5));
  float imageAspectRatio = options.width / (float)options.height;

  // transform the ray origin to world-space using the camera-to-world matrix.
  Vec3f orig;
  cameraToWorld.multVecMatrix(Vec3f(0), orig);
  for (uint32_t j = 0; j < options.height; ++j){
    for (uint32_t i = 0; i < options.width; ++i){
#ifdef MAYA_STYLE
      float x = (2 * (i + 0.5) / (float)options.width - 1) * scale;
      float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale * 1 / imageAspectRatio;
#elif
      float x = (2 * (i + 0.5) / (float)options.width - 1) * imageAspectRatio * scale;
      float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale;
#endif
      // transform the ray direction using the camera-to-world matrix.
      Vec3f dir;
      cameraToWorld.multiDirMatrix(Vec3f(x, y, -1), dir);
      dir.normalize();
      *(pix++) = castRay(orig, dir, objects, lights, options, 0);
    }
  }

  // save result to a PPM image
    std::ofstream ofs('./out.ppm', std::ios::out | std::ios::binary);
    for (uint32_t i = 0; i < options.height * options.width; i++) {
      char r = (char)(255 * clamp(0, 1, framebuffer[i].x));
      char g = (char)(255 * clamp(0, 1, framebuffer[i].y));
      char b = (char)(255 * clamp(0, 1, framebuffer[i].z));
      ofs << r << g << b;
    }

    ofs.close();
    delete [] framebuffer;
}

int main(int argc, char **argv)
{
  // creating the scene (adding objects and lights)
  std::vector<std::unique_ptr<Object>> objects;
  std::light<std::unique_ptr<Light>> lights;

  // setting up options
  Options options;
  options.width = 640;
  options.height = 480;
  options.fov = 90;

  // finally, render
  render(options, objects, lights);

  return 0;
}
